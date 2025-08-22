#!/usr/bin/env python3

# Author: Samson Olofinsae
# Purpose: End-to-end variant calling from FASTQ, separating SNVs and indels
# Usage: Called per-sample by scripts/run_batch_pipeline.sh

import sys
import os
import subprocess
from subprocess import call

if len(sys.argv) != 4:
    print("Usage: python split_snv_indel_variants.py <R1.fastq.gz> <R2.fastq.gz> <ref_genome.fa>")
    sys.exit(1)

# Input arguments
fq1 = sys.argv[1]
fq2 = sys.argv[2]
ref_genome_file = sys.argv[3]

# Derive base sample name
sample_base = os.path.basename(fq1).split('_R1.fastq.gz')[0]

# Create working directory references
wd = os.getcwd()
results_dir = os.path.join(wd, "results")
bam_dir = os.path.join(results_dir, "bam")
bcf_dir = os.path.join(results_dir, "bcf")
vcf_dir = os.path.join(results_dir, "vcf")
snv_dir = os.path.join(vcf_dir, "snvs")
indel_dir = os.path.join(vcf_dir, "indels")

# Create folders if not already there
for folder in [bam_dir, bcf_dir, snv_dir, indel_dir]:
    os.makedirs(folder, exist_ok=True)

# Define all output file names
sorted_bam = os.path.join(bam_dir, f"{sample_base}.aligned.sorted.bam")
raw_bcf = os.path.join(bcf_dir, f"{sample_base}_raw.bcf")
variants = os.path.join(vcf_dir, f"{sample_base}_variants.vcf")
final_vcf_gz = os.path.join(vcf_dir, f"{sample_base}_final_variants.vcf.gz")
snv_file = os.path.join(snv_dir, f"{sample_base}_snvs.vcf")
indel_file = os.path.join(indel_dir, f"{sample_base}_indels.vcf")

# Run alignment, variant calling, and separation
commands = [
    f"bwa mem {ref_genome_file} {fq1} {fq2} | samtools sort -o {sorted_bam}",
    f"samtools index -b {sorted_bam}",
    f"bcftools mpileup -O b -o {raw_bcf} -f {ref_genome_file} {sorted_bam}",
    f"bcftools call --ploidy 1 -m -v -o {variants} {raw_bcf}",
    f"bgzip -c {variants} > {final_vcf_gz}",
    f"bcftools index {final_vcf_gz}",
    f"bcftools view -v snps {final_vcf_gz} > {snv_file}",
    f"bcftools view -v indels {final_vcf_gz} > {indel_file}"
]

print(f" Starting variant pipeline for: {sample_base}")
for cmd in commands:
    print(f" Running: {cmd}")
    ret = call(cmd, shell=True)
    if ret != 0:
        print(f" Command failed: {cmd}")
        sys.exit(1)

# ---- MultiQC Custom Content (NEW) ------------------------------------------
# Emit a small TSV with header comments so MultiQC auto-detects & renders it.
mqc_dir = os.path.join(results_dir, "multiqc_cc")
os.makedirs(mqc_dir, exist_ok=True)
mqc_tsv = os.path.join(mqc_dir, f"{sample_base}_r1_stats_mqc.tsv")

def count_records(vcf_path):
    """Count non-header records in a VCF (0 if missing/unreadable)."""
    try:
        out = subprocess.check_output(f"grep -vc '^#' {vcf_path} || echo 0",
                                      shell=True, stderr=subprocess.DEVNULL)
        return int(out.decode().strip())
    except Exception:
        return 0

snv_count = count_records(snv_file)
indel_count = count_records(indel_file)

with open(mqc_tsv, "w") as fh:
    fh.write("# id: r1_variant_splitter\n")
    fh.write("# section_name: R1 Variant Splitter — variant summary\n")
    fh.write("# description: SNV/INDEL counts per sample (from bcftools outputs)\n")
    fh.write("# plot_type: table\n")
    fh.write("# file_format: tsv\n")
    fh.write("Sample\tsnvs\tindels\n")
    fh.write(f"{sample_base}\t{snv_count}\t{indel_count}\n")

print(f" Pipeline completed successfully for: {sample_base}")
