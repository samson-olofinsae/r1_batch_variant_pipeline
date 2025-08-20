# R1 Batch Variant Calling Pipeline (SNV/INDEL Splitter)

A lightweight, teaching-friendly batch pipeline that performs end-to-end variant calling from paired-end FASTQ files, then **separates SNVs and INDELs** into dedicated VCFs. Designed to be **fault-tolerant** (skips samples with missing mates) and **transparent** (timestamped logs for auditability).

---

## Why this exists
- **Realistic training**: demonstrates alignment → variant calling → post-processing in a few, readable steps.
- **Reproducibility**: consistent folder structure, deterministic commands, and log files per run.
- **Fault tolerance**: gracefully skips incomplete pairs instead of crashing entire batches.

---

## Features
- Batch processing of all `*_R1.fastq.gz` files in `inputs/` with automatic `_R2` pairing.
- Error-tolerant: logs **“Missing R2… skipping.”** and continues.
- Clear audit trail: a global `run_batch.log` with start/finish timestamps plus per-sample stdout/stderr logs.
- Clean output tree: BAM → BCF → VCF → split **SNVs** and **INDELs**.

---

## Requirements
- **Linux/macOS shell** (Bash)
- **Python 3.8+**
- **Tools**: [`bwa`](http://bio-bwa.sourceforge.net/), [`samtools`](http://www.htslib.org/), [`bcftools`](http://www.htslib.org/), `bgzip` (part of htslib)
- **Reference genome** (example path used below): `ref/hg19_chr8.fa`
  - Make sure the reference is indexed for `bwa` (`.bwt` etc.) and `samtools` (`.fai`).

> **Note on ploidy**
> The example command uses `bcftools call --ploidy 1` for simplicity. For diploid human autosomes you may prefer `--ploidy 2`. Adjust to your study design.

---

### Download the reference FASTA

GitHub does not track files >100 MB. Fetch the small demo FASTA used in examples:

```bash
bash scripts/download_reference.sh
```
This saves `ref/hg19_chr8.fa` (with `.fai` and BWA index files) under `ref/`.

> **Required for this pipeline:** Alignment and variant calling need a reference FASTA. If you already have a reference, update the path in the wrapper:
> ```bash
> # scripts/run_batch_pipeline.sh
> REF_GENOME="ref/hg19_chr8.fa"   # change this to your own reference path if needed
> ```
> To use a different reference (including a stable Dropbox link with `?dl=1`), pass a URL and basename:
> ```bash
> bash scripts/download_reference.sh ref <YOUR_URL> <BASENAME>
> ```

---

## Repository layout
```
project_root/
├── inputs/                   # Paired FASTQs: sample_R1.fastq.gz + sample_R2.fastq.gz
├── ref/
│   └── hg19_chr8.fa         # Reference FASTA (+ indexes)
├── results/ (auto-generated at runtime)
│   ├── bam/
│   ├── bcf/
│   ├── logs/                # run_batch.log + per-sample .log/.err
│   └── vcf/
│       ├── indels/
│       └── snvs/
└── scripts/
    ├── run_batch_pipeline.sh
    └── split_snv_indel_variants.py
```

**Note:** `results/` and all its subfolders are created automatically by the pipeline at runtime. You only need to prepare `inputs/` and `ref/` ahead of time.

---

## Getting the reference (without committing large files)

GitHub does **not** track large files. Use the helper script to fetch and index a small demo reference (hg19 chr8 by default). This is **required** for alignment/variant calling:

```bash
bash scripts/download_reference.sh
```
This will create:

- `ref/hg19_chr8.fa`
- `ref/hg19_chr8.fa.fai` and `ref/hg19_chr8.fa.*` (BWA indexes)

To use a different reference, pass a URL and basename:

```bash
bash scripts/download_reference.sh ref <YOUR_URL> <BASENAME>
```

### Keep references out of Git (one-time cleanup)

If you accidentally committed reference files, remove them from version control and ignore them going forward:

```bash
git rm --cached -r ref
printf "results/\ninputs/*\n!inputs/README.md\nref/*\n!ref/README.md\n*.log\n*.err\n" >> .gitignore
mkdir -p inputs ref
echo "# inputs/\nPlace paired FASTQs here (sample_R1.fastq.gz, sample_R2.fastq.gz)." > inputs/README.md
echo "# ref/\nReference FASTA and indexes live here; fetched via scripts/download_reference.sh." > ref/README.md
git add .gitignore inputs/README.md ref/README.md
git commit -m "Stop tracking large files; add ignore rules and placeholders"
git push
```

---

## Quick start
1. **Clone/copy** this repository structure.
2. Place paired FASTQs in `inputs/` with names like `sample1_R1.fastq.gz` and `sample1_R2.fastq.gz`; `sample2_R1.fastq.gz` and `sample2_R2.fastq.gz`; `sample3_R1.fastq.gz` and `sample3_R2.fastq.gz`.
3. Put your reference FASTA (and indexes) in `ref/` (or update the path in the wrapper script).
4. Run:

```bash
bash scripts/run_batch_pipeline.sh
```

The wrapper discovers `*_R1.fastq.gz` files, checks each has a matching `_R2`, and then calls the Python pipeline script per sample.

---

## What the pipeline does (per sample)
1. **Align & sort**: `bwa mem` → `samtools sort` → `*.aligned.sorted.bam`
2. **Index BAM**: `samtools index`
3. **Pileup**: `bcftools mpileup` → `*_raw.bcf`
4. **Call variants**: `bcftools call -m -v` → `*_variants.vcf` → `bgzip` + index
5. **Split by type**: `bcftools view -v snps` → `snvs/*.vcf`; `bcftools view -v indels` → `indels/*.vcf`

Commands are issued from `scripts/split_snv_indel_variants.py` in the following sequence:
```
bwa mem {ref} {R1} {R2} | samtools sort -o results/bam/{sample}.aligned.sorted.bam
samtools index -b results/bam/{sample}.aligned.sorted.bam
bcftools mpileup -O b -o results/bcf/{sample}_raw.bcf -f {ref} results/bam/{sample}.aligned.sorted.bam
bcftools call --ploidy 1 -m -v -o results/vcf/{sample}_variants.vcf results/bcf/{sample}_raw.bcf
bgzip -c results/vcf/{sample}_variants.vcf > results/vcf/{sample}_final_variants.vcf.gz
bcftools index results/vcf/{sample}_final_variants.vcf.gz
bcftools view -v snps   results/vcf/{sample}_final_variants.vcf.gz > results/vcf/snvs/{sample}_snvs.vcf
bcftools view -v indels results/vcf/{sample}_final_variants.vcf.gz > results/vcf/indels/{sample}_indels.vcf
```

---

## Running in batch mode
The wrapper script `scripts/run_batch_pipeline.sh` controls batching, logging, and error handling.

### Example log (with a missing mate)
This demonstrates what happens if one sample (e.g., `sample1`) has no matching `_R2` file.
```
Batch processing started: Wed Aug 20 18:38:15 BST 2025
Missing R2 pair for sample1, skipping.
Processing sample2...
Completed sample2
Processing sample3...
Completed sample3
Batch run finished: Wed Aug 20 18:50:09 BST 2025
```

### Default log (all samples run successfully)
If all three pairs are present (`sample1`, `sample2`, and `sample3`), the batch log would look like this:
```
Batch processing started: Wed Aug 20 18:38:15 BST 2025
Processing sample1...
Completed sample1
Processing sample2...
Completed sample2
Processing sample3...
Completed sample3
Batch run finished: Wed Aug 20 18:50:09 BST 2025
```

Per-sample logs are written to `results/logs/{sample}.log` with corresponding errors (if any) in `{sample}.err`.

---

## Usage details
- **Input naming**: samples must follow `NAME_R1.fastq.gz` / `NAME_R2.fastq.gz`.
- **Reference path**: edit `REF_GENOME="ref/hg19_chr8.fa"` in the wrapper if your reference lives elsewhere.
- **Exit behaviour**: the wrapper uses `set -e` and checks return codes; on a per-sample failure, it prints a helpful message and continues to the next sample.

---

## Output artefacts (per sample)
- `results/bam/{sample}.aligned.sorted.bam` (+ `.bai`)
- `results/bcf/{sample}_raw.bcf`
- `results/vcf/{sample}_final_variants.vcf.gz` (+ `.tbi`)
- `results/vcf/snvs/{sample}_snvs.vcf`
- `results/vcf/indels/{sample}_indels.vcf`

---

## Troubleshooting
- **“Command not found”**: ensure `bwa`, `samtools`, `bcftools`, and `bgzip` are installed and on your `PATH`.
- **“.fai not found”**: create `samtools faidx ref.fa` and `bwa index ref.fa`.
- **Ploidy/organism**: adjust `--ploidy` to your organism and sample type.
- **Performance**: for speed on large datasets, consider multithreading `bwa mem -t N` and `samtools sort -@ N`.

---

## Reproducibility & teaching notes
- This pipeline emphasises **clarity over cleverness**: minimal flags, explicit files, readable logs.
- The fault-tolerant batch runner is ideal for demonstrating **resilience patterns** (skip bad inputs, keep the run alive, inspect logs).

---

## Roadmap (nice-to-haves)
- Configurable ploidy/parameters via CLI flags.
- Optional container recipes (Docker/Singularity) and a Makefile or Nextflow wrapper.
- Basic unit tests on tiny synthetic FASTQs.

---

## License
MIT - © 2025 Samson Olofinsae. See the [LICENSE](./LICENSE) file for details.

---

## Acknowledgements
Authored by **Samson Olofinsae**. Built as part of hands-on training in robust, reproducible genomics workflows.
