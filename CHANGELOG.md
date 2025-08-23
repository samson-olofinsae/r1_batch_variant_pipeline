# Changelog

## v0.1.0 — Initial public release
- End-to-end R1 → BAM → BCF → VCF variant calling.
- Split variants into **SNVs** and **INDELs** (separate VCFs).
- Fault-tolerant batch wrapper: skips missing `_R2` mates and keeps running.
- Timestamped logs: global `run_batch.log` + per-sample `.log/.err`.
- **MultiQC custom content**: writes `<sample>_r1_stats_mqc.tsv` and can auto-render HTML.
- Clear README with repo layout and quick start; helper to fetch demo reference FASTA.
