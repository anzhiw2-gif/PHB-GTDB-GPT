# P07 Domain Annotation

This stage prepares independent domain-architecture and localization evidence
for P06 sequence candidates. It is a review stage: InterProScan, SignalP, and
optional Phobius evidence can support candidate interpretation, but cannot by
itself prove PHB/PHA degradation phenotype.

## Current entry point

First generate the compact P06 reasonableness audit on T141:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p06_candidate_reasonableness.py \
  --candidate-table 05_hmmer_scan/p06_hmmer_candidates.tsv \
  --outdir 06_domain_annotation/p06_reasonableness \
  --total-predicted-genes 615969593 \
  --total-genomes 199923
```

Use the P07 preparation script on T141 after confirming the P06 candidate table
and scan manifest are present in the isolated P06 worktree:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p07_prepare_domain_annotation.py \
  --candidate-table 05_hmmer_scan/p06_hmmer_candidates.tsv \
  --scan-manifest 05_hmmer_scan/p06_hmmer_scan_manifest.tsv \
  --outdir 06_domain_annotation \
  --include-tier High-confidence \
  --sequences-per-shard 500
```

The script writes FASTA shards under `06_domain_annotation/input/fasta_shards/`
and planned annotation commands under
`06_domain_annotation/manifests/p07_domain_annotation_command_manifest.tsv`.
Generated candidate FASTA, raw InterPro/localization output, review tables, and
P07 command manifests remain machine-local unless a compact summary is later
accepted for Git tracking.

Before running annotation commands, preflight the command manifest:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p07_run_domain_annotation.py \
  --manifest 06_domain_annotation/manifests/p07_domain_annotation_command_manifest.tsv \
  --status-dir 06_domain_annotation/run_status \
  --workers 1 \
  --preflight-only
```

The initial 2026-07-28 preflight found `missing_executable=152` in the checked
environment, but that startup issue was resolved on T141 and the full P07
batch completed on 2026-07-29. The final run-status table reported
`completed=103`, `skipped_existing=49`, `failed_exit_code=0`, and
`failed_missing_output=0`.

## Evidence boundary

- InterProScan evidence is domain/family/site architecture evidence.
- SignalP6 evidence is signal-peptide/localization support.
- Optional Phobius evidence is signal-peptide and transmembrane-topology
  support after the local CLI has been verified.
- None of these outputs are phenotype proof; they are filters for P08/P09
  review and reporting.
