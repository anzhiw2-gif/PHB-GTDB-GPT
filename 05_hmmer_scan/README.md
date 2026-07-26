# P06 HMMER Scan

P06 turns the curated P05 HMMs into raw `domtblout` scan outputs over the P03
GTDB proteome shards, then converts those raw hits into a candidate catalog.
Sequence hits are still sequence evidence only; they do not prove phenotype.

As of 2026-07-26, the GTDB-wide P06 scan is pending. The previous P03 FAA
outputs had empty protein translations, so the full scan must wait until the
P03 translation-fix rerun finishes and passes QC.

## Expected Outputs

Under `05_hmmer_scan/`, the scan stage writes:

- `p06_hmmer_scan_manifest.tsv`
- `p06_hmmer_scan_summary.tsv`
- `p06_hmmer_candidates.tsv`
- `p06_hmmer_candidate_summary.tsv`
- `raw_domtblout/` for the preserved HMMER `domtblout` files
- `hmmer_logs/` for the main `hmmsearch` output streams

## Repository Entry Point

Plan the scan jobs from the current HMM and proteome directories:

```powershell
python scripts/p06_scan_family_profiles.py
```

On T141, the P03 proteomes are nested `.faa.gz` files. Use chunked streaming
instead of one process per genome:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
export PATH=/home/data/haoyu/miniconda3/envs/phb_gtdb/bin:$PATH
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p06_scan_family_profiles.py \
  --proteomes-per-job 200
```

Parse an existing manifest after the raw `domtblout` files have been created:

```powershell
python scripts/p06_scan_family_profiles.py --parse-only
```

The scanner uses `hmmsearch --domtblout` and keeps the raw output separate
from derived candidate tables so later review can trace every row back to the
exact family HMM and GTDB shard.
