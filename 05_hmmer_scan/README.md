# P06 HMMER Scan

P06 turns the curated P05 HMMs into raw `domtblout` scan outputs over the P03
GTDB proteome shards, then converts those raw hits into a candidate catalog.
Sequence hits are still sequence evidence only; they do not prove phenotype.

The P03 translation-fix rerun has completed and passed its nonempty-protein
QC. P06 remains frozen because a new scan must use only user-confirmed,
checksum-locked models. See [P05 HMM model catalog](../docs/P05_HMM_MODEL_CATALOG.md).

## Expected Outputs

Under `05_hmmer_scan/`, the scan stage writes:

- `p06_hmmer_scan_manifest.tsv`
- `p06_hmmer_scan_summary.tsv`
- `p06_hmmer_candidates.tsv`
- `p06_hmmer_candidate_summary.tsv`
- `raw_domtblout/` for the preserved HMMER `domtblout` files
- `hmmer_logs/` for the main `hmmsearch` output streams

## Repository Entry Point

Plan the scan jobs from the P06-approved HMM rows in the model registry:

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
exact family HMM and GTDB shard. It refuses to plan a scan unless the registry
marks a row `approved_for_p06=yes` and `scan_permission=approved`, and the
local HMM SHA256 matches the registered value. The scan manifest carries that
model SHA256 for every job.
