# P06 HMMER Scan

P06 turns the curated P05 HMMs into raw `domtblout` scan outputs over the P03
GTDB proteome shards, then converts those raw hits into a candidate catalog.
Sequence hits are still sequence evidence only; they do not prove phenotype.

The P03 translation-fix rerun has completed and passed its nonempty-protein
QC. P06 may use only the four user-confirmed, checksum-locked models approved
in the P05 registry: `archaeal_patatin_like_pha_dep`,
`intracellular_mcl_pha_dep`, `intracellular_phaZ_no_lipase_box`, and
`extracellular_pha_depolymerase_core`. The extracellular core identifies
sequence candidates only; it must not assign mcl/scl or type-I/type-II labels.
See [P05 HMM model catalog](../docs/P05_HMM_MODEL_CATALOG.md).

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

Run the resulting manifest through the resumable executor. It skips existing
nonempty `domtblout` files, records every job in an atomically rewritten status
TSV, and leaves failed or empty-output jobs available for a later retry:

```bash
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p06_run_family_profiles.py \
  --manifest 05_hmmer_scan/p06_hmmer_scan_manifest.tsv \
  --status-dir 05_hmmer_scan/run_status \
  --workers 4
```

Parse an existing manifest after the raw `domtblout` files have been created:

```powershell
python scripts/p06_scan_family_profiles.py --parse-only
```

Only parse after the run-status table reports no `failed_exit_code` or
`failed_empty_domtblout` records. A missing raw output is reported in the
candidate summary and is not interpreted as biological absence.

HMMER 3.4 cannot process an individual target protein longer than 100,000 aa
with its standard comparison pipeline. The manifest therefore streams each
proteome chunk through `p06_stream_proteomes.py`. Eligible proteins are passed
to `hmmsearch`; every overlong target is retained in the ignored
`overlong_protein_exclusions/<family>/<chunk>.tsv` audit table with its source
proteome path, target identifier, length, and tool-limit reason. These records
are not HMMER negatives and are not evidence of biological absence.

The scanner uses `hmmsearch --domtblout` and keeps the raw output separate
from derived candidate tables so later review can trace every row back to the
exact family HMM and GTDB shard. It refuses to plan a scan unless the registry
marks a row `approved_for_p06=yes` and `scan_permission=approved`, and the
local HMM SHA256 matches the registered value. The scan manifest carries that
model SHA256 for every job.
