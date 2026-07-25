# PHB-GTDB-GPT Agent Notes

This file is the current repo-level memory for future sessions.
Keep it aligned with the code, docs, and runtime state.

## Current state

- Project: PHB/PHA depolymerase-related gene analysis on GTDB Release 11 / R232.
- Working repo: `D:\PHB-GTDB-GPT`
- GitHub remote: `https://github.com/anzhiw2-gif/PHB-GTDB-GPT`
- Execution host: T141 at `/home/data/haoyu/PHB-GTDB-GPT`

## Stage status

- P01: complete
- P02: complete
- P03: complete
- P04: next

## Verified outputs

- P01 raw GTDB copy and reconciliation are complete.
- P02 benchmark locked the production predictor to `Pyrodigal GeneFinder(meta=True)`.
- P03 full proteome prediction is complete.
- P03 outputs:
  - `03_gtdb_proteomes/faa/`
  - `03_gtdb_proteomes/gff/`
  - `03_gtdb_proteomes/qc/p03_prediction_qc.tsv`
  - `03_gtdb_proteomes/manifests/p03_prediction_manifest.tsv`
- P03 run summary on T141:
  - genomes predicted: `199,923`
  - predicted genes: `615,969,593`

## Key files

- `docs/DATA_PROVENANCE.md`
- `docs/P02_BENCHMARK_DECISION.md`
- `docs/PREDICTION_POLICY.md`
- `scripts/p01_audit_gtdb.py`
- `scripts/p02_select_benchmark_genomes.py`
- `scripts/p02_compare_predictors.py`
- `scripts/p03_predict_proteomes.py`
- `scripts/p03_monitor_progress.py`

## Working rules

- Keep Chinese as the main language for commit messages, PR text, and high-level summaries.
- Keep English technical terms when precision matters: `GTDB`, `Pyrodigal`, `HMMER`, `InterProScan`, `Nextflow`, `COBRApy`, and similar.
- Do not commit raw GTDB data, full predicted proteomes, HMMER outputs, InterPro outputs, or Nextflow `work/` directories.
- Update docs when inputs, thresholds, reference sets, or interpretation rules change.
- Use the existing stage names `P01` through `P10` only.
- Prefer small, testable changes and keep provenance explicit.

## Useful commands

```powershell
python -m unittest tests/test_repository_layout.py -v
python scripts/validate_repository.py
git diff --check
```

On T141:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_progress.py
```

Update this file whenever the stage status or the project-wide operating rules change.
