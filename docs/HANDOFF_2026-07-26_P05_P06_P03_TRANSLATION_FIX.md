# Handoff 2026-07-26: P05/P06 and P03 Translation Fix

This handoff records the current repository state after the P05 seed
reinforcement, the first P06 scanner scaffold, and the P03 translation-fix
rerun on T141.

> Historical status only: this 2026-07-26 handoff predates the approved
> four-model P06 execution. For the current P06 outcome, use
> `docs/P06_COMPLETION_2026-07-28.md`.

Claim boundary: HMM, domain, localization, and tree evidence are sequence
evidence only. They do not by themselves prove PHB/PHA degradation phenotype.

## Source Of Truth

- Local repository: `D:\PHB-GTDB-GPT`
- T141 repository: `/home/data/haoyu/PHB-GTDB-GPT`
- GitHub remote: `https://github.com/anzhiw2-gif/PHB-GTDB-GPT`
- T141 Python: `/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python`
- Snapshot date: 2026-07-26

## Stage Status

| Stage | Current status | Evidence and next action |
| --- | --- | --- |
| P01 | Complete | GTDB R232 raw copy and reconciliation were completed on T141. Raw genomes remain ignored by Git and are referenced through manifests and compact docs. |
| P02 | Complete | The production predictor is fixed to `Pyrodigal GeneFinder(meta=True)` based on the deterministic benchmark. |
| P03 | Complete | The original full run reached `199,923` FAA/GFF pairs and `615,969,593` predicted genes, but the FAA records contained empty translations. The translation fix reran P03 on T141, regenerated QC and manifest files, and completed with nonzero FAA residues. |
| P04 | Complete, with P05 seed reinforcement | The normalized reference library now has `42` seed rows: `34` Bacteria and `8` Archaea, with evidence levels `E1=20`, `E2=15`, `E3=7`, `Excluded=0`. |
| P05 | In progress | Six active family branches have enough independent qualifying seed accessions and generated HMM artifacts under ignored machine-local output directories. Final model use remains tied to calibration and provenance checks. |
| P06 | Scaffold ready, full scan blocked | `scripts/p06_scan_family_profiles.py` can plan HMMER jobs and parse existing `domtblout` outputs. Earlier planning output is invalid for a new scan: P06 waits for checksum-locked, calibration-approved models in the registry. |

## P03 Translation Fix

Problem detected:

- The original P03 outputs had the expected file count, but protein FASTA
  records contained headers with empty sequences.
- The concrete implementation issue was that Pyrodigal gene objects can expose
  translation through a `translate()` method rather than a plain translation
  attribute.

Repository fix:

- `scripts/p03_predict_proteomes.py` now calls `gene.translate(include_stop=False)`
  when available.
- Empty protein translations now raise an explicit `ValueError` instead of
  silently writing unusable FAA output.
- `tests/test_p03_predict_proteomes.py` covers the Pyrodigal-like
  `translate(include_stop=False)` path.

T141 monitor:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_translation_fix.py --once
```

Observed T141 snapshot at 2026-07-27 16:30:45 Asia/Shanghai:

- fix-run pidfile:
  `03_gtdb_proteomes/run_logs/p03_translation_fix_20260726.pid`
- process: stopped, pid `2873210`
- rewritten FAA files since fix start: `199,923 / 199,923`
- old-by-mtime files: `0`
- stable rewritten sample had nonzero residue counts
- stderr tail was empty

Regenerated P03 QC/manifest snapshot:

- `03_gtdb_proteomes/qc/p03_prediction_qc.tsv`
- `03_gtdb_proteomes/manifests/p03_prediction_manifest.tsv`
- QC rows: `199,923`
- manifest rows: `199,923`

Acceptance gate before P06 full scan:

- P03 rerun completed without errors.
- `p03_prediction_qc.tsv` and `p03_prediction_manifest.tsv` were regenerated.
- A representative sample of final FAA files has nonzero residue counts.
- No empty-translation failure is present in the P03 run logs.

## P04/P05 Seed Library

Current tracked summary:

| Family | Qualifying accessions |
| --- | ---: |
| `archaeal_patatin_like_pha_dep` | 8 |
| `extracellular_mcl_pha_dep` | 6 |
| `extracellular_scl_pha_dep_type_I` | 6 |
| `extracellular_scl_pha_dep_type_II` | 5 |
| `intracellular_mcl_pha_dep` | 8 |
| `intracellular_phaZ_no_lipase_box` | 5 |

The active six-family subset has `38` qualifying seed rows. The broader
reference library preserves 10 family labels, including deferred or anchor
branches, so later work can recover narrower families without mixing them into
the main P05 models.

The latest P05 reinforcement added these bacterial UniProtKB seed FASTA files:

- `A6EXA3`
- `O82950`
- `Q5YEW3`
- `B7UCC9`
- `Q8VV57`
- `Q9AGB5`
- `Q9Z3Y0`
- `Q5Q135`

The corresponding provenance is recorded in:

- `01_reference_library/manifests/reference_library.seed_manifest.tsv`
- `01_reference_library/reference_library.normalized.tsv`
- `01_reference_library/reference_library.bacteria.normalized.tsv`
- `01_reference_library/reference_library_summary.tsv`
- `01_reference_library/retrieval_logs/p05_seed_retrieval_log.tsv`

## P05 HMM Outputs

Generated P05 outputs are intentionally machine-local and ignored by Git:

- `04_family_profiles/seed_bundles/`
- `04_family_profiles/alignments/`
- `04_family_profiles/hmms/`
- `04_family_profiles/calibration/`

On the current workspace, six HMM files are present under
`04_family_profiles/hmms/`, one per active family. They are not committed
because generated HMM artifacts and calibration output are treated as derived
analysis outputs. Git keeps the scripts, README files, small seed sequences,
and compact provenance tables needed to regenerate and audit them.

## P06 Scanner Scaffold

Tracked P06 files:

- `scripts/p06_scan_family_profiles.py`
- `tests/test_p06_scan_family_profiles.py`
- `05_hmmer_scan/README.md`

The scaffold creates deterministic `hmmsearch --domtblout` command manifests,
supports recursive `.faa.gz` P03 proteome discovery, and can group proteomes
into streamed chunks with `--proteomes-per-job`.

Recommended T141 command after P03 repair acceptance:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
export PATH=/home/data/haoyu/miniconda3/envs/phb_gtdb/bin:$PATH
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p06_scan_family_profiles.py \
  --proteomes-per-job 200
```

The P06 planning output written before P05 model approval is not a valid new
scan manifest. Rebuild it only after every selected model has
`approved_for_p06=yes`, `scan_permission=approved`, and a matching HMM SHA256.

P06 output directories and bulk tables are ignored by Git. Commit only compact
accepted summaries or curated candidate tables after review.

## Verification

Commands run locally on 2026-07-26:

```powershell
python -m unittest tests/test_repository_layout.py tests/test_p03_predict_proteomes.py tests/test_p06_scan_family_profiles.py -v
python scripts/validate_repository.py
```

Result:

- 17 tests passed.
- Repository layout validation passed.

Useful final pre-push checks:

```powershell
git diff --check
git status --short --branch
```

## Open Items

- P03 translation-fix rerun completed on T141; keep the regenerated QC and
  manifest under review if downstream assumptions change.
- P06 GTDB-wide HMMER scan remains blocked until the P05 model registry records
  calibrated approval. Raw domtblout and candidate tables stay ignored until
  compact summaries are accepted.
- P05 generated HMMs exist as ignored machine-local artifacts; final model
  calibration/provenance should be checked before downstream biological
  interpretation.
