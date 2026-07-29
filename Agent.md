# PHB-GTDB-GPT Agent Notes

This file is the current repo-level memory for future sessions.
Keep it aligned with the code, docs, and runtime state.

## Project Context

`PHB-GTDB-GPT` is a reproducible PHB/PHA depolymerase-related gene analysis of
GTDB Release 11 R232. The workflow combines experimentally supported reference
sequences, DED family definitions, TIGRFAM profiles, custom HMMs, domain
architecture review, localization evidence, phylogenetics, GTDB taxonomy, and
optional COBRApy validation.

Keep the project claim boundary explicit: homology, domain architecture, and
tree placement are sequence evidence only; they do not by themselves prove a
PHB/PHA degradation phenotype.

## Current State

- Project: PHB/PHA depolymerase-related gene analysis on GTDB Release 11 / R232.
- Working repo: `D:\PHB-GTDB-GPT`
- GitHub remote: `https://github.com/anzhiw2-gif/PHB-GTDB-GPT`
- Execution host: T141 at `/home/data/haoyu/PHB-GTDB-GPT`

## Stage Status

- P01: complete
- P02: complete
- P03: translation-fix rerun complete on T141; QC and manifest regenerated
- P04: complete; seed-admission and boundary-control decisions are recorded
- P05: complete; four scan models are checksum-locked and calibration-approved. The three extracellular subtype artifacts remain blocked and are replaced for P06 entry by `extracellular_pha_depolymerase_core`
- P06: complete for the four-model GTDB R232 HMMER scan and candidate parsing in isolated T141 worktree `/home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728` from commit `80991a7`; 4,000/4,000 raw jobs are accepted and parsing found `missing_domtblout=0`. See `docs/P06_COMPLETION_2026-07-28.md`.
- P07: completed on T141; 37,912/37,912 selected candidate sequences were extracted into 76 FASTA shards, and the full InterProScan/SignalP6 batch completed with 6 workers and `-cpu 10`. The final run-status table reported `completed=103`, `skipped_existing=49`, `failed_exit_code=0`, and `failed_missing_output=0`. See `docs/P07_COMPLETION_2026-07-29.md`.
- P06/P07 audit bridge: compact P06 candidate-table reasonableness audit completed on T141; High-confidence unique targets are `37,912`, with `high_confidence_overlap_targets=0`. Generated audit outputs remain machine-local. See `docs/P06_REASONABLENESS_AUDIT_PLAN_2026-07-28.md`.
- P08: local phylogeny/taxonomy scaffold and `--preflight-only` runner implementation are complete. T141 preflight has not run; no MAFFT alignment, FastTree exploration, IQ-TREE inference, rooting, or biological P08 analysis has been executed. See `docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md`.

## Verified Outputs

- P01 raw GTDB copy and reconciliation are complete.
- P02 benchmark locked the production predictor to `Pyrodigal GeneFinder(meta=True)`.
- The original P03 full-proteome run reached the expected file count, but the
  FAA protein records were later found to contain empty translations. Do not
  use the old FAA files for HMMER scanning.
- The P03 translation-fix rerun completed on T141 on 2026-07-27. The fix uses
  `gene.translate(include_stop=False)` when Pyrodigal exposes translation by
  method, rejects empty protein translations, and is monitored by
  `scripts/p03_monitor_translation_fix.py`.
- P03 machine-local output paths:
  - `03_gtdb_proteomes/faa/`
  - `03_gtdb_proteomes/gff/`
  - `03_gtdb_proteomes/qc/p03_prediction_qc.tsv`
  - `03_gtdb_proteomes/manifests/p03_prediction_manifest.tsv`
- P03 run summary on T141:
  - genomes predicted: `199,923`
  - predicted genes: `615,969,593`
  - current acceptance status: complete; `03_gtdb_proteomes/qc/p03_prediction_qc.tsv`
    and `03_gtdb_proteomes/manifests/p03_prediction_manifest.tsv` were
    regenerated on 2026-07-27 and the rewritten FAA sample checks show
    nonzero residues
- P04 outputs:
  - `01_reference_library/manifests/reference_library.seed_manifest.tsv`
  - `01_reference_library/reference_library.normalized.tsv`
  - `01_reference_library/reference_library.bacteria.normalized.tsv`
  - `01_reference_library/reference_library.archaea.normalized.tsv`
  - `01_reference_library/reference_library_summary.tsv`
  - `01_reference_library/retrieval_logs/p04_seed_retrieval_log.tsv`
  - `01_reference_library/retrieval_logs/p05_seed_retrieval_log.tsv`
  - curated seed FASTA files under `01_reference_library/seeds/`
- P04 manifest summary on 2026-07-27:
  - seed rows: `47`
  - bacterial rows: `35`
  - archaeal rows: `12`
  - evidence levels: `E1=20`, `E2=16`, `E3=11`, `Excluded=0`
  - profile admission: `approved=15`, `boundary_candidate=11`; remaining rows are outside the three revised-family decisions
- P05 planner outputs:
  - `04_family_profiles/manifests/p05_family_profile_plan.tsv`
  - `04_family_profiles/manifests/p05_family_hmm_build_queue.tsv`
  - `04_family_profiles/manifests/p05_family_anchor_set_queue.tsv`
  - `04_family_profiles/manifests/p05_family_profile_summary.tsv`
- P05 scaffold outputs:
  - `04_family_profiles/manifests/p05_family_hmm_build_scaffold_queue.tsv`
  - `04_family_profiles/manifests/p05_family_hmm_build_scaffold_summary.tsv`
  - deterministic eligible-family seed bundles under `04_family_profiles/seed_bundles/`
  - five HMM-ready bacterial families plus one archaeal branch
- P05 command-manifest outputs:
  - `04_family_profiles/manifests/p05_family_profile_command_manifest.tsv`
  - `04_family_profiles/manifests/p05_family_profile_command_summary.tsv`
  - `planned_not_run` command records for the six eligible families
- P05 generated outputs are intentionally ignored by Git:
  - `04_family_profiles/seed_bundles/`
  - `04_family_profiles/alignments/`
  - `04_family_profiles/hmms/`
- P05 Git-tracked compact model provenance:
  - `04_family_profiles/manifests/p05_hmm_model_registry.tsv`
  - `04_family_profiles/manifests/p05_hmm_seed_registry.tsv`
  - `04_family_profiles/manifests/p05_hmm_proposed_seed_updates.tsv`
  - `docs/P05_HMM_MODEL_CATALOG.md`
  - `scripts/p05_catalog_hmm_models.py`
- `archaeal_patatin_like_pha_dep` (6 seeds), `intracellular_mcl_pha_dep`
  (4 seeds), and `intracellular_phaZ_no_lipase_box` (5 seeds) were rebuilt on
  T141 on 2026-07-27 using MAFFT `v7.525` L-INS-i and HMMER `3.4` in the
  isolated `/home/data/haoyu/PHB-GTDB-GPT-p05-rebuild-20260727` worktree.
  Their model, bundle, and alignment SHA256 values are in the tracked model
  registry. All six models remain `approved_for_p06=no`,
  `scan_permission=blocked`, and `calibration_status=not_complete`.
- The seed-selection decision is recorded in
  `docs/P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md`.
- P05 calibration controls and the 2026-07-27 full-model T141 smoke are
  recorded in `docs/P05_HMM_CALIBRATION_PROTOCOL.md`. The smoke found zero
  cross-family hits for the archaeal and two intracellular models, but strong
  extracellular cross-family overlap; this is sequence evidence only and is
  not a P06 approval.
- P05 leave-one-out evidence was completed in isolated
  `/home/data/haoyu/PHB-GTDB-GPT-p05-calibration-r3-20260727` on 2026-07-27:
  31/32 held-out positives were recovered. `AAB40611.1` failed recovery for
  extracellular scl type II; extracellular mcl and scl type I retain hard
  cross-family overlap. Archaeal patatin, intracellular mcl, and intracellular
  no-lipase-box meet the mechanical rule but remain human-review-only. Compact
  outputs are `p05_hmm_leave_one_out_positive_results.tsv`,
  `p05_hmm_control_smoke_results.tsv`, and
  `p05_hmm_calibration_decision_summary.tsv`; raw calibration artifacts remain
  ignored. These archived six-row results are superseded by the approved
  four-model registry decision below.
- The three extracellular subtype HMMs cannot support mutually exclusive P06
  classification: mcl retains seven hard cross-family passes, type I retains
  two, and type II misses `AAB40611.1` in leave-one-out. Do not loosen a cutoff
  or remove an experimentally supported seed to force subtype separation.
  Build and calibrate `extracellular_pha_depolymerase_core` from all 17
  existing extracellular bacterial seeds instead. The hard panel must contain
  the 15 non-extracellular active P05 seeds plus five accessioned close
  non-target hydrolases; see `docs/P05_EXTRACELLULAR_CORE_DECISION_2026-07-27.md`.
- `extracellular_pha_depolymerase_core` was built in isolated T141 r7 with
  MAFFT `v7.525` L-INS-i and HMMER `3.4`: `NSEQ=17`, HMM SHA256
  `74c4b69a2d845f0725d0bc348402e6a51ba3c17a9f67f8cabed3b63df6a6e2f4`.
  It recovered `17/17` held-out positives. Its final rule is full score
  `>=159.9` and HMM coverage `>=0.405498`; `0/20` hard controls pass it.
  Approved P06 models are this core plus archaeal patatin, intracellular mcl,
  and intracellular no-lipase-box. P06 must propagate these per-model
  thresholds from the registry into its manifest and High-confidence tier.
- P06 tracked scaffold:
  - `scripts/p06_scan_family_profiles.py`
  - `scripts/p06_run_family_profiles.py`
  - `scripts/p06_stream_proteomes.py`
  - `tests/test_p06_scan_family_profiles.py`
  - `tests/test_p06_run_family_profiles.py`
- P06 launch record:
  - `docs/P06_LAUNCH_2026-07-28.md`
  - `docs/P06_COMPLETION_2026-07-28.md`
  - P06 accepted 4,000 model-chunk jobs (four approved models x 1,000 chunks); raw outputs, candidate tables, checkpoints, and overlong-target audit tables remain machine-local in r8
  - `05_hmmer_scan/README.md`
- P06 raw HMMER outputs and candidate tables stay ignored after compact
  summaries are accepted; retain them machine-local for P07 traceability.
- P06 planning outputs written before the model-approval gate on 2026-07-27 are not valid for a new scan:
  - `05_hmmer_scan/p06_hmmer_scan_manifest.tsv`
  - `05_hmmer_scan/p06_hmmer_scan_summary.tsv`
- P07 scaffold:
  - `scripts/p06_candidate_reasonableness.py`
  - `scripts/p07_prepare_domain_annotation.py`
  - `scripts/p07_run_domain_annotation.py`
  - `tests/test_p06_candidate_reasonableness.py`
  - `tests/test_p07_prepare_domain_annotation.py`
  - `tests/test_p07_run_domain_annotation.py`
  - `docs/P06_REASONABLENESS_AUDIT_PLAN_2026-07-28.md`
  - `docs/P07_DOMAIN_LOCALIZATION_PLAN_2026-07-28.md`
  - `docs/P07_COMPLETION_2026-07-29.md`
  - run the compact P06 reasonableness audit before P07 sequence extraction; generated audit outputs stay machine-local unless an accepted summary is promoted
  - default input is P06 `High-confidence` only; `Review` requires an explicit secondary pass
  - P07 runner supports `--preflight-only` and resumable status tracking; rerun preflight after installing or loading InterProScan/Java and SignalP6, then run a one-shard smoke before any full batch
  - generated FASTA shards, InterPro/localization outputs, review tables, and P07 generated manifests stay machine-local unless compact summaries are later accepted
- P08 scaffold:
  - `scripts/p08_prepare_phylogeny.py` and `scripts/p08_run_phylogeny.py`
  - `tests/test_p08_prepare_phylogeny.py` and `tests/test_p08_run_phylogeny.py`
  - `docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md` and `07_phylogeny/README.md`
  - preparation accepts only approved P05 models, defaults to P06 `High-confidence`, permits `Review` only when explicit, preserves nonexclusive multi-model targets, and fail-closes unless every P07 row matches the explicit GTDB release and normalized P06 candidate-table/scan-manifest sources; it also rejects P03 FAA-source mismatches, P07 duplicate status keys, family mismatch, or same-stem/different-path status association
  - each candidate retains P06 score/coverage thresholds, verified P03/P06/P07 source paths and whole-file SHA-256 contracts, plus actual InterProScan/SignalP6 terminal statuses and validated input/output paths; this is traceability only, not phenotype evidence
  - the pooled extracellular core authority always derives 17 seeds plus a 15+5 hard panel and rejects direct core seed/control rows in ordinary P05 tables; large-family plans declare an unmaterialized representative-selection algorithm/version/parameters/mapping contract and require future separate approval
  - Bac120/Ar53 taxonomy sources are domain-validated and hash-locked; tree paths/SHA-256 are recorded in the machine-local input provenance manifest as preflight-only inputs, without reading tree topology
  - preparation never runs MAFFT, FastTree, or IQ-TREE and writes command rows as `planned_not_run`
  - the runner CLI remains `--manifest --status-dir --workers --preflight-only`; status binds the full command-manifest SHA-256 and nonexecuted tool-resolution provenance. `skipped_existing` during preflight is only a nonempty-artifact integrity/resume state, not a completed biological analysis. T141 preflight and every actual alignment/tree inference remain not run and need separate authorization
  - P08 family FASTA, mappings, manifests, reviews, statuses, logs, alignments, and trees stay machine-local unless a compact, accepted summary is later promoted

## Key Files

- `docs/DATA_PROVENANCE.md`
- `docs/HANDOFF_2026-07-26_P05_P06_P03_TRANSLATION_FIX.md`
- `docs/P02_BENCHMARK_DECISION.md`
- `docs/PREDICTION_POLICY.md`
- `scripts/p01_audit_gtdb.py`
- `scripts/p02_select_benchmark_genomes.py`
- `scripts/p02_compare_predictors.py`
- `scripts/p03_predict_proteomes.py`
- `scripts/p03_monitor_progress.py`
- `scripts/p03_monitor_translation_fix.py`
- `scripts/p05_plan_family_profiles.py`
- `scripts/p05_family_profile_commands.py`
- `scripts/p05_catalog_hmm_models.py`
- `scripts/p06_scan_family_profiles.py`
- `scripts/p06_candidate_reasonableness.py`
- `docs/P05_FAMILY_PROFILE_PLAN.md`
- `docs/P05_HMM_MODEL_CATALOG.md`
- `docs/P05_HMM_CONSTRUCTION_VALIDATION_RECORD_2026-07-27.md`
- `scripts/p07_prepare_domain_annotation.py`
- `docs/P06_REASONABLENESS_AUDIT_PLAN_2026-07-28.md`
- `docs/P07_DOMAIN_LOCALIZATION_PLAN_2026-07-28.md`
- `scripts/p08_prepare_phylogeny.py`
- `scripts/p08_run_phylogeny.py`
- `docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md`

## Operating Rules

### Completed-Work Check

Before every operation, inspect the current stage record, relevant compact
manifest, Git history/status, and live runtime state to determine whether the
same work has already completed. Do not repeat a completed retrieval, model
build, calibration, validation, copy, or scan merely to recreate prior output.
Reuse the verified artifact and record that it was reused. Repeat work only
when the user explicitly requests it, an input/model checksum changed, a
documented integrity check failed, or a new stage contract requires a distinct
output; state the specific reason before running it. Prefer the narrowest
missing check or next uncompleted stage over broad reruns.

### Before Any Operation

1. Check for a suitable local skill before starting the operation.
   - Read the relevant `SKILL.md` fully before acting.
   - Choose the smallest and fastest skill set that directly fits the task.
   - Prefer domain skills for biological work, such as `nextflow`,
     `bioinformatics-workflows`, `biopython`, `phylogenetics`,
     `protein-sequence-similarity-search`, `protein-sequence-msa`,
     `cobrapy`, `statistical-analysis`, or literature/database skills when
     they match the task.
   - Announce the selected skill briefly and explain why it is suitable.

2. Before running code, search online documentation and GitHub for existing
   reproducible implementations.
   - Prefer official documentation, maintained GitHub projects, workflow
     modules, published supplementary code, and versioned releases.
   - Record the source URL, version/tag/commit, license when relevant, and why
     it is appropriate.
   - Use existing reproducible code or documented command patterns when they
     fit the project. Write new code only after checking that no suitable
     existing implementation is available, or when custom glue code is clearly
     simpler and safer.
   - Do not copy unverifiable snippets into the repository without provenance.

3. Inspect local project context before editing.
   - Read nearby scripts, configs, tests, and docs first.
   - Preserve existing stage names (`P01` through `P10`), path conventions, and
     provenance policies.
   - Do not commit raw GTDB genomes, full predicted proteomes, bulk HMMER or
     InterPro outputs, Nextflow work directories, or other large generated
     artifacts.

### Biological Evidence Requirements

Every step and operation must have biological support.

- For each stage or script change, state the biological rationale: what
  biological signal is being measured, filtered, predicted, annotated, or
  compared.
- For every sequence used, provide a clear source:
  accession, database, release/version, organism/taxon, retrieval date, and
  file path or manifest location.
- Prefer literature-supported reference sequences. For archaeal coverage,
  annotation-supported E3 rows are acceptable when the accession, organism,
  and support scope are explicit. Record DOI, PMID, PMCID, or a stable
  publisher/preprint URL when available.
- Bacterial HMM profile seeds require experimental support tied to the exact
  accession: direct biochemical activity, purified-protein assay, or an
  unambiguous knockout/complementation or physiological experiment. Annotation-
  only E3 bacterial records are not profile seeds.
- Archaeal E3 records may be used for family coverage when direct experiments
  are sparse, but their accession, organism, architecture, and exact support
  scope must be explicit. E3 is not phenotype evidence.
- `intracellular_mcl_pha_dep` is expected to be Pseudomonas-enriched as an
  intrinsic biological property of the family. Do not interpret this pattern
  as seed-selection bias or force artificial cross-genus balance.
- HMM seed inclusion also requires architecture coherence. For
  `archaeal_patatin_like_pha_dep`, retain the PhaZh1-like patatin architecture
  and demote conflicting PHB-synthase-like, HHH, PKD, or AxeA-like records to
  boundary candidates.
- For GTDB-derived sequences, record the GTDB release, assembly accession,
  protein identifier, genome path or manifest reference, and any transformation
  applied by the workflow.
- For UniProt, NCBI, PDB, TIGRFAM, Pfam, InterPro, SignalP, Phobius, HMMER,
  MAFFT, IQ-TREE, FastTree, Pyrodigal, or COBRApy inputs, record database/tool
  versions and retrieval or execution dates.
- If a sequence source is ambiguous, stop and resolve provenance before using
  it in analysis.
- Treat legacy data under `/home/data/haoyu/PHB_gtdb` as read-only historical
  evidence. It may inform checks only after independent verification.

### Reproducibility Rules

- Keep commands deterministic when practical: fixed random seeds, pinned tool
  versions, stable input manifests, and explicit output paths.
- Add or update tests for code changes that affect parsing, filtering,
  sequence selection, metrics, reports, or workflow control.
- Validate lightweight repository changes with:

```powershell
python -m unittest tests/test_repository_layout.py -v
python scripts/validate_repository.py
```

- For stage-specific code, also run the narrowest relevant unit tests.
- If compute-heavy validation must run on T141, document the exact command,
  environment, expected outputs, and why it was not run locally.
- P06 planning must read `p05_hmm_model_registry.tsv`; it may select only
  rows with `approved_for_p06=yes` and `scan_permission=approved`, and must
  verify the local HMM SHA256 before writing a scan-manifest row.

### Documentation And Provenance

- Update `docs/DATA_PROVENANCE.md`, stage decision docs, or compact manifests
  when inputs, databases, reference sequences, thresholds, or biological
  interpretation rules change.
- Write enough provenance that another researcher can reproduce the source
  sequence set and understand why each sequence or model was included.
- Keep generated reports compact in Git; large raw and derived files stay on
  T141 and are referenced through checksums or tracked summaries.

### Git And GitHub

- Before committing or pushing, check `git status --short` and avoid reverting
  unrelated user changes.
- Commit messages, push notes, PR titles, and PR descriptions should be written
  mainly in Chinese. Keep English only for precise technical terms, tool names,
  stage names, file names, and common workflow keywords. Do not write them in
  all English.
- Use a Chinese-first summary and keep English technical keywords only when
  they improve precision.
  Examples:
  - `P03: 加速 GTDB proteome gzip output / speed up gzip output`
  - `P04: 固化 reference provenance / lock sequence provenance`
  - `P06: 调整 HMMER family gates / tune family-specific filters`
- PR descriptions should include:
  - 中文概述为主，必要时保留 English technical terms.
  - 生物学依据和证据来源。
  - 可复现性说明：commands, versions, manifests, and tests.
  - 未运行的验证及原因。

### When Unsure

Prefer evidence over speed. If a biological claim, sequence source, tool
choice, or threshold cannot be justified from documentation, database records,
or literature, pause that part of the work and document the uncertainty instead
of silently proceeding.

## Useful Commands

```powershell
python -m unittest tests/test_repository_layout.py -v
python scripts/validate_repository.py
git diff --check
```

On T141:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_progress.py
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_translation_fix.py --once
```

Update this file whenever the stage status or the project-wide operating rules
change. Treat this file as the live repo-level memory.
