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
- P03: complete
- P04: complete
- P05: in progress

## Verified Outputs

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
- P04 outputs:
  - `01_reference_library/manifests/reference_library.seed_manifest.tsv`
  - `01_reference_library/reference_library.normalized.tsv`
  - `01_reference_library/reference_library.bacteria.normalized.tsv`
  - `01_reference_library/reference_library.archaea.normalized.tsv`
  - `01_reference_library/reference_library_summary.tsv`
  - `01_reference_library/retrieval_logs/p04_seed_retrieval_log.tsv`
  - `01_reference_library/retrieval_logs/p05_seed_retrieval_log.tsv`
  - curated seed FASTA files under `01_reference_library/seeds/`
- P04 manifest summary on 2026-07-26:
  - seed rows: `27`
  - bacterial rows: `19`
  - archaeal rows: `8`
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

## Key Files

- `docs/DATA_PROVENANCE.md`
- `docs/P02_BENCHMARK_DECISION.md`
- `docs/PREDICTION_POLICY.md`
- `scripts/p01_audit_gtdb.py`
- `scripts/p02_select_benchmark_genomes.py`
- `scripts/p02_compare_predictors.py`
- `scripts/p03_predict_proteomes.py`
- `scripts/p03_monitor_progress.py`
- `scripts/p05_plan_family_profiles.py`
- `scripts/p05_family_profile_commands.py`
- `docs/P05_FAMILY_PROFILE_PLAN.md`

## Operating Rules

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
```

Update this file whenever the stage status or the project-wide operating rules
change. Treat this file as the live repo-level memory.
