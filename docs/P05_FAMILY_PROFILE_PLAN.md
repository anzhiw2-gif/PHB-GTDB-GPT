# P05 Family Profile Plan

**Snapshot date:** 2026-07-26

P05 turns the curated P04 seed library into a family-by-family profile plan.
The project boundary stays explicit: sequence evidence can justify a family
assignment or homology signal, but it does not by itself prove PHB/PHA
degradation phenotype.

## Current Snapshot

The current normalized P04 manifest has 42 seed rows across 10 families after
the current seed-extension pass, including fresh extracellular scl type I / type
II reinforcement and a five-sequence intracellular mcl-PHA reinforcement, but
the default P05 planner now reads
[`04_family_profiles/manifests/p05_family_keep_now.tsv`](../04_family_profiles/manifests/p05_family_keep_now.tsv)
and keeps only the 6 review-level main branches active.

Under the current rule, five active bacterial families and the archaeal branch
now reach the threshold of three independent qualifying source accessions
needed to justify a custom HMM. The active six-family subset contains 38
qualifying seed rows after the intracellular mcl-PHA family was raised from 3
to 8 independent UniProtKB accessions.

For the keep/defer split and its literature basis, see
[P05 family classification priorities](P05_FAMILY_CLASSIFICATION_PRIORITIES.md).

The detailed 2026-07-27 seed review is recorded in
[P05 HMM seed selection decision](P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md).
That review supersedes the earlier assumption that all eight current archaeal
rows are suitable members of one patatin-like alignment. It also records the
rule that bacterial profile seeds require experimental support, while
accessioned archaeal E3 rows may be used for explicitly labelled coverage.

The GitHub-trackable model metadata, current 38-seed registry, proposed seed
updates, local-artifact checksums, and P06 approval gate are recorded in
[P05 HMM model catalog](P05_HMM_MODEL_CATALOG.md). The raw HMM, MAFFT
alignment, and seed-bundle artifacts remain machine-local; the catalog records
their SHA256 values instead of committing the generated files.

That means the first P05 output is a planning queue:

- HMM-ready families
- anchor-set families
- a compact summary of the current evidence balance

## Threshold Rule

- `build_hmm` only when a family has at least three independent qualifying
  source accessions
- `anchor_set` when the family is below threshold
- `E1` and `E2` rows count as qualifying evidence for bacteria
- archaeal E3 annotation-supported rows can also count as qualifying evidence
- `Excluded` rows stay out of profile building

The qualifying-row threshold does not override architecture coherence. Rows
that pass the accession/evidence threshold but fail the family architecture
check are boundary candidates and must not be used to inflate a custom HMM.

## Planned Outputs

The planner writes these ignored TSVs under `04_family_profiles/manifests/`:

- `p05_family_profile_plan.tsv`
- `p05_family_hmm_build_queue.tsv`
- `p05_family_anchor_set_queue.tsv`
- `p05_family_profile_summary.tsv`

## Build Scaffold

The next deterministic step consumes both the normalized P04 manifest and the
saved P05 plan. If the saved plan is missing, the scaffold first materializes
it from the current manifest and keep-now family classification, then rechecks
that the manifest and plan agree on the qualifying source accessions for each
family before writing a build queue for the eligible families.

The scaffold writes these ignored TSVs under `04_family_profiles/manifests/`:

- `p05_family_hmm_build_scaffold_queue.tsv`
- `p05_family_hmm_build_scaffold_summary.tsv`

For each eligible family, the scaffold also writes one ignored, unaligned seed
FASTA bundle under `04_family_profiles/seed_bundles/` by default. The queue
records the exact `seed_bundle_path` and `bundled_sequence_count`. Bundle
records are sorted by `seed_id`, `source_accession`, and `sequence_path`; each
header is the stable `seed_id|source_accession` identifier, and sequence lines
use a fixed 80-character width. The input sequence paths are resolved from the
P04 manifest location or, when needed, from the repository root, so the bundle
is a direct, auditable materialization of the saved reference library rather
than a new evidence source.

The scaffold only prepares these inputs. It does not run MAFFT, `hmmbuild`, or
`hmmsearch`. A custom HMM remains gated by the three-independent-accession
rule and must be calibrated against close non-target hydrolases before use.

The scaffold records the planned MAFFT alignment mode and HMMER build/search
tools in the queue so later stages can run them without re-deriving the
decision.

The command-manifest step consumes the saved scaffold queue and writes these
ignored TSVs under `04_family_profiles/manifests/`:

- `p05_family_profile_command_manifest.tsv`
- `p05_family_profile_command_summary.tsv`

Each eligible family gets one deterministic record containing the seed bundle,
future alignment and HMM paths, a MAFFT L-INS-i command, and an `hmmbuild`
command with the HMM output path before the alignment input path. The status is
`planned_not_run`; this step never executes either tool. The current real queue
now has five bacterial families plus one archaeal branch, so the command
manifest contains six command records and the summary reports six prepared
commands. The summary records `eligible_families` alongside the scaffold queue
and command row counts.
Serialized command paths are written in slash-only form so the TSV stays
portable across Windows and POSIX hosts.

The command forms follow the MAFFT manual's L-INS-i definition
([MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html))
and the HMMER documentation page for `hmmbuild`, where the HMM output file
precedes the multiple-alignment input file
([HMMER documentation](https://hmmer.org/documentation.html)).

Run it from the repo root with:

```powershell
python scripts/p05_plan_family_profiles.py --build-scaffold
```

Prepare the command manifest from the saved scaffold queue with:

```powershell
python scripts/p05_family_profile_commands.py
```

Use `--bundle-dir PATH` to place the ignored seed bundles in an explicit
directory, for example in a temporary test workspace.

## Biological Rationale

Custom HMMs can overfit badly when the seed set is too small. A three-seed
minimum keeps the model from collapsing onto one quirky sequence and makes it
possible to calibrate against close non-target hydrolases instead of treating a
single accession as a family definition.
