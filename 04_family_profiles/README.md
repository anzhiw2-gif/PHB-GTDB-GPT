# P05 Family Profiles

Family seed sets, alignments, calibrated HMMs, and model-quality records will be generated here.

Current P05 planner and build scaffold live in
[scripts/p05_plan_family_profiles.py](../scripts/p05_plan_family_profiles.py).

Current rule boundary:

- build custom HMMs only when a family has at least three independent qualifying source accessions
- plan active P05 work from `p05_family_keep_now.tsv`, not from every historical seed-family label
- keep smaller families as explicit anchor sets
- for the archaeal branch, annotation-supported E3 rows can count toward the qualifying accession total when they are explicitly admitted in the seed manifest
- calibrate any future HMM against close non-target hydrolases before it is used downstream

Planner outputs remain under `04_family_profiles/manifests/`:

- `p05_family_keep_now.tsv`
- `p05_family_deferred.tsv`
- `p05_family_profile_plan.tsv`
- `p05_family_hmm_build_queue.tsv`
- `p05_family_anchor_set_queue.tsv`
- `p05_family_profile_summary.tsv`

The build scaffold adds:

- `p05_family_hmm_build_scaffold_queue.tsv`
- `p05_family_hmm_build_scaffold_summary.tsv`
- one ignored unaligned seed bundle per eligible family under
  `04_family_profiles/seed_bundles/`

The scaffold queue records `seed_bundle_path` and `bundled_sequence_count`.
Bundle headers are stable `seed_id|source_accession` identifiers and records
are written in deterministic seed order. These FASTA files are inputs for
future MAFFT/HMMER work only; the scaffold does not run those tools.

If the saved P05 plan is missing, the scaffold first materializes it from the
current manifest and keep-now family classification, then writes the build
queue from that saved plan.

The command-manifest step adds these ignored files:

- `p05_family_profile_command_manifest.tsv`
- `p05_family_profile_command_summary.tsv`

For a literature-backed keep/defer split of the current seed families, see
[P05 family classification priorities](../docs/P05_FAMILY_CLASSIFICATION_PRIORITIES.md).

The manifest contains one deterministic MAFFT L-INS-i and `hmmbuild` command
record per eligible family. Its `planned_not_run` status is explicit, and HMM
calibration against close non-target hydrolases remains required before
downstream use. With the current six-family active P05 snapshot, five bacterial
families plus one archaeal branch are now command-ready. The command summary
also reports `eligible_families`, and the serialized paths use slash-only form
so the manifest remains portable between Windows and POSIX hosts.

When an eligible family is present, the command-manifest step also prepares the
ignored `04_family_profiles/alignments/` and `04_family_profiles/hmms/`
directories for the future MAFFT and HMMER outputs.

As of 2026-07-26, the normalized P04 reference library contains 31 seed rows
across 10 families after the latest cross-genus bacterial expansion. P05 keeps
6 active main branches, and all 6 now meet the three-seed HMM threshold, so
the real-world build scaffold queue now contains 6 eligible families.
