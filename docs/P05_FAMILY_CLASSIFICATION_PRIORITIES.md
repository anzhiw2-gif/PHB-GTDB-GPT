# P05 Family Classification Priorities

**Snapshot date:** 2026-07-26

This note records a literature-backed simplification of the current P04/P05
seed-family set. The goal is not to delete evidence, but to keep only the
major review-level branches active for the next round of curation and push the
rest into an explicit deferred list.

## What the recent reviews emphasize

Recent PHA/PHB depolymerase reviews keep returning to the same two main axes:

- localization: intracellular, extracellular, and periplasmic depolymerases
- substrate preference: short-chain-length (scl) versus medium-chain-length
  (mcl) PHA

That means the broadest useful operational classes are:

- extracellular scl-PHA depolymerases
- extracellular mcl-PHA depolymerases
- intracellular PHA depolymerases
- periplasmic depolymerases
- archaeal patatin-like PHA depolymerases as a separate archaeal branch

The more specific historical family names remain biologically real, but they
are better treated as refinement branches unless we have enough independent
seed evidence to justify a separate HMM.

## Why this matches the GTDB/HMM side

GTDB and GTDB-Tk rely on curated HMM marker sets and model-specific cutoffs,
not on maximal seed inflation. That is the same design principle we want here:
keep the family set small enough to stay interpretable and calibratable, then
expand only when a branch has enough independent experimental support.

## Keep now

These are the families I would keep active for the next pass:

- `intracellular_phaZ_no_lipase_box`
- `intracellular_mcl_pha_dep`
- `extracellular_scl_pha_dep_type_I`
- `extracellular_scl_pha_dep_type_II`
- `extracellular_mcl_pha_dep`
- `archaeal_patatin_like_pha_dep`

Reason in one line: they cover the main review-level trunks, preserve the
subcellular and substrate splits that matter most, and retain the archaeal
branch without mixing it into bacterial HMM space.

## Defer for later

These families are still valid evidence records, but I would postpone them for
now:

- `phaZ7_like`
- `phaZd_like`
- `rhodospirillum_periplasmic_like`
- `tigr02240_aromatic_pha_related`

Reason in one line: they are narrower branches, boundary anchors, or
historically named subtypes that are useful later, but they are not needed to
represent the major review-level structure of the current seed library.

## Decision rule

- keep only the major trunks in the near-term curation set
- preserve the deferred families in explicit files rather than dropping them
- do not promote a deferred family into a custom HMM until it has enough
  independent qualifying seed evidence

## P05 planner integration

The default P05 planner consumes
`04_family_profiles/manifests/p05_family_keep_now.tsv` through
`--family-classification`. Rows with `priority_status=keep_now` are retained
for the active family-profile plan; deferred rows remain recorded separately
and are not included in the normal P05 build scaffold.

## Source links

- GTDB methods: https://gtdb.ecogenomic.org/methods
- GTDB-Tk classify workflow:
  https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html
- 2024 review on PHA synthase and depolymerase proteins:
  https://eprints.whiterose.ac.uk/223456/1/s10924-024-03474-4.pdf
- 2025 Frontiers review on PHA biodegradation:
  https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2025.1542468/full
