# P04 Reference Library

P04 creates the auditable reference layer for the PHB/PHA depolymerase review.
It keeps sequence evidence, literature provenance, and family-definition
anchors separate from downstream HMM calibration.

## What P04 Produces

- curated seed sequences with provenance
- a normalized seed manifest
- separate normalized bacterial and archaeal seed manifests
- retrieval logs for the source databases and literature
- separate family-definition assets for DED/TIGRFAM anchors

## Reference Libraries

P04 builds two seed libraries with different admission rules.

- `bacteria_high_confidence`: bacterial PHB/PHA depolymerase seeds with direct
  experimental activity or strong functional evidence with clear limitations.
  Annotation-only `E3` rows are excluded from this library.
- `archaea_literature_supported`: archaeal PHB/PHA depolymerase-related seeds
  that are supported by current literature. These rows must cite a PMID, DOI,
  or PMCID and must state what the literature supports. Automatic database
  annotation alone is not enough.

The split is intentional. Bacterial PHB/PHA depolymerases have multiple
literature-supported families, while archaeal candidates are sparse and should
not be used to inflate bacterial positive seeds.

Within the bacterial library, keep the PHB-specific core strict and allow only
a small PHA bridge set when it is needed to cover major DED family branches
that lack a cleaner PHB-only seed. Bridge rows should stay explicitly labeled
in notes and must not replace the PHB core.

## Canonical Family Categories

- `intracellular_phaZ_no_lipase_box`
- `extracellular_scl_pha_dep_type_I`
- `extracellular_scl_pha_dep_type_II`
- `phaZ7_like`
- `phaZd_like`
- `rhodospirillum_periplasmic_like`
- `intracellular_mcl_pha_dep`
- `extracellular_mcl_pha_dep`
- `tigr02240_aromatic_pha_related`
- `archaeal_patatin_like_pha_dep`
- `auxiliary_mobilization_context`

## Evidence Levels

- `E1`: direct experimental activity
- `E2`: functional evidence with limitations
- `E3`: annotation or homology only
- `Excluded`: tracked but not used as a reference seed

For `bacteria_high_confidence`, use only `E1` or `E2` in non-excluded rows.
For `archaea_literature_supported`, use only `E1` or `E2` in non-excluded rows
and record the literature support scope explicitly. Archaeal proteins inferred
only from InterPro, Pfam, TIGRFAM/NCBIfam, or NCBI/UniProt automatic names
remain candidates for later review, not seed rows.

## Current Archaeal Literature Anchor

The initial archaeal anchor is `PhaZh1/HFX_6464` from
*Haloferax mediterranei* ATCC 33500. The NCBI Protein accession is
`AFK21580.1`, and the UniProtKB accession is `I3RBH0`. PMID `25710370` reports
that a patatin-like protein associated with *H. mediterranei* PHA granules acts
as an efficient depolymerase in native PHA degradation. This supports an
archaeal literature-supported seed, but it should remain separate from the
bacterial PHB/PHA depolymerase seed library.

## Stage Documents

- [P04 retrieval contract](P04_RETRIEVAL_CONTRACT.md)
- [P04 manifest schema](P04_MANIFEST_SCHEMA.md)

## Implementation Entry

The current normalizer is [scripts/p04_build_reference_library.py](../scripts/p04_build_reference_library.py).
