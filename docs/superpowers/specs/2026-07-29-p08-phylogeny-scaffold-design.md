# P08 Phylogeny and GTDB Taxonomy Scaffold Design

**Date:** 2026-07-29, Asia/Shanghai
**Project stage:** P08
**Status:** approved design; implementation not started

## Goal

Create a reproducible, fail-closed P08 scaffold that turns the completed P06
High-confidence candidate set and completed P07 annotation outputs into
family-level phylogeny/taxonomy execution manifests. The scaffold must prepare
and validate inputs; it must not start MAFFT, IQ-TREE 2, FastTree, or make a
PHB/PHA-degradation phenotype claim.

## Scope and biological boundary

P08 evaluates evolutionary placement and GTDB taxonomic distribution of
protein candidates. P06 HMM matches, P07 domain architecture, P07 signal
peptide/localization predictions, family trees, and GTDB species-tree context
are all sequence or annotation evidence. None of them alone or in combination
proves PHB/PHA degradation phenotype.

The scaffold applies only to the P06 `High-confidence` input tier by default.
P06 `Review` rows require a documented secondary-pass decision. P06 `Rejected`
rows are excluded. A multi-model match remains an overlap review flag; the
scaffold must retain every relevant model assignment rather than forcing a
mutually exclusive family call.

## Inputs and provenance contract

The implementation will require explicit paths for the following machine-local
or tracked inputs:

- P06 candidate table and scan manifest from the accepted four-model R232 scan;
- P07 candidate-sequence manifest, FASTA shards, annotation command manifest,
  and final run-status table from the completed High-confidence pass;
- the tracked P05 HMM model registry and seed registry, which identify the
  approved family models and checksum-locked reference sequences;
- GTDB R232 bacterial and archaeal taxonomy files, plus the corresponding
  species-tree files when available;
- P03 prediction manifest/QC and FAA paths for source-protein traceability.

For every P08 candidate-manifest row, retain: GTDB release; assembly accession;
protein identifier; source FAA path; P06 model, score, coverage, threshold and
tier; P07 output/status paths; complete GTDB lineage; family-task identifier;
and source/input SHA-256 values. The manifest must preserve the original P06
candidate identifier and never substitute an untraceable FASTA header.

## Components and data flow

`p08_prepare_phylogeny.py` will be the deterministic preparation and
validation entry point. It will read the declared input manifests, validate
row-level joins and sequence lengths, calculate family/task sizes, and write:

- a candidate-to-family task manifest, including multi-model overlap flags;
- a taxonomy-join manifest and an explicit unmapped-taxonomy review table;
- one sequence-input manifest per family task, including curated seed and
  close non-target-control provenance;
- a deterministic execution-command manifest and a compact preparation
  summary;
- a blocker/review table for every integrity failure.

`p08_run_phylogeny.py` will consume only the prepared command manifest. Its
first supported use is `--preflight-only`: resolve executable paths, capture
versions, verify input checksums and required output locations, and write a
resumable status table. Actual command execution is a later, separately
approved P08 operation.

The scripts will write distinct per-task paths and atomic status updates. They
will never append concurrent jobs into a shared output file. Generated FASTA,
alignments, trees, tool logs, raw GTDB joins, and large command manifests remain
machine-local under `07_phylogeny/`; Git will retain only source code, tests,
compact tracked summaries, and accepted documentation.

## Family-size routing and tree policy

The preparation manifest will calculate the size of every complete candidate
pool before selecting a command route:

| Complete family-task size | Planned alignment/tree route |
| --- | --- |
| fewer than 200 proteins | MAFFT L-INS-i; review alignment; IQ-TREE 2 only after final subset/outgroup approval |
| 200 to 2,000 proteins | MAFFT `--auto`; review alignment; IQ-TREE 2 only after final subset/outgroup approval |
| more than 2,000 proteins | deterministic dereplication/representative-selection plan plus exploratory FastTree route; no final interpretation without an approved subset and outgroup decision |

The first preparation pass must not silently discard sequences. For the
more-than-2,000 route, it will record the deterministic representative-selection
policy, parameters, input checksum, and excluded-to-representative mapping for
review before use. Reference seeds and accessioned close non-target hydrolases
are explicit controls, not GTDB candidate calls.

Final interpretable trees must use IQ-TREE 2 with model testing and 1,000
ultrafast bootstrap replicates. Rooting requires an explicit accessioned
outgroup; if no justified outgroup exists, only midpoint display rooting may be
used and must be labelled as such. GTDB species trees provide contextual
comparison only; gene-tree/species-tree discordance is not, by itself, evidence
for phenotype or horizontal transfer.

## Fail-closed behavior

A family task must be blocked rather than run when any selected candidate has a
missing P06/P07 join, missing/empty source sequence, P06-to-FASTA length
disagreement, failed checksum verification, absent taxonomy mapping, unknown
model-registry row, missing required executable, or missing expected output.
Blocked records must be written to a review table with a machine-readable
reason. They are data-integrity or execution states, not biological negatives.

P07 annotations are linked as independent evidence; absence or ambiguity in a
localization/domain prediction is recorded as unknown/review rather than
converted to a family exclusion without a documented biological rule.

## Tests and acceptance criteria

Unit tests will use tiny synthetic P06, P07, P05, taxonomy, and FASTA fixtures
to verify deterministic ordering, accession/provenance preservation,
multi-model overlap retention, domain-separated taxonomy joins, size-route
selection, checksum/length failure closure, missing-executable preflight, and
resumable status behavior.

Before any T141 calculation, acceptance requires:

1. all selected P08 rows resolve to a P06 source and P07 sequence record;
2. every resolved FASTA sequence length equals the P06 target length;
3. every family task has a declared model/seed/control provenance record;
4. taxonomy mapping and species-tree availability are separately reported for
   bacterial and archaeal records;
5. executable versions, command manifest checksums, and preflight status are
   captured; and
6. no blocked row is silently emitted as a runnable tree command.

Repository validation after implementation will include the P08-specific unit
tests, `python -m unittest tests/test_repository_layout.py -v`,
`python scripts/validate_repository.py`, and `git diff --check`. Any
compute-heavy T141 preflight or later tree run will record its exact command,
environment path, tool versions, input/output paths, and why it was not run
locally.

## Non-goals

This design does not authorize a P08 full run, P06/P07 rerun, seed-model
recalibration, Review-tier expansion, candidate deletion, profile subtype
classification, phenotype inference, or P09 reporting. Each requires its own
documented decision when reached.
