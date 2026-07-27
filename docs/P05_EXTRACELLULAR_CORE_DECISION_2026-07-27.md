# P05 Extracellular PHA-Depolymerase Core Decision

**Decision date:** 2026-07-27
**Stage:** P05 model optimization before P06
**Status:** calibration complete and approved for P06

## Why The Three Extracellular Models Cannot Enter P06 As Subtype Classifiers

The archived extracellular mcl, scl type-I, and scl type-II HMMs were tested
against each other's accessioned seeds before any GTDB-wide search. They do
not provide a mutually exclusive sequence boundary:

| Archived model | Leave-one-out recovery | Strictest retained rule | Hard cross-family seeds passing |
| --- | ---: | ---: | ---: |
| `extracellular_mcl_pha_dep` | 6/6 | score >= 35.5; HMM coverage >= 0.222222 | 7 |
| `extracellular_scl_pha_dep_type_I` | 6/6 | score >= 105.7; HMM coverage >= 0.716495 | 2 |
| `extracellular_scl_pha_dep_type_II` | 4/5 | none | not evaluated |

The low mcl cutoff is caused by the experimentally supported but short,
divergent `CAQ46066.1` recovery (score 35.5; coverage 0.222222). Relaxing or
deleting it merely to create subtype separation would contradict the
accession-level experimental seed rule. Type-II also fails leave-one-out for
the experimentally supported 635-aa `AAB40611.1`, whereas the other type-II
references are 491-495 aa. These are data-supported architecture boundaries,
not errors that can be solved by lowering a threshold.

Therefore, a raw HMM hit cannot be used to assign mcl versus scl substrate
preference, or type I versus type II, in a GTDB protein. Such a claim would
exceed sequence evidence.

## Approved Optimization

The P06 entry model is a single top-level
`extracellular_pha_depolymerase_core` HMM. It pools the current, experimentally
supported bacterial seeds from all three extracellular branches:

- mcl: `BBB62377.1`, `BBB68334.1`, `CAQ46066.1`, `Q51718`, `Q6UFW4`, `Q84C08`
- scl type I: `A6EXA3`, `BAA32541.1`, `BAF35850.1`, `O82950`, `P12625`, `ZP_01169502.1`
- scl type II: `AAB40611.1`, `BAA35137.1`, `BAA92354.1`, `O05527`, `Q5YEW3`

All 17 records retain their existing P04/P05 accession, database, retrieval,
literature, and SHA256 provenance. No control is promoted to a seed. The
three archived subtype HMMs remain retained, blocked reference artifacts for
downstream sequence architecture and phylogenetic review only.

## Close Non-Target Hydrolase Challenges

The core calibration adds five experimentally characterized bacterial
hydrolases that are not PHA-depolymerase references. Their exact FASTAs,
accession-level provenance, and residue SHA256 values are in
[`p05_extracellular_core_close_controls.tsv`](../04_family_profiles/manifests/p05_extracellular_core_close_controls.tsv).

| Accession | Function used as challenge | Key experimental support |
| --- | --- | --- |
| `A6WFI5` | *Kineococcus* cutinase KrCUT | Hydrolyzes cutin/esters but binds rather than hydrolyzes PHB (PMID `34705546`). |
| `P0DX29` | *Amycolatopsis* cutinase/lipase | Ester, triglyceride and synthetic-polyester activity (PMIDs `19806375`, `21145735`, `33598102`). |
| `Q47RJ6` | *Thermobifida* cutinase/PET hydrolase | Cutin, ester and PET-hydrolase assays (PMIDs `18658138`, `20729325`, `21594592`, `21751386`, `23603671`, `25545638`). |
| `O53581` | Mycobacterial Culp6 esterase/phospholipase | Direct ester and phospholipid hydrolysis assays (PMIDs `19169353`, `19225166`, `20656688`, `29247008`). |
| `P9WP43` | Mycobacterial Culp1 carboxylesterase | Direct short-/medium-chain ester and weak lipase assays (PMIDs `16716602`, `19225166`, `20103719`, `23843969`). |

These are hard specificity challenges, not phenotype-negative training labels.
An accepted core rule must reject every one, as well as every seed from the
three non-extracellular active P05 branches.

## Required Acceptance Evidence

The reproducible helper is
[`p05_extracellular_core.py`](../scripts/p05_extracellular_core.py). It uses
MAFFT L-INS-i and HMMER 3.4 command forms already reviewed in
[`P05_HMM_CALIBRATION_PROTOCOL.md`](P05_HMM_CALIBRATION_PROTOCOL.md).

Before the core can change the P06 registry, the isolated T141 worktree must
record all of the following:

1. HMM, pooled bundle, and alignment SHA256 values, with HMM header `NSEQ=17`.
2. `17/17` independent leave-one-out positive recoveries.
3. A full-model control smoke against 20 hard challenges: 15 seeds from the
   three non-extracellular P05 branches plus the five close non-target
   hydrolases above.
4. No hard challenge may pass the retained score-and-HMM-coverage rule.
5. A human-reviewed registry update that approves only
   `extracellular_pha_depolymerase_core`; the three subtype models remain
   blocked from direct P06 scanning.

P06 candidates from the core model remain sequence candidates. Localization,
signal peptide, architecture, and phylogenetic placement are subsequent
evidence layers; none establishes a PHA-degradation phenotype or a substrate
subtype by itself.

## T141 Calibration Result

The isolated T141 r7 worktree built the core with MAFFT `v7.525` L-INS-i and
HMMER `3.4`. The HMM header reports `LENG=582`, `NSEQ=17`, and
`EFFN=2.618896`. Artifact SHA256 values are bundle
`208ba02978158968563ca4de0a41dabebb9ae0e8bde5cbc16cacc1cbeb77f9e7`,
alignment `3f5a475dfd4de86a03c2f4f009009302765410adae7f043299cdedb80bd036da`,
and HMM `74c4b69a2d845f0725d0bc348402e6a51ba3c17a9f67f8cabed3b63df6a6e2f4`.

All `17/17` held-out positives were recovered. The final rule is full score
`>=159.9` and HMM coverage `>=0.405498`. Of 20 hard controls, only `A6WFI5`
reported any HMMER hit (score `1.3`, coverage `0.292096`), so `0/20` controls
pass the rule. The compact accession-level results are
`p05_extracellular_core_leave_one_out_positive_results.tsv` and
`p05_extracellular_core_control_smoke_results.tsv`; the approved decision is
in `p05_extracellular_core_calibration_decision_summary.tsv`.
