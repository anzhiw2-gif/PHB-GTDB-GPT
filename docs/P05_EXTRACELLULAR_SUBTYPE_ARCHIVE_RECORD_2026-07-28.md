# P05 Extracellular Subtype Archive And Core HMM Record

**Record date:** 2026-07-28
**Stage:** P05 family-profile construction, calibration, and P06 handoff
**Scope:** three archived extracellular subtype HMMs and their validated pooled replacement

## Purpose And Claim Boundary

This record preserves the full construction and validation trail for the three
initial extracellular PHA-depolymerase branches:
`extracellular_mcl_pha_dep`, `extracellular_scl_pha_dep_type_I`, and
`extracellular_scl_pha_dep_type_II`. It also records their accession-preserving
handover to the approved `extracellular_pha_depolymerase_core` P06 model.

An HMM hit supports sequence-level relatedness only. It does not establish
PHA-depolymerase activity, extracellular localization, mcl versus scl substrate
preference, or type-I/type-II subtype for a GTDB protein.

## Seed Admission And Archived Model Construction

All 17 extracellular seeds are bacterial and were admitted only after
accession-level experimental support. The source database/release, organism,
retrieval date, literature identifiers, FASTA path, and residue SHA256 are
retained in `p05_hmm_seed_registry.tsv` and the P04 reference manifest.

| Archived subtype | Seed accessions | Evidence | HMM header | HMM SHA256 |
|---|---|---|---|---|
| `extracellular_mcl_pha_dep` | `BBB62377.1`, `BBB68334.1`, `CAQ46066.1`, `Q51718`, `Q6UFW4`, `Q84C08` | E1=6 | LENG=337; NSEQ=6; EFFN=0.893555 | `75bf6bfabb660b775ba382b60cd487ab80175b88aa9500263f8a086f834218e8` |
| `extracellular_scl_pha_dep_type_I` | `A6EXA3`, `BAA32541.1`, `BAF35850.1`, `O82950`, `P12625`, `ZP_01169502.1` | E1=5; E2=1 | LENG=569; NSEQ=6; EFFN=1.133789 | `8d38501c7ecc860bfb416dd46b1a6bca4151e46a133b14b20ec5ca47bf537744` |
| `extracellular_scl_pha_dep_type_II` | `AAB40611.1`, `BAA35137.1`, `BAA92354.1`, `O05527`, `Q5YEW3` | E1=5 | LENG=502; NSEQ=5; EFFN=0.913086 | `d90af7a8879dce760b4b80ef18691ce1f6a25f59c852affb333cf0929cc89e0b` |

The archived profiles were generated as machine-local P05 artifacts with MAFFT
L-INS-i (`mafft --localpair --maxiterate 1000 --inputorder`) followed by HMMER
3.4 `hmmbuild --amino`. Their recorded build dates are 2026-07-26. The exact
seed-bundle SHA256, alignment SHA256, model path, and recorded `hmmbuild`
command for each branch are authoritative in `p05_hmm_model_registry.tsv`.
Raw HMMs, bundles, and alignments remain ignored because they are generated
artifacts; the hashes make their identity independently checkable.

## Artifact Verification And Subtype Calibration

On 2026-07-27, each archived HMM was verified in isolated T141 calibration
worktree `PHB-GTDB-GPT-p05-calibration-r3-20260727`: `sha256sum` and HMM header
fields (`NAME`, `LENG`, `NSEQ`, and `EFFN`) had to match the tracked registry.
Each profile was evaluated by leave-one-out recovery and by hard cross-family
seed challenges. A retained subtype rule must recover every held-out seed and
reject every other active-family seed.

| Archived subtype | Leave-one-out recovery | Strictest retained rule | Hard cross-family seeds passing | P06 outcome |
|---|---:|---|---:|---|
| `extracellular_mcl_pha_dep` | 6/6 | score >=35.5; HMM coverage >=0.222222 | 7 | blocked |
| `extracellular_scl_pha_dep_type_I` | 6/6 | score >=105.7; HMM coverage >=0.716495 | 2 | blocked |
| `extracellular_scl_pha_dep_type_II` | 4/5 | none | not evaluated | blocked |

The experimentally supported `CAQ46066.1` is a short, divergent mcl seed, so
its recovery fixes the low mcl retained threshold. Removing it merely to force
subtype separation would violate the bacterial accession-level experimental
seed rule. Type II cannot recover the experimentally supported 635-aa
`AAB40611.1` when held out; its remaining references are 491-495 aa. These
are observed architecture and family-overlap limits, not a threshold-tuning
problem.

The accession-level leave-one-out and hard-control results are retained in
`p05_hmm_leave_one_out_positive_results.tsv`,
`p05_hmm_control_smoke_results.tsv`, and
`p05_hmm_calibration_decision_summary.tsv`.

## Pooled Core Construction And Approval

The 17 unchanged, experimentally supported seed accessions above were pooled
in isolated T141 r7 worktree
`/home/data/haoyu/PHB-GTDB-GPT-p05-extracellular-core-r7-20260727` to build
`extracellular_pha_depolymerase_core` with the same MAFFT and HMMER versions.
Its final HMM header is `LENG=582`, `NSEQ=17`, and `EFFN=2.618896`.

| Artifact | SHA256 |
|---|---|
| Seed bundle | `208ba02978158968563ca4de0a41dabebb9ae0e8bde5cbc16cacc1cbeb77f9e7` |
| Alignment | `3f5a475dfd4de86a03c2f4f009009302765410adae7f043299cdedb80bd036da` |
| HMM | `74c4b69a2d845f0725d0bc348402e6a51ba3c17a9f67f8cabed3b63df6a6e2f4` |

All 17 leave-one-out positives were recovered. The final parser rule is full
score `>=159.9` and HMM coverage `>=0.405498`. The 20 hard challenges include
the 15 non-extracellular P05 seeds and five accessioned close non-target
hydrolases. No challenge passed the final rule; `A6WFI5` had only a low raw hit
(score 1.3; coverage 0.292096). The final core seed registry now binds every
one of the 17 accessions to the final HMM SHA256, and
`tests/test_p05_extracellular_provenance.py` checks this relationship.

## P06 Permission And Downstream Use

`p05_hmm_model_registry.tsv` is the P06 authority. The three archived subtype
rows have `approved_for_p06=no` and `scan_permission=blocked`; they remain
subtype-reference artifacts for P07 architecture review and P08 phylogenetic
interpretation. The core row alone has `approved_for_p06=yes` and
`scan_permission=approved` for extracellular discovery.

P06 may identify extracellular-core sequence candidates only. Candidate
subtype, localization, and phenotype interpretation require the downstream
evidence layers and must not be inferred from the core HMM hit alone.

## Related Tracked Records

- `docs/P05_EXTRACELLULAR_CORE_DECISION_2026-07-27.md`: approval rationale and close-hydrolase controls.
- `docs/P05_HMM_CONSTRUCTION_VALIDATION_RECORD_2026-07-27.md`: all six original branches, validation protocol, and four-model P06 gate.
- `docs/P05_HMM_MODEL_CATALOG.md`: complete model, alignment, bundle, and seed-registry index.
- `04_family_profiles/manifests/p05_extracellular_core_*.tsv`: core seeds, controls, leave-one-out results, and calibration summary.
