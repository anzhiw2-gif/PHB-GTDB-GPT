# P05 HMM Construction And Validation Record

**Record date:** 2026-07-27
**Stage:** P05 family-profile construction and calibration
**GTDB analysis target:** Release 11 R232
**Execution host:** T141 (`/home/data/haoyu`)
**P06 state:** four calibration-approved models are ready for a new scan; no
GTDB-wide scan had been launched at the time this P05 validation record was updated.

## Purpose And Claim Boundary

This record is the single audit entry point for the six custom P05 HMMs. It
connects accession-level seed provenance, biological inclusion rules, generated
artifact identities, calibration inputs, execution environments, and the
resulting scan gate.

An HMM match, a motif, a domain architecture, or tree placement is sequence
evidence only. None of these records establish PHB/PHA depolymerase activity
for a GTDB protein. In particular, archaeal E3 architecture-coverage seeds are
not phenotype evidence, and boundary/control sequences are not biological
negatives.

## Evidence And Seed Admission Rules

The user-approved decision document is
[P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md](P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md).
It locks the following rules before a sequence may become a profile seed:

1. Bacterial seeds require accession-level experimental support. Annotation-only
   bacterial E3 records are excluded from profile alignment.
2. Archaeal E3 records are admissible only as explicitly labelled
   architecture-coverage seeds when direct experimental coverage is sparse.
3. Each profile must retain a coherent architecture; an automatic product name
   cannot override conflicting motifs, domains, or alignment placement.
4. Taxonomic balance must not be manufactured with weak seeds. The
   *Pseudomonas*-enriched intracellular mcl-PHA branch is retained because the
   observed biology supports this enrichment.
5. Boundary candidates are tracked as report-only observations and never used
   as negative labels or as a way to choose a more stringent rejection cutoff.

The exact source database/release, retrieval date, organism, evidence level,
DOI/PMID/PMCID, local seed FASTA location, and sequence SHA256 for all current
seeds are in
[`p05_hmm_seed_registry.tsv`](../04_family_profiles/manifests/p05_hmm_seed_registry.tsv).
The 26 seed-admission, demotion, and boundary decisions are retained in
[`p05_hmm_proposed_seed_updates.tsv`](../04_family_profiles/manifests/p05_hmm_proposed_seed_updates.tsv).

## Current Seed Bundles

| Family | Current profile-seed accessions | Evidence composition | Biological scope |
|---|---|---|---|
| `archaeal_patatin_like_pha_dep` | `AFK21580.1`, `AGN01047.1`, `CCQ32286.1`, `CCQ36014.1`, `ELY43313.1`, `KYH27761.1` | E1=1; E3=5 | `AFK21580.1` is the experimental PhaZh1 anchor. The five E3 records provide coherent patatin-like `GxSxG` architecture coverage only. |
| `extracellular_mcl_pha_dep` | `BBB62377.1`, `BBB68334.1`, `CAQ46066.1`, `Q51718`, `Q6UFW4`, `Q84C08` | E1=6 | Experimentally supported extracellular mcl-PHA branch. |
| `extracellular_scl_pha_dep_type_I` | `A6EXA3`, `BAA32541.1`, `BAF35850.1`, `O82950`, `P12625`, `ZP_01169502.1` | E1=5; E2=1 | Experimentally supported extracellular scl-PHA type-I branch. |
| `extracellular_scl_pha_dep_type_II` | `AAB40611.1`, `BAA35137.1`, `BAA92354.1`, `O05527`, `Q5YEW3` | E1=5 | Experimentally supported extracellular scl-PHA type-II branch. |
| `intracellular_mcl_pha_dep` | `B7UCC9`, `Q5Y152`, `Q88D24`, `Q9R9W3` | E1=2; E2=2 | *Pseudomonas*-enriched intracellular mcl-PHA family; `Q5Y152` is the primary direct experimental anchor. |
| `intracellular_phaZ_no_lipase_box` | `O87189`, `Q0K7T2`, `Q71KW6`, `Q92TD3`, `Q9WX79` | E1=1; E2=4 | Intracellular no-lipase-box branch; `Q9WX79` provides an experimentally supported cross-genus seed. |

The following records remain report-only boundaries rather than alignment
members: six non-patatin archaeal records, `Q5Q135`, `Q8VV57`, `Q9AGB5`,
`Q9Z3Y0` for intracellular mcl-PHA, and `Q0K4D5` for the no-lipase-box model.
Their accession-level rationale is preserved in the seed decision and proposed
update manifest.

## Artifact Construction And Identity

The three revised models were rebuilt in the isolated
`/home/data/haoyu/PHB-GTDB-GPT-p05-rebuild-20260727` worktree with MAFFT
`v7.525` L-INS-i (`--localpair --maxiterate 1000 --inputorder`) followed by
HMMER `3.4` `hmmbuild --amino`. The three extracellular models are retained
archived artifacts and were not represented as newly rebuilt models.

All models are `HMMER3/f` HMMER `3.4` artifacts. The model registry records
the exact bundle, alignment, and model SHA256 values and the exact recorded
`hmmbuild` command. Raw HMMs, seed bundles, and alignments remain machine-local
because they are generated artifacts; their content identity is preserved by
these hashes.

| Family | Build status | HMM length | NSEQ | EFFN | Model SHA256 |
|---|---|---:|---:|---:|---|
| `archaeal_patatin_like_pha_dep` | rebuilt | 324 | 6 | 0.612305 | `4d0fd5a38e8465834e1e559e99d66162c432c4ccc6eaa36b0953f00137e51582` |
| `extracellular_mcl_pha_dep` | archived | 337 | 6 | 0.893555 | `75bf6bfabb660b775ba382b60cd487ab80175b88aa9500263f8a086f834218e8` |
| `extracellular_scl_pha_dep_type_I` | archived | 569 | 6 | 1.133789 | `8d38501c7ecc860bfb416dd46b1a6bca4151e46a133b14b20ec5ca47bf537744` |
| `extracellular_scl_pha_dep_type_II` | archived | 502 | 5 | 0.913086 | `d90af7a8879dce760b4b80ef18691ce1f6a25f59c852affb333cf0929cc89e0b` |
| `intracellular_mcl_pha_dep` | rebuilt | 284 | 4 | 0.386719 | `d7628183e88204ecac3f2d165eff5ad78fa087aa2ca09eccf1809ff53e72b933` |
| `intracellular_phaZ_no_lipase_box` | rebuilt | 411 | 5 | 0.756836 | `9cb33d27edce3af5266d3d80d9a7b6965eeef20deb9f5028c1e427ae5d033457` |

### Artifact Verification

On 2026-07-27, the six r3 HMM files in
`/home/data/haoyu/PHB-GTDB-GPT-p05-calibration-r3-20260727` were checked with
`sha256sum` and their HMM headers were checked for `NAME`, `LENG`, `NSEQ`, and
`EFFN`. Every value matched the tracked
[`p05_hmm_model_registry.tsv`](../04_family_profiles/manifests/p05_hmm_model_registry.tsv),
including `NSEQ` equal to the profile-seed count above. A filename alone is
never accepted as model identity.

## Calibration Design

The calibration input panel is
[`p05_hmm_calibration_control_panel.tsv`](../04_family_profiles/manifests/p05_hmm_calibration_control_panel.tsv).
It is generated deterministically from the P04 reference manifest, current
seed registry, and model registry. It contains 171 target-model/control rows:

- 160 `cross_family_challenge` records: seeds from another active family that
  must not pass a final family rule. These are specificity challenges, not
  phenotype-negative labels.
- 11 `boundary_observation` records: accessioned records deliberately excluded
  from a target HMM. They are report-only and cannot define a rejection
  threshold.

Full-model smoke searches were run serially in isolated r2 with:

```bash
hmmsearch --noali --acc --seed 42 --cpu 1 --domtblout <output.domtblout> <model.hmm> <controls.faa>
```

The six full-model control `domtblout` files were copied read-only to isolated
r3 for joint parsing; each source/target file pair had matching SHA256. The
control panel was regenerated in r3 and compared byte-for-byte with r2 before
this parse.

Leave-one-out evaluation then produced 32 variants. Each variant removed
exactly one current seed, retained at least three training sequences, built a
temporary MAFFT/HMMER model, and searched only the held-out positive:

```bash
mafft --localpair --maxiterate 1000 --inputorder <training.faa> > <training.aligned.faa>
hmmbuild --amino <leave_one_out.hmm> <training.aligned.faa>
hmmsearch --noali --acc --seed 42 --cpu 1 --domtblout <holdout.domtblout> <leave_one_out.hmm> <held_out_positive.faa>
```

The parser verifies the expected held-out FASTA target name, parses all HMMER
domain rows, and computes HMM coverage as the union of domain HMM-coordinate
intervals divided by HMM length. A file with no data row is an explicit missing
positive recovery, not an omitted result. It also refuses to recommend a model
when its family lacks a hard cross-family challenge.

## Calibration Results And Scan Gate

The three tracked calculation outputs are:

- [`p05_hmm_leave_one_out_positive_results.tsv`](../04_family_profiles/manifests/p05_hmm_leave_one_out_positive_results.tsv): 32 accession-level held-out positive results.
- [`p05_hmm_control_smoke_results.tsv`](../04_family_profiles/manifests/p05_hmm_control_smoke_results.tsv): 171 control observations with score, HMM coverage, and hit state.
- [`p05_hmm_calibration_decision_summary.tsv`](../04_family_profiles/manifests/p05_hmm_calibration_decision_summary.tsv): family-level proposed score/coverage rule and status.

The proposed rule is the strictest conjunction that still retains every
leave-one-out positive: `full score >= minimum positive score AND HMM coverage
>= minimum positive coverage`. The rule is evaluated against full-model hard
controls; boundary observations do not affect it.

| Model | Leave-one-out recovery | Proposed score / coverage | Hard controls passing | Calibration result |
|---|---:|---:|---:|---|
| `archaeal_patatin_like_pha_dep` | 6/6 | 322.9 / 0.807453 | 0 | Eligible for human review only |
| `extracellular_mcl_pha_dep` | 6/6 | 35.5 / 0.222222 | 7 | Blocked: cross-family overlap |
| `extracellular_scl_pha_dep_type_I` | 6/6 | 105.7 / 0.716495 | 2 | Blocked: cross-family overlap |
| `extracellular_scl_pha_dep_type_II` | 4/5 | not proposed | not evaluated | Blocked: `AAB40611.1` was not recovered |
| `intracellular_mcl_pha_dep` | 4/4 | 617.6 / 0.992933 | 0 | Eligible for human review only |
| `intracellular_phaZ_no_lipase_box` | 5/5 | 415.0 / 0.829384 | 0 | Eligible for human review only |

Boundary observations behaved as expected for unresolved near-family records:
the four intracellular mcl-PHA boundaries and `Q0K4D5` had strong full-model
matches but were not used to set a threshold. The six non-patatin archaeal
boundary records had no smoke hit.

## Current P06 Decision

Human review approved the three individually calibrated models with no hard
challenge passing their retained rule: `archaeal_patatin_like_pha_dep`,
`intracellular_mcl_pha_dep`, and `intracellular_phaZ_no_lipase_box`. The
extracellular subtype models remain blocked because their seed classes overlap.
They are replaced for P06 entry by the checksum-locked
`extracellular_pha_depolymerase_core`; its 17/17 leave-one-out recovery and
0/20 hard-control passing result are recorded in
[P05_EXTRACELLULAR_CORE_DECISION_2026-07-27.md](P05_EXTRACELLULAR_CORE_DECISION_2026-07-27.md).

The model registry now permits exactly four scan models. P06 must carry each
model's calibrated score and HMM-coverage thresholds into its manifest and
High-confidence parsing tier; an HMM hit remains sequence evidence only.

The archived extracellular subtype construction and the accession-preserving
handover to the pooled core are expanded in
[P05_EXTRACELLULAR_SUBTYPE_ARCHIVE_RECORD_2026-07-28.md](P05_EXTRACELLULAR_SUBTYPE_ARCHIVE_RECORD_2026-07-28.md).

## Reproduction And Verification Commands

Run lightweight repository checks locally:

```powershell
python -m unittest tests/test_p05_hmm_calibration.py -v
python -m unittest tests/test_p05_catalog_hmm_models.py -v
python -m unittest tests/test_p06_scan_family_profiles.py -v
python scripts/validate_repository.py
git diff --check
```

On T141, use the project Conda environment and an isolated worktree. Do not
reset, pull, stash, or otherwise modify the live primary worktree:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT-p05-calibration-r3-20260727
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p05_hmm_calibration.py --build-commands
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p05_hmm_calibration.py --build-leave-one-out
```

The tracked control panel and model registry are the reproducible inputs; raw
FASTA targets, alignments, HMMs, HMMER logs, and `domtblout` files remain
ignored machine-local outputs. HMMER source and documentation provenance is
recorded in [P05_HMM_CALIBRATION_PROTOCOL.md](P05_HMM_CALIBRATION_PROTOCOL.md).
