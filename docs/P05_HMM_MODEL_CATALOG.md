# P05 HMM Model Catalog

**Catalog date:** 2026-07-27
**Purpose:** publish checksum-locked provenance for the machine-local P05 HMMs without committing the raw HMMs, seed bundles, or alignments.

## Scope

The tracked registries record the exact local profile, seed bundle, and alignment used for each archived or rebuilt model. HMMER scores, domains, and tree placement remain sequence evidence only; they do not prove a PHB/PHA degradation phenotype for a GTDB hit.

The governing biological decisions are in [P05 HMM seed selection decision](P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md). The raw artifacts remain ignored by Git, while their provenance is published in:

- [`p05_hmm_model_registry.tsv`](../04_family_profiles/manifests/p05_hmm_model_registry.tsv): HMMER headers, HMM/bundle/alignment SHA256 values, recorded `hmmbuild` command, calibration state, and P06 permission.
- [`p05_hmm_seed_registry.tsv`](../04_family_profiles/manifests/p05_hmm_seed_registry.tsv): all 32 current bundle members with accession, organism, evidence level, source/version/date, DOI/PMID/PMCID, sequence FASTA SHA256, and rebuild role.
- [`p05_hmm_proposed_seed_updates.tsv`](../04_family_profiles/manifests/p05_hmm_proposed_seed_updates.tsv): 26 implemented profile-seed or retained boundary-control decisions, each with a source FASTA SHA256.

`scripts/p05_catalog_hmm_models.py` recreates these files and verifies that each HMM `NSEQ` equals its seed-bundle header count.

## Current Profiles

All six profiles are HMMER `HMMER3/f` version `3.4` artifacts. Three affected
families were rebuilt with MAFFT `v7.525` L-INS-i followed by `hmmbuild
--amino` on T141 in the isolated
`/home/data/haoyu/PHB-GTDB-GPT-p05-rebuild-20260727` worktree. The three
extracellular models are retained archived artifacts. No profile contains
model-specific `GA`, `TC`, or `NC` thresholds, so calibration against close
non-target hydrolases remains incomplete. All rows are `approved_for_p06=no`
and `scan_permission=blocked`.

| Family | Status | NSEQ | HMM length | EFFN | Model SHA256 |
|---|---|---:|---:|---:|---|
| `archaeal_patatin_like_pha_dep` | `rebuilt_pending_calibration` | 6 | 324 | 0.612305 | `4d0fd5a38e8465834e1e559e99d66162c432c4ccc6eaa36b0953f00137e51582` |
| `extracellular_mcl_pha_dep` | `provisional_archived_pending_calibration` | 6 | 337 | 0.893555 | `75bf6bfabb660b775ba382b60cd487ab80175b88aa9500263f8a086f834218e8` |
| `extracellular_scl_pha_dep_type_I` | `provisional_archived_pending_calibration` | 6 | 569 | 1.133789 | `8d38501c7ecc860bfb416dd46b1a6bca4151e46a133b14b20ec5ca47bf537744` |
| `extracellular_scl_pha_dep_type_II` | `provisional_archived_pending_calibration` | 5 | 502 | 0.913086 | `d90af7a8879dce760b4b80ef18691ce1f6a25f59c852affb333cf0929cc89e0b` |
| `intracellular_mcl_pha_dep` | `rebuilt_pending_calibration` | 4 | 284 | 0.386719 | `d7628183e88204ecac3f2d165eff5ad78fa087aa2ca09eccf1809ff53e72b933` |
| `intracellular_phaZ_no_lipase_box` | `rebuilt_pending_calibration` | 5 | 411 | 0.756836 | `9cb33d27edce3af5266d3d80d9a7b6965eeef20deb9f5028c1e427ae5d033457` |

The TSV registry is authoritative for the matching bundle and alignment hashes; a raw HMM is never accepted merely because its filename matches a family name.

## Calibration Control Panel

[`p05_hmm_calibration_control_panel.tsv`](../04_family_profiles/manifests/p05_hmm_calibration_control_panel.tsv)
is the accession-level input contract for P05 calibration. For every target
model it records two deliberately different categories:

- `cross_family_challenge`: a current seed from another active profile family.
  It is a hard sequence-specificity challenge and must fail the final calibrated
  threshold. This does not describe its source sequence as phenotype-negative.
- `boundary_observation`: a target-family `boundary_candidate` retained for
  architecture/evidence review. It is reported but never used to choose a
  rejection threshold, because its biological family membership remains
  unresolved by design.

The panel carries the exact source accession, organism, database/release,
retrieval date, literature identifiers, FASTA path, and sequence SHA256. It is
generated deterministically from the P04 reference manifest plus the current
P05 model/seed registries:

```powershell
python scripts/p05_hmm_calibration.py
```

For an executable-on-T141 but deliberately non-executed calibration plan, use:

```powershell
python scripts/p05_hmm_calibration.py --build-commands
```

This creates one checksum-verified control FASTA and one deterministic
`hmmsearch --noali --acc --seed 42 --cpu 1 --domtblout` command per model under
the ignored `04_family_profiles/calibration/` directory. The initial commands
are a control-panel smoke check; they do not themselves set a P06 threshold.

Leave-one-out seed models and the full-model control smoke have now been
parsed into compact tracked results:
[`p05_hmm_leave_one_out_positive_results.tsv`](../04_family_profiles/manifests/p05_hmm_leave_one_out_positive_results.tsv),
[`p05_hmm_control_smoke_results.tsv`](../04_family_profiles/manifests/p05_hmm_control_smoke_results.tsv),
and [`p05_hmm_calibration_decision_summary.tsv`](../04_family_profiles/manifests/p05_hmm_calibration_decision_summary.tsv).
Three models are eligible only for human review; the other three are blocked
by positive-recovery failure or cross-family overlap. No result changes
`approved_for_p06=no` or `scan_permission=blocked`.

## Family Decisions

`intracellular_mcl_pha_dep` remains *Pseudomonas*-enriched because this is an intrinsic feature of the known intracellular mcl-PHA mobilization family, not seed bias. Its rebuilt bundle is `Q5Y152`, `B7UCC9`, `Q88D24`, and `Q9R9W3`; `Q5Y152` is the primary experimental anchor. `Q5Q135`, `Q8VV57`, `Q9AGB5`, and `Q9Z3Y0` are retained only as boundary controls. Do not add weak non-*Pseudomonas* sequences merely to force cross-genus balance.

The rebuilt `archaeal_patatin_like_pha_dep` bundle contains `AFK21580.1`, `CCQ36014.1`, `CCQ32286.1`, `AGN01047.1`, `KYH27761.1`, and `ELY43313.1`. The four E3 additions are architecture coverage only, not phenotype evidence. Six longer PHB-synthase-like/HHH/PKD/AxeA-like accessions remain boundary controls outside the alignment.

For `intracellular_phaZ_no_lipase_box`, retain `O87189`, `Q0K7T2`, `Q71KW6`, and `Q92TD3`; add *Paracoccus denitrificans* `Q9WX79` (PMID `11267773`; DOI `10.1111/j.1574-6968.2001.tb10558.x`); and keep `Q0K4D5` only as a labelled boundary/control sequence.

## Scan Enforcement

`scripts/p06_scan_family_profiles.py` requires the model registry. It accepts only rows with `approved_for_p06=yes` and `scan_permission=approved`, recomputes the local HMM SHA256, and writes that value into every P06 scan-manifest row. It rejects missing or empty approval, a model outside `--hmm-dir`, and a hash mismatch.

Before a new scan: calibrate every model against close non-target hydrolases; record thresholds and the calibration decision; then mark only approved rows for P06 and pass single-file plus multi-file streaming smoke tests.

## Reproduction

```powershell
python scripts/p05_catalog_hmm_models.py
python -m unittest tests/test_p05_catalog_hmm_models.py -v
python -m unittest tests/test_p06_scan_family_profiles.py -v
```

Current artifacts used MAFFT L-INS-i followed by HMMER `hmmbuild --amino`. HMMER 3.4 source/user guide: https://github.com/EddyRivasLab/hmmer (BSD-3-Clause).
