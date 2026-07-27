# P05 HMM Model Catalog

**Catalog date:** 2026-07-27
**Purpose:** publish checksum-locked provenance for the machine-local P05 HMMs without committing the raw HMMs, seed bundles, or alignments.

## Scope

The tracked registries record the exact local profile, seed bundle, and alignment used for each archived model. HMMER scores, domains, and tree placement remain sequence evidence only; they do not prove a PHB/PHA degradation phenotype for a GTDB hit.

The governing biological decisions are in [P05 HMM seed selection decision](P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md). The raw artifacts remain ignored by Git, while their provenance is published in:

- [`p05_hmm_model_registry.tsv`](../04_family_profiles/manifests/p05_hmm_model_registry.tsv): HMMER headers, HMM/bundle/alignment SHA256 values, recorded `hmmbuild` command, calibration state, and P06 permission.
- [`p05_hmm_seed_registry.tsv`](../04_family_profiles/manifests/p05_hmm_seed_registry.tsv): all 38 current bundle members with accession, organism, evidence level, source/version/date, DOI/PMID/PMCID, sequence FASTA SHA256, and rebuild role.
- [`p05_hmm_proposed_seed_updates.tsv`](../04_family_profiles/manifests/p05_hmm_proposed_seed_updates.tsv): explicit candidates that are not yet retrieved or added to a profile.

`scripts/p05_catalog_hmm_models.py` recreates these files and verifies that each HMM `NSEQ` equals its seed-bundle header count.

## Current Archived Profiles

All six profiles are HMMER `HMMER3/f` version `3.4` artifacts. No profile contains model-specific `GA`, `TC`, or `NC` thresholds, so calibration against close non-target hydrolases remains incomplete. All rows are `approved_for_p06=no` and `scan_permission=blocked`.

| Family | Status | NSEQ | HMM length | Model SHA256 |
|---|---|---:|---:|---|
| `archaeal_patatin_like_pha_dep` | `blocked_for_rebuild` | 8 | 454 | `3fc97d9feb4de35906e921aa53607bc0691d8d4ef9416f0a1133e05a23842a21` |
| `extracellular_mcl_pha_dep` | `provisional_archived_pending_calibration` | 6 | 337 | `75bf6bfabb660b775ba382b60cd487ab80175b88aa9500263f8a086f834218e8` |
| `extracellular_scl_pha_dep_type_I` | `provisional_archived_pending_calibration` | 6 | 569 | `8d38501c7ecc860bfb416dd46b1a6bca4151e46a133b14b20ec5ca47bf537744` |
| `extracellular_scl_pha_dep_type_II` | `provisional_archived_pending_calibration` | 5 | 502 | `d90af7a8879dce760b4b80ef18691ce1f6a25f59c852affb333cf0929cc89e0b` |
| `intracellular_mcl_pha_dep` | `provisional_archived_needs_row_audit` | 8 | 285 | `14ac60e28f78cdf54d6276b5368038702b1c35f2c4f09cb429e4e5634be5bd6f` |
| `intracellular_phaZ_no_lipase_box` | `proposed_seed_update_pending_user_confirmation` | 5 | 408 | `fa936a7b38dc30b5aa50a1dcd715074d8b4909501d5963e2c7ffbe662076abc7` |

The TSV registry is authoritative for the matching bundle and alignment hashes; a raw HMM is never accepted merely because its filename matches a family name.

## Family Decisions

`intracellular_mcl_pha_dep` remains *Pseudomonas*-enriched because this is an intrinsic feature of the known intracellular mcl-PHA mobilization family, not seed bias. `Q5Y152` (*Pseudomonas putida* KT2442; PMID `17170116`; DOI `10.1074/jbc.M608119200`) is the primary experimental anchor. The other current rows require accession-level experimental-evidence audit; do not add weak non-*Pseudomonas* sequences merely to force cross-genus balance.

The archived `archaeal_patatin_like_pha_dep` model is architecture-mixed and blocked for rebuild. The coherent proposal retains `AFK21580.1` (PhaZh1; PMID `25710370`; DOI `10.1128/AEM.04269-14`) and `CCQ36014.1`, adds E3 patatin-like coverage candidates `CCQ32286.1`, `AGN01047.1`, `KYH27761.1`, and `ELY43313.1`, and demotes six longer PHB-synthase-like/HHH/PKD/AxeA-like rows. Those archaeal E3 rows are coverage evidence, not phenotype evidence.

For `intracellular_phaZ_no_lipase_box`, retain `O87189`, `Q0K7T2`, `Q71KW6`, and `Q92TD3`; add *Paracoccus denitrificans* `Q9WX79` (PMID `11267773`; DOI `10.1111/j.1574-6968.2001.tb10558.x`); and keep `Q0K4D5` only as a labelled boundary/control sequence.

## Scan Enforcement

`scripts/p06_scan_family_profiles.py` requires the model registry. It accepts only rows with `approved_for_p06=yes` and `scan_permission=approved`, recomputes the local HMM SHA256, and writes that value into every P06 scan-manifest row. It rejects missing or empty approval, a model outside `--hmm-dir`, and a hash mismatch.

Before a new scan: confirm the revised seed set; retrieve and log all approved new accessions; rebuild or explicitly retain profiles; record HMM/bundle/alignment hashes; make the calibration decision; then mark only approved rows for P06 and pass single-file plus multi-file streaming smoke tests.

## Reproduction

```powershell
python scripts/p05_catalog_hmm_models.py
python -m unittest tests/test_p05_catalog_hmm_models.py -v
python -m unittest tests/test_p06_scan_family_profiles.py -v
```

Current artifacts used MAFFT L-INS-i followed by HMMER `hmmbuild --amino`. HMMER 3.4 source/user guide: https://github.com/EddyRivasLab/hmmer (BSD-3-Clause).
