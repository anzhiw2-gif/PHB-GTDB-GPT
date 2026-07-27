# P05 HMM Calibration Protocol

**Protocol date:** 2026-07-27

## Purpose And Claim Boundary

P05 calibration asks whether a checksum-locked family HMM separates its
experimentally anchored/architecture-coherent seed set from close sequence
challenges. The resulting score and coverage limits are sequence-classification
criteria only. They do not demonstrate PHB/PHA depolymerase activity for a
GTDB hit.

P06 remains blocked until each of the six profiles has a documented calibration
decision and the model registry is explicitly updated after review.

## Input Contract

The tracked
[`p05_hmm_calibration_control_panel.tsv`](../04_family_profiles/manifests/p05_hmm_calibration_control_panel.tsv)
is generated from the P04 reference manifest, current P05 seed registry, and
checksum-locked HMM registry by:

```powershell
python scripts/p05_hmm_calibration.py --build-commands
```

It creates two categories for every target family:

| Category | Threshold role | Biological interpretation |
|---|---|---|
| `cross_family_challenge` | hard challenge; must fail the final rule | A seed from another active P05 family, used to test sequence-family separation. It is not called phenotype-negative. |
| `boundary_observation` | report only; excluded from threshold selection | A target-family record intentionally excluded from the profile for weak evidence or incompatible architecture. Its true membership is unresolved. |

Thus the six non-patatin archaeal records, four evidence-limited
*Pseudomonas* mcl-PHA records, and `Q0K4D5` are not silently treated as
negative training labels.

## Current Smoke Panel

The initial generated panel has 171 target-model/control records:

- 160 `cross_family_challenge` records
- 11 `boundary_observation` records
- six planned `hmmsearch` jobs, one per checksum-locked full HMM

The generated FASTAs, command manifest, raw `domtblout`, and logs remain under
ignored `04_family_profiles/calibration/`. The compact, accessioned input panel
is tracked in Git.

## Planned Acceptance Decision

The full calibration will use leave-one-out seed models so that a retained seed
is scored by a profile that did not include it. A family may be considered for
P06 only when all conditions below are documented:

1. Every leave-one-out positive has a retained, full-length-enough HMM match.
2. A score-and-HMM-coverage rule separates every hard cross-family challenge.
3. Boundary observations are reported separately and never used to manufacture
   a more stringent cutoff.
4. The final full model also passes its control-panel smoke check using the
   derived rule.
5. The result, threshold provenance, tool version, model SHA256, and failed or
   retained boundary observations are recorded before a human review changes
   `approved_for_p06` or `scan_permission`.

If a family fails separation, it remains blocked; a weaker threshold is not
used to force P06 approval.

## Software And Implementation Review

The command form follows HMMER `hmmsearch` HMMER 3.4, using `--domtblout` for
machine-readable domain coordinates plus `--noali --acc --seed 42 --cpu 1` for
deterministic, compact output. The source and documentation review was
performed on 2026-07-27:

- HMMER source/User Guide: https://github.com/EddyRivasLab/hmmer
  - reviewed release commit: `9acd8b6758a0ca5d21db6d167e0277484341929b`
    (`release-3.4`, 2023-08-15)
  - HMMER version: `3.4` (the repository `configure.ac` identifies August 2023)
  - license: BSD open-source license as stated by HMMER source `configure.ac`
  - rationale: authoritative implementation for the profile-HMM build/search
    format already used by the checksum-locked P05 models.
- HMMER documentation: https://hmmer.org/documentation.html
  - rationale: official command and User Guide landing page; used for the
    `hmmsearch`/`hmmbuild` argument order and output convention.

GitHub repository discovery for `hmmsearch calibration negative controls` did
not identify a maintained, directly reusable calibration workflow on
2026-07-27. The small project-specific planner is therefore deliberately
limited to provenance validation, deterministic FASTA materialization, and
HMMER command generation; it does not copy unverifiable third-party code.

## Reproduction Checks

```powershell
python -m unittest tests/test_p05_hmm_calibration.py -v
python scripts/p05_hmm_calibration.py --build-commands
```
