# PHB-GTDB-GPT

Reproducible, family-resolved analysis of PHB/PHA depolymerase-related genes in GTDB Release 11 R232.

The project combines experimentally supported reference sequences, DED family definitions, TIGRFAM profiles, custom HMMs, domain-architecture review, and family-specific phylogenetic analysis. A detected homolog is reported as sequence evidence; it is not treated as proof of PHB-degradation phenotype.

## Workflow Stages

| Stage | Purpose |
| --- | --- |
| P01 | Audit and physically copy unfiltered GTDB R232 inputs on T141. |
| P02 | Benchmark Pyrodigal meta mode and Prodigal modes on a deterministic GTDB subset. |
| P03 | Predict and quality-control complete GTDB proteomes. |
| P04 | Build an auditable experimental reference library and family definitions. |
| P05 | Build and calibrate family HMMs. |
| P06 | Scan GTDB proteomes with HMMER and apply family-specific gates. |
| P07 | Review InterPro domain architecture and localization evidence. |
| P08 | Build family trees and join candidates to GTDB taxonomy and species trees. |
| P09 | Produce tables, figures, reports, and provenance manifests. |
| P10 | Optionally validate compatible metabolic models with COBRApy. |

## Execution Environment

The source repository is maintained locally and mirrored to GitHub. Compute-intensive execution occurs on T141 under `/home/data/haoyu/PHB-GTDB-GPT`; see `config/paths.example.yaml` before creating the server-side configuration.

Raw GTDB files, predicted proteomes, HMMER tables, InterPro results, Nextflow work directories, and generated reports are deliberately ignored by Git. They remain auditable through manifests produced by the workflow.

## Initial Validation

```powershell
python -m unittest tests/test_repository_layout.py -v
python scripts/validate_repository.py
```

The Nextflow workflow requires Java 17+ and Nextflow 24.10.0 or newer on its execution host.
