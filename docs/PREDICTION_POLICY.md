# Prediction Policy

P02 uses a single predictor route:

- `Pyrodigal GeneFinder(meta=True)`

This is the only predictor route to be used for the GTDB-wide proteome prediction stage. The project does not benchmark or switch between `Prodigal -p meta` and `Prodigal -p single` for the production pipeline.

Why this choice:

- It is the metagenomic mode exposed by the Python-native interface.
- It is the most consistent choice for a large GTDB collection that contains assemblies of mixed completeness and fragmentation.
- It keeps the implementation path uniform and easy to wire into the later Nextflow stages.

P02 still records benchmark evidence and quality metrics, but only to confirm that the chosen route meets the project's acceptance criteria.
