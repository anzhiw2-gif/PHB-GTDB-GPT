# P02 Prediction Benchmark

Deterministic stratified genome selection and Pyrodigal meta-mode quality metrics will be stored here.

P02 is locked to one production route:

```text
Pyrodigal GeneFinder(meta=True)
```

The benchmark scripts confirm that this route can parse the selected GTDB assemblies and produce expected quality metrics. They do not compare against `Prodigal -p meta` or `Prodigal -p single`.

Generated P02 tables and metric reports should be written under ignored subdirectories such as `02_prediction_benchmark/qc/` or `02_prediction_benchmark/predictions/`.
