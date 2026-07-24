# P02 Benchmark Decision

**Run date:** 2026-07-24 Asia/Shanghai

P02 confirms the single production predictor route:

```text
Pyrodigal GeneFinder(meta=True)
```

## Inputs

- Raw genomes: `/home/data/haoyu/PHB-GTDB-GPT/00_raw_gtdb_r232/genomes`
- Bacterial taxonomy: `/home/data/haoyu/PHB-GTDB-GPT/00_raw_gtdb_r232/bac120_taxonomy_r232.tsv`
- Archaeal taxonomy: `/home/data/haoyu/PHB-GTDB-GPT/00_raw_gtdb_r232/ar53_taxonomy_r232.tsv`
- Selection output: `/home/data/haoyu/PHB-GTDB-GPT/02_prediction_benchmark/qc/benchmark_genomes.tsv`
- Metrics output: `/home/data/haoyu/PHB-GTDB-GPT/02_prediction_benchmark/qc/pyrodigal_meta_metrics.tsv`

Generated P02 TSV outputs remain ignored by Git.

## Benchmark Set

- Selected genomes: 240
- Selection method: deterministic GTDB-taxonomy-stratified sampling
- Random seed: 20260724
- Bacterial strata: phylum-level
- Archaeal strata: class-level

## Pyrodigal Meta-Mode Metrics

- Genomes processed: 240
- Successful genomes: 240
- Failed genomes: 0
- Coding density minimum: 0.730487
- Coding density mean: 0.899989
- Coding density maximum: 0.976991
- Predicted gene count minimum: 481
- Predicted gene count mean: 2454.30
- Predicted gene count maximum: 9286
- Short ORFs total: 2257
- Internal stops reported by the P02 metric script: 0
- Illegal amino acids reported by the P02 metric script: 0

## Decision

The P02 benchmark accepts `Pyrodigal GeneFinder(meta=True)` as the single production route for P03 complete GTDB proteome prediction.

P03 may proceed using this policy:

```text
config/prediction_policy.yaml
```

Every P03 protein identifier must remain reversible:

```text
{gtdb_genome_accession}|{contig_id}|{orf_index}
```
