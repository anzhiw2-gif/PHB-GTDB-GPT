# P07 Completion Record

**Date:** 2026-07-29, Asia/Shanghai  
**Stage:** P07, independent domain-architecture and localization review for
P06 GTDB R232 sequence candidates  
**Status:** completed on T141

## What was completed

P07 used the machine-local P06 candidate table from isolated T141 worktree
`/home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728` and completed the full
High-confidence review pass.

Scope of the completed run:

- selected unique candidate sequences: `37,912`
- FASTA shards: `76`
- planned InterProScan command rows: `76`
- planned SignalP6 command rows: `76`
- runner workers: `6`
- InterProScan CPU per job: `10`

The final runner status table reported:

- total rows: `152`
- `completed: 103`
- `skipped_existing: 49`
- `failed_exit_code: 0`
- `failed_missing_output: 0`

InterProScan evidence and SignalP6 localization evidence were produced for the
selected candidate set. These outputs remain machine-local unless a later
compact summary is explicitly promoted for Git tracking.

## Operational notes

The first launch attempt exposed a shell-environment issue and was restarted
with a preserved system `PATH`. That resolved the transient startup failure and
the batch then completed normally.

## Evidence boundary

- InterProScan evidence is domain/family/site architecture evidence.
- SignalP6 evidence is signal-peptide/localization support.
- These outputs are sequence and annotation evidence only. They do not prove
  PHB/PHA degradation phenotype on their own.
