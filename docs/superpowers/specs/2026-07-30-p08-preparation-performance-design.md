# P08 preparation performance design

## Purpose

Make planning-only P08 preparation practical for the accepted full P06/P07
machine-local inputs without relaxing provenance checks. This does not
authorize alignment, tree inference, rooting, or biological interpretation.

## Observed cause

The candidate manifest calculates SHA-256 repeatedly for the same P06 ledger,
P06 scan manifest, P03 manifests, and P07 status table while constructing
every candidate row. The stopped full T141 run was CPU-active for nearly two
hours and had issued roughly 1.68 TB of logical reads without materializing
output.

## Approved design

1. An invocation-scoped input-digest cache binds each canonical path to its
   SHA-256 and file identity (`st_dev`, `st_ino`, `st_size`, `st_mtime_ns`). A
   changed identity is an error, not a cache refresh.
2. `prepare_p08_inputs(..., workers=1)` and CLI `--workers` accept only 1–60.
3. Up to `min(workers, distinct P07 FASTA shards)` threads preload independent
   P07 FASTA records and digests. The main thread consumes sorted results and
   is the only artifact writer, preserving deterministic outputs.
4. Requested/effective workers are recorded in compact P08 input provenance.

## Boundaries and verification

The candidate table remains fully validated. Inputs remain machine-local and
must be immutable during one invocation. New tests cover cache reuse, mutation
failure, invalid worker counts, deterministic one-vs-sixty-worker candidate
rows, and worker provenance. Existing P08 tests and repository validation must
pass before a fresh T141 output directory is launched with `--workers 60`.
