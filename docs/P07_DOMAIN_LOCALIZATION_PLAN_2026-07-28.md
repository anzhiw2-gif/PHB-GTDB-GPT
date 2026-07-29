# P07 domain-architecture and localization plan

**Date:** 2026-07-28, Asia/Shanghai  
**Stage:** P07, independent annotation review for P06 GTDB R232 sequence candidates  
**Status:** candidate extraction complete on T141; annotation execution
completed on 2026-07-29 after resolving the local toolchain

## Entry condition

P06 is technically accepted for the four checksum-locked HMMER models. The P07
input is the machine-local P06 candidate table and scan manifest from isolated
T141 worktree `/home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728`.

Default P07 input is restricted to P06 `High-confidence` rows. `Review` rows may
be added only as an explicit secondary review pass. `Rejected` rows are not P07
annotation inputs unless a later documented audit requires a narrow rescue set.

## Biological rationale

P07 measures independent protein-level evidence that was not part of the P06
HMM hit threshold itself:

- domain/family/site architecture, to distinguish coherent PhaZ/PhaZh-like
  candidates from near-neighbor hydrolases or conflicting architectures;
- signal-peptide and membrane-topology evidence, to separate extracellular,
  cell-associated, and intracellular interpretation routes;
- sequence extraction provenance, to preserve the link from each P07 protein
  back to the P06 hit, P03 proteome file, GTDB R232 shard, and source candidate
  manifest.

These signals remain sequence and annotation evidence. They do not prove an
experimental PHB/PHA degradation phenotype.

## Tool sources checked before scaffolding

- InterProScan official documentation:
  `https://interproscan-docs.readthedocs.io/`
  - appropriate because P07 needs reproducible standalone InterProScan command
    patterns and multi-format output for protein signature/domain review.
- SignalP 6.0 official service/software page:
  `https://services.healthtech.dtu.dk/services/SignalP-6.0/`
  - appropriate because signal-peptide prediction is a direct localization
    support signal for bacterial and archaeal candidate interpretation.
- Phobius official page:
  `https://phobius.sbc.su.se/`
  - appropriate as an optional independent signal-peptide/transmembrane-topology
    check, but local CLI availability and license/access must be verified on
    T141 before running it.
- nf-core modules InterProScan module:
  `https://github.com/nf-core/modules/tree/master/modules/nf-core/interproscan`
  - appropriate as a maintained workflow-module reference if this stage is later
    promoted from Python command manifests into Nextflow DSL2.

No external code was copied into the repository. The committed P07 script is
custom glue code because the immediate need is smaller: extract P06 candidate
sequences with provenance and prepare planned-not-run annotation commands.

## Prepared local scaffold

Preparation script:

```bash
scripts/p07_prepare_domain_annotation.py
```

Resumable runner script:

```bash
scripts/p07_run_domain_annotation.py
```

Default T141 command:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p07_prepare_domain_annotation.py \
  --candidate-table 05_hmmer_scan/p06_hmmer_candidates.tsv \
  --scan-manifest 05_hmmer_scan/p06_hmmer_scan_manifest.tsv \
  --outdir 06_domain_annotation \
  --include-tier High-confidence \
  --sequences-per-shard 500
```

Expected outputs:

- `06_domain_annotation/input/fasta_shards/p07_candidates_*.faa`
- `06_domain_annotation/manifests/p07_candidate_sequence_manifest.tsv`
- `06_domain_annotation/manifests/p07_domain_annotation_command_manifest.tsv`
- `06_domain_annotation/manifests/p07_domain_annotation_summary.tsv`
- if any selected sequence is missing:
  `06_domain_annotation/review/p07_missing_candidate_sequences.tsv`

The script fails closed when a selected P06 target is absent from the P06
manifest proteome path or when P06 target length disagrees with the extracted
FASTA sequence. Missing sequences should be treated as provenance/integrity
blockers, not as biological negatives.

## Planned command policy

The command manifest records commands as `planned_not_run`.

- InterProScan command rows write TSV, JSON, and GFF3 outputs and include GO and
  pathway annotations when the local InterProScan installation supports them.
- SignalP6 command rows use `--organism other` as the cross-domain default for
  bacterial and archaeal candidates; the exact installed CLI should be checked
  with `signalp6 --help` on T141 before execution.
- Phobius command rows are emitted only when `--phobius-exe` is supplied,
  because local standalone availability is not guaranteed.

Raw annotation outputs stay ignored by Git. A later P07 completion record should
commit only compact summaries, tool versions, commands, checksums, acceptance
metrics, and interpretation boundaries.

## T141 preflight result

The High-confidence P07 input extraction was run in isolated T141 worktree
`/home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728` and produced:

- selected unique candidate sequences: `37,912`
- extracted sequences: `37,912`
- missing sequences: `0`
- FASTA shards: `76`
- planned InterProScan command rows: `76`
- planned SignalP6 command rows: `76`

The P07 runner preflight was then executed with:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p07_run_domain_annotation.py \
  --manifest 06_domain_annotation/manifests/p07_domain_annotation_command_manifest.tsv \
  --status-dir 06_domain_annotation/run_status \
  --workers 1 \
  --preflight-only
```

Actual execution initially used the local script piped over SSH because the
isolated T141 worktree did not yet contain the new local script. The preflight
wrote `06_domain_annotation/run_status/p07_domain_annotation_run_status.tsv`
and reported:

- `preflight_ok: 0`
- `missing_executable: 152`
- `completed: 0`
- `skipped_existing: 0`
- `failed_exit_code: 0`
- `failed_missing_output: 0`

`interproscan.sh`, `signalp6`, and `java` were not found in the checked T141
environment or common local paths. That startup issue was resolved on T141
with a preserved system `PATH`, after which the full P07 batch completed
normally.

## T141 execution result

- selected unique candidate sequences: `37,912`
- FASTA shards: `76`
- runner workers: `6`
- InterProScan CPU per job: `10`
- final run-status rows: `152`
- `completed`: `103`
- `skipped_existing`: `49`
- `failed_exit_code`: `0`
- `failed_missing_output`: `0`

The completion record is now `docs/P07_COMPLETION_2026-07-29.md`.
