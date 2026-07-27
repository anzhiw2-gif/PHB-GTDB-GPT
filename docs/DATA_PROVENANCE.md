# Data Provenance

P01 is responsible for creating the first auditable input snapshot for this project.

## P01 implementation record

Tracked P01 source files:

- `scripts/p01_audit_gtdb.py`
- `tests/test_p01_audit_gtdb.py`

Representative T141 command:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
python scripts/p01_audit_gtdb.py \
  --paths config/paths.yaml \
  --source /home/data/haoyu/GTDB/gtdb_genomes_reps_r232/database \
  --target /home/data/haoyu/PHB-GTDB-GPT/00_raw_gtdb_r232/genomes \
  --threads 60
```

The audit performs a physical copy with `rsync -aHAX --info=progress2`, then validates the copy by:

- comparing source and target file counts
- comparing source and target total byte counts
- comparing per-top-level-directory byte counts
- verifying at least 1,000 sampled SHA256 hashes, which is also at least 1% of copied files
- copying the support files into `00_raw_gtdb_r232/`

The P01 tracked outputs on T141 are:

- `00_raw_gtdb_r232/genomes/`
- `00_raw_gtdb_r232/bac120_taxonomy_r232.tsv`
- `00_raw_gtdb_r232/ar53_taxonomy_r232.tsv`
- `00_raw_gtdb_r232/bac120_r232.tree`
- `00_raw_gtdb_r232/ar53_r232.tree`
- `00_raw_gtdb_r232/manifests/raw_genomes_manifest.tsv`

Observed completion state:

- raw genome manifest rows: `199,924`
- source and target reconciliation: passed
- sampled checksum verification: passed

P01 implementation history:

- `a0e1ccc` `feat: extend P01 audit tooling`
- `80c7930` `feat: make P01 audit executable`

Required provenance for the raw GTDB copy:

- Source release and source paths from `config/paths.yaml`
- Audit date in Asia/Shanghai time
- Source and target file counts
- Source and target byte counts
- At least 1% SHA256 verification, with a minimum of 1,000 copied files
- Tool versions for Python, Nextflow, Java, Slurm, Prodigal, Pyrodigal, HMMER, MAFFT, IQ-TREE 2, FastTree, InterProScan, SignalP, and Phobius

The Git repository must not contain the copied GTDB genomes or the full manifest for the raw tree copy. Those artifacts stay on T141 and are referenced through compact tracked reports and checksums.

The old project at `/home/data/haoyu/PHB_gtdb` is read only and may be used only as historical reference evidence after independent verification.

## P05 seed and model update, 2026-07-27

The approved P05 seed decision is recorded in
`docs/P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md`. The primary manifest is
`01_reference_library/manifests/reference_library.seed_manifest.tsv`; it
records accession, source database, version, organism, taxon ID, retrieval
date, evidence scope, FASTA path, and `profile_seed_status`.

New archaeal coverage sequences `CCQ32286.1`, `AGN01047.1`, `KYH27761.1`, and
`ELY43313.1` were retrieved from NCBI Protein E-utilities `efetch` on
2026-07-27. Exact-accession retrieval returned one FASTA record per accession;
local acceptance required the reported length and an N-terminal `GxSxG` motif.
These are E3 architecture-coverage sequences, not phenotype evidence. The
new bacterial `Q9WX79` sequence was retrieved from UniProtKB REST on
2026-07-27; it is an E2 *Paracoccus denitrificans* `phaZ` seed linked to
PMID `11267773` and DOI `10.1111/j.1574-6968.2001.tb10558.x`.

`archaeal_patatin_like_pha_dep`, `intracellular_mcl_pha_dep`, and
`intracellular_phaZ_no_lipase_box` were rebuilt in isolated T141 worktree
`/home/data/haoyu/PHB-GTDB-GPT-p05-rebuild-20260727` with MAFFT `v7.525`
L-INS-i (`--localpair --maxiterate 1000`) followed by HMMER `3.4`
`hmmbuild --amino`. The compact model, alignment, bundle, and sequence hashes
are in the tracked P05 registries. Raw profiles remain ignored and machine
local. All six P05 models are blocked from P06 pending calibration against
close non-target hydrolases; homology remains sequence evidence only.
