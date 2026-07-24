# Data Provenance

P01 is responsible for creating the first auditable input snapshot for this project.

Required provenance for the raw GTDB copy:

- Source release and source paths from `config/paths.yaml`
- Audit date in Asia/Shanghai time
- Source and target file counts
- Source and target byte counts
- At least 1% SHA256 verification, with a minimum of 1,000 copied files
- Tool versions for Python, Nextflow, Java, Slurm, Prodigal, Pyrodigal, HMMER, MAFFT, IQ-TREE 2, FastTree, InterProScan, SignalP, and Phobius

The Git repository must not contain the copied GTDB genomes or the full manifest for the raw tree copy. Those artifacts stay on T141 and are referenced through compact tracked reports and checksums.

The old project at `/home/data/haoyu/PHB_gtdb` is read only and may be used only as historical reference evidence after independent verification.
