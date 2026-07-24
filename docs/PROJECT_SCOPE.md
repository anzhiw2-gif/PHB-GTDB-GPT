# Project Scope

This repository controls a new PHB/PHA degradation-gene analysis of the complete, unfiltered GTDB R232 genome collection on T141.

It will use a conservative combined reference framework: DED family definitions, TIGRFAM models including TIGR01849 and TIGR02240, and accessioned experimentally supported sequence seeds. Family HMM hits will be reviewed with domain architecture, localization evidence, family trees, and GTDB taxonomy.

The repository stores source code, configurations, tests, documentation, manifests, and small curated reference inputs. It never stores raw GTDB genomes, full predicted proteomes, bulk HMMER/InterPro outputs, or other large generated data.

The legacy project at `/home/data/haoyu/PHB_gtdb` is read only. Its candidate sequences can serve as historical/reference evidence after independent verification, but its output is not a complete GTDB proteome set and is not the comparison baseline.
