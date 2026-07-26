# P03 GTDB Proteomes

Complete predicted proteomes, GFF files, and per-genome QC will be generated here on T141 and ignored by Git.

As of 2026-07-26, the original P03 file-count-complete run is not accepted for
downstream HMMER scanning because its FAA records contained empty protein
translations. A translation-fix rerun is active on T141.

Monitor the repair run with:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_translation_fix.py --once
```

P06 must wait for the repaired FAA files, regenerated QC, and regenerated
manifest before GTDB-wide scanning.
