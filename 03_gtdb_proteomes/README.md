# P03 GTDB Proteomes

Complete predicted proteomes, GFF files, and per-genome QC are generated on
T141 and ignored by Git.

The original file-count-complete P03 run was not accepted because its FAA
records contained empty protein translations. The translation-fix rerun was
accepted before P06: 199,923 FAA/QC/manifest rows have populated protein
metrics and P06 used those repaired proteomes.

For a later provenance check, inspect the accepted rerun with:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_translation_fix.py --once
```

P06 completed its four-model GTDB-wide scan and candidate parsing on 2026-07-28.
The current scan acceptance record is `docs/P06_COMPLETION_2026-07-28.md`.
