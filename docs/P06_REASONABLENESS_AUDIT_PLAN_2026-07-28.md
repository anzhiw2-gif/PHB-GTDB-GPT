# P06 candidate reasonableness audit plan

**Date:** 2026-07-28, Asia/Shanghai  
**Stage:** P06 audit before P07 annotation  
**Status:** compact audit completed on T141; generated outputs remain
machine-local

## Purpose

This audit converts the large machine-local P06 candidate table into compact
scale checks before P07 interpretation:

- row counts by family and tier;
- de-duplicated target counts by family and tier;
- High-confidence candidate rates against the P03 total predicted-gene and
  genome denominators;
- targets that pass more than one P06 model at High-confidence level.

The audit measures candidate-table reasonableness only. HMM hits are sequence
classification evidence and are not phenotype proof.

## Source review before implementation

Before adding project-specific code, public documentation and code patterns
were checked for reusable alternatives:

- SeqKit official documentation: `https://bioinf.shenwei.me/seqkit/`
  - useful for FASTA statistics and filtering, but not sufficient for the
    project-specific combination of P06 tier, family, de-duplicated target, and
    interpretation-boundary summaries.
- csvkit official documentation: `https://csvkit.readthedocs.io/`
  - useful for generic delimited-table inspection, but not enough for stable
    P06 family/tier overlap audit without custom joins.
- GTDB R232 release information:
  `https://gtdb.ecogenomic.org/stats/r232`
  - denominator reference for the GTDB R232 representative-genome scale.

No third-party code was copied. The committed helper uses Python standard
library only because T141 project rules prefer dependency-light, auditable
stage glue code.

## Prepared script

```bash
scripts/p06_candidate_reasonableness.py
```

Default T141 command:

```bash
cd /home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p06_candidate_reasonableness.py \
  --candidate-table 05_hmmer_scan/p06_hmmer_candidates.tsv \
  --outdir 06_domain_annotation/p06_reasonableness \
  --total-predicted-genes 615969593 \
  --total-genomes 199923
```

Expected compact outputs:

- `06_domain_annotation/p06_reasonableness/p06_candidate_reasonableness_summary.tsv`
- `06_domain_annotation/p06_reasonableness/p06_family_tier_reasonableness.tsv`
- `06_domain_annotation/p06_reasonableness/p06_high_confidence_overlap_targets.tsv`
- `06_domain_annotation/p06_reasonableness/P06_REASONABLENESS_AUDIT.md`

These generated audit outputs stay machine-local by default. A later P06/P07
completion record may promote a compact accepted summary into Git after review.

## T141 audit result

The compact audit was run in isolated T141 worktree
`/home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728` against the machine-local P06
candidate table. It reported:

- total candidate rows: `3,080,953`
- unique candidate targets: `2,330,939`
- High-confidence rows: `37,915`
- High-confidence unique targets: `37,912`
- High-confidence overlap targets: `0`
- High-confidence rows per predicted gene: `0.000062`
- High-confidence unique targets per genome: `0.189633`

High-confidence unique targets by family:

- `archaeal_patatin_like_pha_dep`: `678`
- `extracellular_pha_depolymerase_core`: `10,453`
- `intracellular_mcl_pha_dep`: `349`
- `intracellular_phaZ_no_lipase_box`: `26,432`

The audit found no model-overlap conflict among High-confidence targets, but
the large no-lipase-box and extracellular pools still require P07 architecture
and localization review before interpretation.

## Biological interpretation

The literature-based expectation is that PHA depolymerase-related candidates
are heterogeneous across extracellular, intracellular, no-lipase-box, and
archaeal patatin-like routes. Therefore, a reasonable P06 audit should not
expect a single uniform family count. Large no-lipase-box or broad extracellular
candidate pools can be plausible, but they must be reduced by P07 domain
architecture/localization and P08 phylogeny/taxonomy evidence before reporting.

High-confidence overlap between models is not automatically a false positive:
it is a review flag for architecture and tree placement. Conversely, absence
from a model or exclusion by HMMER length limits is not biological absence.
