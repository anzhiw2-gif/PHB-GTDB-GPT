# P08 系统发育与 GTDB 映射

P08 将已接受的 P06 候选、P07 序列/注释状态、P05 参考种子和对照、以及
GTDB R232 Bac120/Ar53 分类信息组织成可审计的家族输入和系统发育计划。家族树、
GTDB 映射、域架构和定位信息均是序列/注释证据；它们不单独证明 PHB/PHA
降解表型。

## 当前状态与输入边界

- 本地 scaffold 和仅预检 CLI 已完成；尚未在 T141 执行 P08 preflight，未执行
  MAFFT、FastTree 或 IQ-TREE，也没有实际比对或树推断结果。
- 仅接受 P05 `approved_for_p06=yes` 且 `scan_permission=approved` 的四个模型。
- P06 默认只纳入 `High-confidence`；`Review` 必须在 CLI 中显式重复传入
  `--include-tier Review`。`Rejected` 不能作为 P08 输入。
- 同一个 target 可保留多个批准模型的命中；P08 不把它们强制改写为互斥家族或
  表型分类。
- 种子和硬对照必须保留 accession、来源路径、模型 SHA-256、序列 SHA-256 和
  P05 manifest 溯源；GTDB 候选必须保留 P06 shard/target、P07 sequence ID 和
  assembly/taxonomy 映射。

## 生成物与执行边界

`scripts/p08_prepare_phylogeny.py` 只写 machine-local 的 FASTA、输入 manifest、
GTDB 映射和 `planned_not_run` 命令 manifest。Bac120/Ar53 树是 provenance/
preflight 输入：CLI 检查它们为非空本地文件，并在
`manifests/p08_input_provenance.tsv` 中固化 taxonomy/tree 的角色、路径和 SHA-256；
绝不读取树拓扑、请求或生成推断树。Bac120/Ar53 taxonomy 也分别核对
`d__Bacteria`/`d__Archaea`，拒绝跨来源 accession 重叠。

所有 `07_phylogeny/family_fastas/`、`alignments/`、`trees/`、`gtdb_mapping/`、
`manifests/`、`review/`、`run_status/` 和 `run_logs/` 生成物均保持 machine-local
并由 `.gitignore` 排除；Git 只保存脚本、测试和紧凑文档。

候选数量决定计划路线：少于 200 条用 MAFFT L-INS-i（`--localpair --maxiterate
1000`）；200–2,000 条用 `mafft --auto`；超过 2,000 条先形成可审计、确定性的
代表序列计划，之后才可把 FastTree 作为探索性路线。FastTree 不替代正式分析；
任何 IQ-TREE 子集与外群都需要独立批准。根定策略固定为：必须使用 accessioned
outgroup；没有已记录外群时只能生成 midpoint display，不得把展示根解释为祖先状态。

## T141 预检（尚未运行）

在已经带有 P08 脚本和所需 machine-local 输入的隔离工作树中，先生成计划，再只做
预检；下列命令是示例，不是本次任务执行记录。

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p08_prepare_phylogeny.py \
  --candidate-table 05_hmmer_scan/p06_hmmer_candidates.tsv \
  --p07-sequence-manifest 06_domain_annotation/manifests/p07_candidate_sequence_manifest.tsv \
  --p07-status-table 06_domain_annotation/run_status/p07_domain_annotation_run_status.tsv \
  --model-registry 04_family_profiles/manifests/p05_hmm_model_registry.tsv \
  --seed-registry 04_family_profiles/manifests/p05_hmm_seed_registry.tsv \
  --control-panel 04_family_profiles/manifests/p05_hmm_calibration_control_panel.tsv \
  --bac120-taxonomy /path/to/bac120_taxonomy.tsv \
  --ar53-taxonomy /path/to/ar53_taxonomy.tsv \
  --bac120-tree /path/to/bac120.tree \
  --ar53-tree /path/to/ar53.tree \
  --outdir 07_phylogeny \
  --include-tier High-confidence

/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p08_run_phylogeny.py \
  --manifest 07_phylogeny/manifests/p08_phylogeny_command_manifest.tsv \
  --status-dir 07_phylogeny/run_status \
  --workers 1 \
  --preflight-only
```

`planned_not_run` 表示计划已写入但从未执行；`preflight_ok` 只表示输入、SHA-256、
命令解析和可执行文件检查通过；`missing_executable`、`missing_input`、
`checksum_mismatch`、`failed_exit_code` 和 `failed_missing_output` 是完整性/运行状态，
不是生物学阴性。`completed` 或 `skipped_existing` 只有在未来另行授权的运行中才可能
出现，且仍不构成 PHB/PHA 表型证明。详细计划见
`docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md`。
