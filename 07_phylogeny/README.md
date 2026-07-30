# P08 系统发育与 GTDB 映射

> 溯源契约更新：P08 对每条 P07 sequence manifest 记录校验显式 GTDB release，以及
> 标准化后的 P06 candidate table/scan manifest 路径；候选记录保留这些已验证来源及 SHA-256。
> P07 status 只能通过实际 FASTA shard stem 和相同 `input_fasta` 关联。pool core 永远从
> 17 个 authority seeds、15 个 cross-family challenges 和 5 个 close controls 构造；普通 P05
> 表中的 direct core seed/control 行一律失败关闭。

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
  P05 manifest 溯源；`source_checksum_kind` 必须显式为 `file_sha256` 或
  `residue_sha256`（仅具名历史 `p05_hmm_calibration_control_panel.tsv` 允许受限
  兼容适配）。核心 authority
  必须为 17 个种子、15 个非胞外挑战和 5 个近缘非靶标对照；不足即失败关闭。
- GTDB 候选必须逐行保留显式 GTDB release、P03 FAA 路径及 manifest/QC SHA-256、P06
  scan-manifest SHA-256、模型分数/coverage 阈值、P07 实际 terminal status、状态表路径和
  InterProScan/SignalP6 输出路径。P07 `family_categories` 必须包含当前 P06 family，重复状态键
  或来源不完整都会失败关闭。
- P06 candidate table 保留 HMMER 的 domain 行；P08 仅对被选择 tier 中同一
  `(family_category, proteome_shard, target_id)` 的多 domain 行作 target-level 汇总。汇总前必须
  核对 accession、target length、full-sequence score 与两个校准阈值完全一致；不一致即失败关闭。
  同一 target 同时出现 `High-confidence` 与 `Review` 时前者优先；同 tier 内选择最高
  `hmm_coverage` 的 domain，完全并列时以全行字典序决胜。候选 manifest 保留 domain 行数、观察到的
  selected tiers、选中 domain index 与选择规则，不能把 domain 数解释为独立候选或表型证据。
- P05 seed/control 的 `file_sha256` 先按原始字节严格核验；仅当 registry 声明 SHA 恰好等于同一
  实际文件的 LF 或 CRLF 规范表示时，才接受 Git 工作树换行转换。P08 同时记录声明 SHA、实际 raw-file
  SHA、canonical-LF SHA、验证模式与声明匹配表示；氨基酸残基、FASTA header、空白或任何非换行字节
  改变均继续失败关闭。`residue_sha256` 对照仍按残基核验，不使用该 EOL 适配器。
- P07 r8 manifest 中的 `fasta_shard`、`candidate_table_path`、`scan_manifest_path` 与 status
  `input_fasta`/`output_path` 可以是 r8 worktree 相对路径；必须用显式
  `--p07-source-root` 解析。P03 manifest/QC 中的相对 `faa_path` 必须由显式
  `--p03-source-root` 解析。P08 不改写这些原始声明，而是在 candidate manifest 中同时写入
  原始字面量、解析后路径和相应 SHA-256；未提供 source root 的相对路径将失败关闭。

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
1000`）；200–2,000 条用 `mafft --auto`；超过 2,000 条先形成可审计、确定性的代表序列契约，
命令 manifest 固定算法/version/参数、未来 mapping 路径、`not_materialized` SHA-256、记录数 0 和
“需单独授权”状态；本步骤不选择或物化代表序列。FastTree 不替代正式分析；
任何 IQ-TREE 子集与外群都需要独立批准。根定策略固定为：必须使用 accessioned
outgroup；没有已记录外群时只能生成 midpoint display，不得把展示根解释为祖先状态。

## T141 预检（尚未运行）

在已经带有 P08 脚本和所需 machine-local 输入的隔离工作树中，先生成计划，再只做
预检；下列命令是示例，不是本次任务执行记录。

```bash
WT=/home/data/haoyu/PHB-GTDB-GPT-p08-preflight-20260730
R8=/home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728
OLD=/home/data/haoyu/PHB-GTDB-GPT
cd "$WT"
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p08_prepare_phylogeny.py \
  --candidate-table "$R8/05_hmmer_scan/p06_hmmer_candidates.tsv" \
  --gtdb-release "GTDB Release 11 R232" \
  --p06-scan-manifest "$R8/05_hmmer_scan/p06_hmmer_scan_manifest.tsv" \
  --p03-prediction-manifest "$OLD/03_gtdb_proteomes/manifests/p03_prediction_manifest.tsv" \
  --p03-prediction-qc "$OLD/03_gtdb_proteomes/qc/p03_prediction_qc.tsv" \
  --p03-source-root "$OLD" \
  --p07-sequence-manifest "$R8/06_domain_annotation/manifests/p07_candidate_sequence_manifest.tsv" \
  --p07-status-table "$R8/06_domain_annotation/run_status/p07_domain_annotation_run_status.tsv" \
  --p07-source-root "$R8" \
  --model-registry "$WT/04_family_profiles/manifests/p05_hmm_model_registry.tsv" \
  --seed-registry "$WT/04_family_profiles/manifests/p05_hmm_seed_registry.tsv" \
  --control-panel "$WT/04_family_profiles/manifests/p05_hmm_calibration_control_panel.tsv" \
  --bac120-taxonomy "$OLD/00_raw_gtdb_r232/bac120_taxonomy_r232.tsv" \
  --ar53-taxonomy "$OLD/00_raw_gtdb_r232/ar53_taxonomy_r232.tsv" \
  --bac120-tree "$OLD/00_raw_gtdb_r232/bac120_r232.tree" \
  --ar53-tree "$OLD/00_raw_gtdb_r232/ar53_r232.tree" \
  --outdir "$WT/07_phylogeny" \
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
不是生物学阴性。预检中的 `skipped_existing` 专指“既有 route 输出非空”的完整性/恢复状态，
不表示已完成生物学分析；未来另行授权的运行才可能产生 `completed`，且仍不构成 PHB/PHA 表型证明。详细计划见
`docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md`。
