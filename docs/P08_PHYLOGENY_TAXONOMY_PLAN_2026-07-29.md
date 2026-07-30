# P08 系统发育与 GTDB 分类映射计划

## 最终溯源失败关闭补充

P08 唯一接受 P07 默认和文档使用的精确 release 字面量 `GTDB Release 11 R232`；不接受
缩写、别名或大小写漂移。它在读取 P07 sequence manifest 时逐行核对 `gtdb_release`，并以标准化路径核对
`candidate_table_path` 和 `scan_manifest_path` 是否分别等于显式 P06 输入；验证后的两条
P06 来源及 SHA-256 写回候选记录。P07 status 仅以实际 FASTA shard stem 关联，且每个工具的
`input_fasta` 必须解析为同一实际 shard，避免同 stem 的异路径复用。

P07 r8 manifest 使用相对于其原始 r8 worktree 的路径。P08 因而提供显式
`--p07-source-root`：它只对 P07 manifest/status 中的相对声明生效；绝对声明保持原义。候选
manifest 同时记录原始声明、解析后的实体路径和 SHA-256，避免以新 P08 worktree 的当前目录
错误解释旧 P07 产物。任一相对声明未给 source root、解析后不等于显式 P06 输入或实体不可读，均失败关闭。

`extracellular_pha_depolymerase_core` 禁止普通 P05 seed/control 表中的 direct core 行。每次
都由 tracked authority 重新构造并验证 17 个 core seeds、15 个活跃非胞外 cross-family
challenges 和 5 个 accessioned close controls（共 20 个 hard-panel records）。这仍是序列与
边界对照证据，不构成 PHB/PHA degradation phenotype 或 subtype 结论。

**日期：** 2026-07-29，Asia/Shanghai
**阶段：** P08，家族系统发育输入、GTDB 分类映射和仅预检命令计划
**当前状态：** 本地 scaffold/preflight 实现完成；T141 preflight 尚未运行；未执行
MAFFT、FastTree 或 IQ-TREE，未形成实际多序列比对或推断树。

## 生物学目的与证据边界

P08 在每个已批准 P05 模型家族中，保留 GTDB R232 候选、accessioned 种子和硬对照，
再把 P07 的序列/独立注释状态与 GTDB Bac120/Ar53 分类映射连接起来。该组织方式可
帮助审查序列相似性、家族内分布、域/定位证据与分类背景是否一致，并使后续树的输入
可追溯。

但是，家族树、树上位置、GTDB 分类、同源性、域架构和定位都只是序列或注释证据。
它们不能单独证明任何 GTDB 基因组、蛋白或家族具有 PHB/PHA 降解表型；运行/完整性
失败同样不能视为生物学阴性。

## 输入准入与溯源

1. P05 模型 registry 只允许 `approved_for_p06=yes` 和
   `scan_permission=approved` 的四个模型。P08 不恢复被阻断的胞外亚型模型，也不把
   `extracellular_pha_depolymerase_core` 拆成互斥亚型。
2. P06 默认仅选 `High-confidence`。需要检查 `Review` 时，必须显式添加
   `--include-tier Review`；`Rejected` 不是 P08 输入。多个批准模型命中的 target 保持
   非互斥记录，不能仅凭树或命中数强行定为单一功能类别。
3. P05 种子和对照都必须由 manifest 中的 accession、source path、模型 SHA-256 和
   序列 SHA-256 验证。种子是已有的实验/注释范围证据，对照用于邻近水解酶边界检查；
   两者都不替代 GTDB 候选的表型实验。
4. P08 必须显式传入 GTDB release、P06 scan manifest、P03 prediction manifest/QC；每个候选
   逐行保留 FAA 源路径、upstream 路径与 SHA-256，以及 P06 模型分数、coverage、阈值和 tier。
   P03 accession/FAA 路径与 P07 源路径不一致即关闭处理，绝不推断 release。
5. P07 只接受已完成或 `skipped_existing` 的 InterProScan 与 SignalP6 状态，并核对
   P06/P07 target、family_categories、长度、FASTA 和 shard 关联；候选记录保留实际 terminal
   status、P07 状态表路径及每个工具的 output path。重复 `(tool, fasta_shard)` 状态键、状态缺失、失败
   或校验和冲突必须关闭处理并写入 review 记录。
6. Bac120 taxonomy 仅用于细菌 assembly 映射，Ar53 taxonomy 仅用于古菌 assembly
   映射；CLI 以两个独立参数加载，分别验证 `d__Bacteria`/`d__Archaea`，并拒绝
   标准化 accession 跨来源重叠。taxonomy join 逐行记录来源角色、路径和 SHA-256。
   Bac120/Ar53 参考树只作为现有的 provenance/preflight 输入：必须为非空文件，且
   路径、角色和 SHA-256 写入 `p08_input_provenance.tsv`；CLI 不读取其拓扑来推断、
   重根或生成任何新树。

## 本地 scaffold 与命令计划

`scripts/p08_prepare_phylogeny.py` 写出候选、taxonomy join、reference、family input、
summary、phylogeny command 和 hash-locked input provenance manifests，并写出每个家族
的 machine-local 输入 FASTA。所有命令状态在该步固定为 `planned_not_run`；脚本不调用
`subprocess.run`，不运行 MAFFT、FastTree 或 IQ-TREE。

计划的 MAFFT 路线按照候选数量（不含种子/对照）确定：

- `<200`：L-INS-i 计划，`--localpair --maxiterate 1000 --inputorder`；
- `200–2000`：`mafft --auto --inputorder` 计划；
- `>2000`：先写入确定性代表序列契约，不物化子集；固定算法/version/参数、未来 mapping 路径、
  `not_materialized` mapping SHA-256、记录数 0 与单独授权边界。将 FastTree 仅保留为后续探索性
  选择，代表输入、mapping 和 SHA-256 必须在未来独立授权后才可物化。

FastTree 用于大型家族的探索性概览，不能替代系统的统计推断。IQ-TREE 仅可在另行批准
的子集、模型选择和 accessioned outgroup 已确定后使用；现有 `-m TEST -B 1000` 模板
只是未执行的待审批记录。根定策略精确固定为
`explicit_accessioned_outgroup_required; otherwise midpoint_display_only`：无已 accessioned
外群时只能做 midpoint 展示根，不能报告祖先方向、获得/丢失事件或功能演化结论。

## 工具命令来源与版本记录

本地脚本复用 Tasks 1–4 已有的计划命令模式，没有复制第三方实现，也没有在本次任务
执行任何工具。后续 T141 preflight 前应重新核对以下官方命令参考，并将实际 `--version`
输出、获取日期和可执行文件路径写入 machine-local status/compact summary：

- MAFFT manual：<https://mafft.cbrc.jp/alignment/software/manual/manual.html>；
- IQ-TREE 2 command reference：<https://iqtree.github.io/doc/Command-Reference>；
- FastTree 官方页：<http://www.microbesonline.org/fasttree/>。

本地对上述网页的读取尝试在当前环境超时，因此本计划不把任何在线页面内容当作已验证
的版本事实；T141 的实际工具版本和可用选项仍是 preflight 接受条件。

## T141 运行顺序（未执行）

在隔离工作树、输入 checksum 和工具版本确认后，先准备命令 manifest；随后只运行
`scripts/p08_run_phylogeny.py --preflight-only`。该 runner CLI 保持严格为
`--manifest --status-dir --workers --preflight-only`，缺少 `--preflight-only` 会直接拒绝。
本任务没有执行下列命令：

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p08_run_phylogeny.py \
  --manifest 07_phylogeny/manifests/p08_phylogeny_command_manifest.tsv \
  --status-dir 07_phylogeny/run_status \
  --workers 1 \
  --preflight-only
```

预检只有在输入 FASTA、whole-file SHA-256、命令解析和可执行文件可用时才给出
`preflight_ok`，并且仍不执行命令。status 逐行绑定 command-manifest 全文件 SHA-256、解析的
可执行文件路径和 `not_queried_preflight_only` 版本状态；不调用工具取得版本。`planned_not_run`
表示准备阶段写出的原始计划；`missing_executable`、`missing_input`、`checksum_mismatch`、
`failed_exit_code` 与 `failed_missing_output` 都是故障/完整性状态。预检的 `skipped_existing`
仅指既有 route 输出非空的完整性/恢复状态，不表示已完成生物学分析；未来经独立授权才可出现
`completed`，且这些执行状态不改变生物学证据边界。

## 接受标准与产物位置

本地 scaffold 的接受条件是：CLI 以 High-confidence 默认值运行时返回 0；显式 Review
选择被保留；Bac120/Ar53 分类映射均可追溯；每个 family input/command manifest 行保留
来源与 SHA-256；所有 command status 均为 `planned_not_run`；prepare 步不执行外部
工具；runner 仅允许 `--preflight-only`。

`07_phylogeny/family_fastas/`、`alignments/`、`trees/`、`gtdb_mapping/`、`annotations/`、
`manifests/`、`review/`、`run_status/` 和 `run_logs/` 均为 machine-local 生成物并被
Git 忽略。若未来需要把紧凑摘要纳入 Git，必须同时记录输入 manifest、P05 模型和种子
checksum、GTDB R232 版本、工具版本、命令、日期、接受指标和上述证据边界。
