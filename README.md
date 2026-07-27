# PHB-GTDB-GPT

这是一个面向 GTDB Release 11 / R232 的 PHB/PHA depolymerase 相关基因分析项目。

项目主线是：
用经过实验支持的参考序列、DED family definitions、TIGRFAM profiles、自建 HMM、domain architecture review 和 family-level phylogenetic analysis，
去做一套可追溯、可复现、分家族的系统分析。

这里检测到的同源蛋白只被当作 sequence evidence，不会被直接当成“已经证明具备 PHB 降解表型”。

## 当前阶段

| 阶段 | 说明 |
| --- | --- |
| P01 | 在 T141 上审计并物理复制未过滤的 GTDB R232 输入。 |
| P02 | 用确定性的 GTDB 子集锁定唯一生产预测路线 `Pyrodigal GeneFinder(meta=True)`。 |
| P03 | 预测全量 GTDB proteomes，并做基础 QC。 |
| P04 | 构建可审计的实验参考库和 family definitions。 |
| P05 | 构建并校准 family HMMs。 |
| P06 | 用 HMMER 扫描 GTDB proteomes，并按 family 规则筛选。 |
| P07 | 做 InterPro domain architecture 和 localization 复核。 |
| P08 | 构建 family trees，并结合 GTDB taxonomy / species tree 分析。 |
| P09 | 输出表格、图、报告和 provenance manifests。 |
| P10 | 可选：在合适的 SBML/JSON model 上做 COBRApy 验证。 |

## 目录说明

主要阶段性文档在 `docs/`：

- [docs/DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md) 记录 P01 的审计、复制和 manifest 轨迹
- [docs/P02_BENCHMARK_DECISION.md](docs/P02_BENCHMARK_DECISION.md) 记录 P02 的 benchmark 结果和预测器锁定
- [docs/PREDICTION_POLICY.md](docs/PREDICTION_POLICY.md) 记录 P03 使用的生产预测策略
- [docs/HANDOFF_2026-07-26_P05_P06_P03_TRANSLATION_FIX.md](docs/HANDOFF_2026-07-26_P05_P06_P03_TRANSLATION_FIX.md) 记录 P05/P06 当前状态和 P03 翻译修复重跑
- [docs/P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md](docs/P05_HMM_SEED_SELECTION_DECISION_2026-07-27.md) 固化细菌/古菌 HMM 种子证据规则和家族优化决定
- [docs/P05_HMM_MODEL_CATALOG.md](docs/P05_HMM_MODEL_CATALOG.md) 发布 GitHub 可追溯的 HMM 元数据、种子登记和 P06 校验门控
- [docs/P05_EXTRACELLULAR_SUBTYPE_ARCHIVE_RECORD_2026-07-28.md](docs/P05_EXTRACELLULAR_SUBTYPE_ARCHIVE_RECORD_2026-07-28.md) 记录三个胞外子模型、校准边界和 pooled core 的 P06 交接

主要脚本在 `scripts/`：

- `p01_audit_gtdb.py`
- `p02_select_benchmark_genomes.py`
- `p02_compare_predictors.py`
- `p03_predict_proteomes.py`
- `p03_monitor_progress.py`
- `p03_monitor_translation_fix.py`
- `p06_scan_family_profiles.py`

## 运行环境

计算密集型任务都在 T141 上执行，项目目录是：

`/home/data/haoyu/PHB-GTDB-GPT`

服务器配置模板见 [config/paths.example.yaml](config/paths.example.yaml)。

Git 不跟踪这些大文件或运行产物：

- 原始 GTDB 文件
- 预测 proteomes
- HMMER / InterPro 结果
- Nextflow work 目录
- 各类生成报告

这些内容会通过 manifest 和 QC 文件保留可追溯性。

## 快速验证

本地常用检查：

```powershell
python -m unittest tests/test_repository_layout.py -v
python scripts/validate_repository.py
git diff --check
```

在 T141 上查看 P03 进度：

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_progress.py
```

只看一眼当前状态：

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_progress.py --once --no-clear
```

查看当前 P03 翻译修复重跑：

```bash
cd /home/data/haoyu/PHB-GTDB-GPT
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python scripts/p03_monitor_translation_fix.py --once
```

当前护栏：P03 翻译修复已完成并通过 QC；但 P06 的新扫描仍被冻结，直到确认优化种子、重建或明确保留模型、记录并校验 HMM/bundle/alignment SHA256、完成校准阈值决定，并由 model registry 将模型标记为 `approved_for_p06=yes`。原始 domtblout 和候选表仍要等 compact summaries 审核后再纳入 Git。

## 运行前提

Nextflow 工作流需要 T141 上具备 Java 17+ 和 Nextflow 24.10.0+。

如果你要先看阶段决策，建议从这两个文件开始：

1. [docs/DATA_PROVENANCE.md](docs/DATA_PROVENANCE.md)
2. [docs/P02_BENCHMARK_DECISION.md](docs/P02_BENCHMARK_DECISION.md)
