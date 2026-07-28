# P06 四模型 GTDB 扫描完成与解析验收

**完成日期：** 2026-07-28（Asia/Shanghai）  
**阶段：** P06，GTDB Release 11 R232 P03 proteomes 的 HMMER 序列候选筛查  
**结论：** 四个已批准模型的原始扫描和候选解析均已技术验收；P07 的结构域、定位、系统发育与分类学复核尚未启动。

## 可复现执行记录

- 主机与隔离工作树：T141，`/home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728`
- 源代码修订：`80991a7`；启动记录：`docs/P06_LAUNCH_2026-07-28.md`
- Python：`/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python`；HMMER：3.4
- 输入：199,923 份已通过 P03 translation/QC 验收的 GTDB R232 `.faa.gz` proteomes
- 调度：每个任务 200 份 proteomes；每个模型 1,000 个 chunk；4 个模型共 4,000 个 HMMER 任务；16 个可恢复 workers；每个 `hmmsearch` 使用 `--cpu 1`。
- 模型、阈值、SHA256 与解释边界固定在启动记录及 `docs/P05_HMM_MODEL_CATALOG.md`。其中 `extracellular_pha_depolymerase_core` 使用批准的 r6 工件（SHA256 `74c4b69a2d845f0725d0bc348402e6a51ba3c17a9f67f8cabed3b63df6a6e2f4`）。

## 原始扫描验收

恢复日志 `05_hmmer_scan/run_status/p06_full_run_r8_resume.log` 记录：

| 指标 | 结果 |
|---|---:|
| `completed` | 3,760 |
| `skipped_existing` | 240 |
| `failed_exit_code` | 0 |
| `failed_empty_domtblout` | 0 |
| 已核对的状态行 | 4,000（3,760 + 240） |
| 非空原始 `domtblout` | 4,000 |

初次运行遇到 HMMER 3.4 对单条目标长度超过 100,000 aa 的工具限制。r8 的流式输入修复会将这些目标从 HMMER 输入中排除，并在 `overlong_protein_exclusions/` 中以 `hmmsearch_target_length_gt_100000` 逐条审计。该标记不是 HMMER 阴性结果，也不能解释为生物学缺失。

## 解析验收与候选分层

在全部 4,000 份原始输出存在且无失败状态后，执行：

```bash
cd /home/data/haoyu/PHB-GTDB-GPT-p06-r8-20260728
/home/data/haoyu/miniconda3/envs/phb_gtdb/bin/python \
  scripts/p06_scan_family_profiles.py \
  --parse-only \
  --parse-manifest 05_hmmer_scan/p06_hmmer_scan_manifest.tsv \
  --outdir 05_hmmer_scan
```

`p06_hmmer_candidate_summary.tsv` 的验收项为 `families=4`、`proteome_shards=1000`、`missing_domtblout=0`。解析产生 3,080,953 条域级 HMMER 记录；它们不是去重后的蛋白数，也不是表型阳性数。

| 模型 | High-confidence | Review | Rejected | P06 可用解释 |
|---|---:|---:|---:|---|
| `archaeal_patatin_like_pha_dep` | 678 | 44,972 | 265,791 | 古菌 patatin-like 序列候选 |
| `intracellular_mcl_pha_dep` | 349 | 258,555 | 1,352,913 | 胞内 mcl-PHA depolymerase 序列候选 |
| `intracellular_phaZ_no_lipase_box` | 26,432 | 5,601 | 46,940 | 无典型 lipase box 的胞内 PhaZ 序列候选 |
| `extracellular_pha_depolymerase_core` | 10,456 | 23,409 | 1,044,857 | 胞外 PHA depolymerase core 序列候选；不能判定 mcl/scl 或 type-I/type-II |
| **总计** | **37,915** | **332,537** | **2,710,501** | **3,080,953 条域级记录** |

`High-confidence` 表示该域记录通过对应模型登记的 calibrated score、HMM coverage 及审查门槛；它不等价于实验验证的 PHB/PHA 降解表型。细菌参考序列的实验支持要求和古菌 E3 架构覆盖限制仍只用于模型构建与解释边界，不能由本轮同源筛查消除。

## 数据保留与下一关

以下工件保持 machine-local 且由 `.gitignore` 排除：原始 `domtblout`、全量 `p06_hmmer_candidates.tsv`、运行状态表、日志、以及超长目标审计表。Git 仅记录本完成报告、启动/模型来源、代码版本、命令和紧凑验收指标。

后续只有在 P07 对高置信和必要的 Review 候选完成独立的结构域架构、分泌/定位、近缘非目标水解酶排除、系统发育和 GTDB taxonomy 审核后，才可进行更精细的候选解释；任何阶段均不得将 HMM 命中直接表述为降解表型证据。
