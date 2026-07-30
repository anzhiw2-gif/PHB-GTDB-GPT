# P08 Domain-to-Target Collapse Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 让 P08 将同一 approved family、proteome shard 和 target 的 P06 HMMER 多 domain 记录确定性地汇总为一条可追溯的 target-level 候选，且不改变 P06/P07 原始账本。

**Architecture:** P06 仍是 domain-row 账本；`prepare_p08_inputs` 只对本次选择的 tier 建立 target key 分组。分组先执行 target-level 一致性门控，随后以 tier、coverage、完整行词典序决定代表 domain，并将汇总审计字段写入 P08 candidate manifest。

**Tech Stack:** Python 3 standard library、`unittest`、TSV manifests、Git。

## Global Constraints

- 只允许 `approved_for_p06=yes` 且 `scan_permission=approved` 的 P05 模型；序列、domain、taxonomy 与树位置均是序列/注释证据，不是 PHB/PHA 降解表型证明。
- 不修改或重写 T141 r8 P06/P07 manifests，也不接触旧主 worktree `/home/data/haoyu/PHB-GTDB-GPT`。
- 汇总 key 固定为 `(family_category, proteome_shard, target_id)`；不同 family 的命中不得合并。
- 必须逐字段一致地核验 `target_accession`、`target_length`、`full_sequence_score`、`calibrated_full_score_threshold` 和 `calibrated_hmm_coverage_threshold`；不一致即写 review block 并失败关闭。
- 当显式选择多个 tier 时，`High-confidence` 优先于 `Review`；同 tier 选择最大 `hmm_coverage`，完全并列以完整行词典序决胜。
- 不运行 MAFFT、FastTree、IQ-TREE、rooting 或 T141 `--preflight-only`；P08 生成物保持 machine-local/Git ignored。

---

### Task 1: 以失败测试固定 domain-row 汇总契约

**Files:**
- Modify: `tests/test_p08_prepare_phylogeny.py: P06 fixture and PrepareP08InputsTests`
- Modify: `scripts/p08_prepare_phylogeny.py: remove any uncommitted target-collapse implementation before RED`

**Interfaces:**
- Consumes: `prepare_p08_inputs(..., include_tiers=...)` and the existing temporary TSV fixture helpers.
- Produces: three regression tests whose expected behavior is independent of the later helper implementation.

- [ ] **Step 1: Remove the uncommitted production helper and its uncommitted regression tests with `apply_patch`.**

  Restore the source to the `f6e36fd` behavior that rejects duplicate selected P06 keys, and retain no test asserting collapse. Do not use `git checkout --`, `git reset`, or modify any r8 input.

- [ ] **Step 2: Write the failing tests in `PrepareP08InputsTests`.**

  Add the following three tests using real temporary TSV inputs:

  ```python
  def test_selected_multidomain_p06_rows_collapse_to_one_target(self) -> None:
      rows = read_tsv(self.p06)
      extra = dict(next(row for row in rows if row["target_id"] == "bac1"))
      extra.update({"domain_index": "2", "hmm_coverage": "0.6"})
      rows.append(extra)
      write_tsv(self.p06, P06_FIELDS, rows)
      candidates = read_tsv(self._prepare()["candidate_manifest"])
      bac1 = next(row for row in candidates if row["target_id"] == "bac1")
      self.assertEqual([row["target_id"] for row in candidates].count("bac1"), 1)
      self.assertEqual(bac1["p06_target_domain_row_count"], "2")
  ```

  ```python
  def test_selected_multidomain_p06_rows_with_conflicting_target_level_values_fail_closed(self) -> None:
      rows = read_tsv(self.p06)
      conflict = dict(next(row for row in rows if row["target_id"] == "bac1"))
      conflict["target_length"] = "9"
      rows.append(conflict)
      write_tsv(self.p06, P06_FIELDS, rows)
      with self.assertRaisesRegex(ValueError, "P06 target-level field mismatch"):
          self._prepare()
  ```

  ```python
  def test_high_confidence_domain_outweighs_review_domain_for_one_target(self) -> None:
      rows = read_tsv(self.p06)
      review = dict(next(row for row in rows if row["target_id"] == "bac1"))
      review.update({"tier": "Review", "domain_index": "2", "hmm_coverage": "0.95"})
      rows.append(review)
      write_tsv(self.p06, P06_FIELDS, rows)
      bac1 = next(row for row in read_tsv(self._prepare(include_tiers=("High-confidence", "Review"))["candidate_manifest"]) if row["target_id"] == "bac1")
      self.assertEqual(bac1["tier"], "High-confidence")
      self.assertEqual(bac1["p06_observed_selected_tiers"], "High-confidence;Review")
  ```

- [ ] **Step 3: Run each test before production implementation.**

  Run:

  ```powershell
  python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_selected_multidomain_p06_rows_collapse_to_one_target -v
  python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_selected_multidomain_p06_rows_with_conflicting_target_level_values_fail_closed -v
  python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_high_confidence_domain_outweighs_review_domain_for_one_target -v
  ```

  Expected: the first and third fail because baseline P08 rejects duplicate target keys; the second fails because baseline error text does not identify the target-level inconsistency.

### Task 2: 实现最小、确定性的 P06 target 汇总

**Files:**
- Modify: `scripts/p08_prepare_phylogeny.py: constants, selected-P06 normalization, candidate manifest projection`
- Test: `tests/test_p08_prepare_phylogeny.py: three Task 1 tests`

**Interfaces:**
- Consumes: selected P06 rows conforming to `P06_FIELDS` plus optional `domain_index`.
- Produces: `_collapse_selected_p06_domain_rows(rows, *, outdir, blocks, source_path) -> list[dict[str, str]]`; each returned row has `_p06_target_domain_row_count` and `_p06_observed_selected_tiers`.

- [ ] **Step 1: Add explicit contract constants.**

  ```python
  P06_TARGET_LEVEL_FIELDS = (
      "target_accession", "target_length", "full_sequence_score",
      "calibrated_full_score_threshold", "calibrated_hmm_coverage_threshold",
  )
  P06_DOMAIN_SELECTION_RULE = "highest_tier_then_max_hmm_coverage_then_full_row_lexical"
  ```

- [ ] **Step 2: Implement the helper and call it immediately after selected-tier filtering.**

  ```python
  def _collapse_selected_p06_domain_rows(rows, *, outdir, blocks, source_path):
      groups = defaultdict(list)
      for row in rows:
          groups[(row["family_category"], row["proteome_shard"], row["target_id"])].append(row)
      collapsed = []
      for key, group in sorted(groups.items()):
          for field in P06_TARGET_LEVEL_FIELDS:
              if len({row[field] for row in group}) != 1:
                  _fail(outdir, blocks, "P06 target-level field mismatch", family_category=key[0], proteome_shard=key[1], target_id=key[2], source_path=str(source_path), notes=f"field={field}")
          rank = {"High-confidence": 2, "Review": 1}
          best_tier = max(rank[row["tier"]] for row in group)
          best = [row for row in group if rank[row["tier"]] == best_tier]
          representative = sorted(best, key=lambda row: (-float(row["hmm_coverage"]), tuple(row[field] for field in sorted(row))))[0]
          observed = ";".join(tier for tier in ("High-confidence", "Review") if any(row["tier"] == tier for row in group))
          collapsed.append({**representative, "_p06_target_domain_row_count": str(len(group)), "_p06_observed_selected_tiers": observed})
      return collapsed
  ```

  Then project `p06_target_domain_row_count`, `p06_observed_selected_tiers`, `p06_selected_domain_index`, and `p06_domain_selection_rule` into every P08 candidate row. Preserve existing failure checks for unknown/blocked models and P06/P07 joins.

- [ ] **Step 3: Run the three regression tests and the full P08 preparation test module.**

  Run:

  ```powershell
  python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_selected_multidomain_p06_rows_collapse_to_one_target -v
  python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_selected_multidomain_p06_rows_with_conflicting_target_level_values_fail_closed -v
  python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_high_confidence_domain_outweighs_review_domain_for_one_target -v
  python -m unittest tests/test_p08_prepare_phylogeny.py -v
  ```

  Expected: all pass; no external phylogeny executable is invoked.

### Task 3: 固化可追溯性说明并完成仓库级验证

**Files:**
- Modify: `07_phylogeny/README.md: P06/P08 evidence and domain-row explanation`
- Modify: `docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md: admissibility and error behavior`
- Create: `docs/superpowers/plans/2026-07-30-p08-domain-target-collapse-implementation.md`

**Interfaces:**
- Consumes: the helper contract and manifest audit fields from Task 2.
- Produces: a Chinese-first explanation that distinguishes target-level sequence records from domain-row evidence.

- [ ] **Step 1: Document the exact key, fail-closed fields, precedence, deterministic tie-breaker, and four audit columns.**

  State that multiple domains neither create multiple candidate sequences nor establish a phenotype. State that P06/P07 r8 source artifacts are unchanged.

- [ ] **Step 2: Run repository-wide validation.**

  Run:

  ```powershell
  python -m unittest discover -s tests -v
  python scripts/validate_repository.py
  git diff --check
  ```

  Expected: 173 tests pass, repository validation passes, and `git diff --check` is silent.

- [ ] **Step 3: Commit only the four implementation/docs/test files and this plan.**

  ```powershell
  git add scripts/p08_prepare_phylogeny.py tests/test_p08_prepare_phylogeny.py 07_phylogeny/README.md docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md docs/superpowers/plans/2026-07-30-p08-domain-target-collapse-implementation.md
  git commit -m "P08：按 target 汇总 P06 多 domain 候选"
  ```

  Do not stage `docs/superpowers/specs/2026-07-29-p08-phylogeny-scaffold-design.md`, `docs/superpowers/plans/` files already present on `main`, or the untracked file `1`.

## Self-Review

- Spec coverage: Task 1 proves the three observed P06 failure modes; Task 2 enforces every target-level safety gate and deterministic representative rule; Task 3 records the evidence boundary and verifies repository integration.
- Placeholder scan: no `TBD`, `TODO`, or unspecified validation remains.
- Type consistency: all paths use `dict[str, str]`, the helper returns `list[dict[str, str]]`, and its two audit keys are consumed by the candidate-manifest projection.
