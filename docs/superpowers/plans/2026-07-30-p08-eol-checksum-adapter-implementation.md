# P08 EOL Checksum Adapter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 允许 P08 在 Git 将受版本控制 FASTA 的 CRLF/LF 工作树换行转换后，仍严格验证 P05 seed/control 的来源；任何非换行字节差异继续失败关闭。

**Architecture:** P08 先验证 registry 声明 SHA 与实际 raw 文件 SHA。仅在不相等时，才由实际字节串确定性生成 LF 与 CRLF 两个候选表示；registry SHA 必须精确匹配其中一个候选。候选/参考 manifest 同时保存声明、观察与 canonical-LF 哈希及验证模式，绝不改写 P05 registry 或输入 FASTA。

**Tech Stack:** Python 3 standard library、`unittest`、TSV manifests、Git。

## Global Constraints

- 只适用于 P05 seed/control 的 file SHA-256；P06/P07 r8 manifests、P05 registry、FASTA 文件、T141 旧主 worktree 均不得改写。
- 放行条件只能是 `actual_raw == declared`、`sha256(LF-normalized actual) == declared` 或 `sha256(CRLF-normalized actual) == declared`；不可接受空白、header、sequence 或任意非 EOL 字节变化。
- 保留 declared SHA、observed raw SHA、canonical-LF SHA、验证模式与 declaration 匹配表示；这些是溯源证据，不是 PHB/PHA 降解表型证据。
- 不运行 MAFFT、FastTree、IQ-TREE、rooting 或 T141 `--preflight-only`。

---

### Task 1: 以失败测试固定 EOL-only 放行边界

**Files:**
- Modify: `tests/test_p08_prepare_phylogeny.py: checksum fixture helpers and PrepareP08InputsTests`

**Interfaces:**
- Consumes: an LF seed/control FASTA plus a registry `sequence_sha256` calculated from its CRLF representation.
- Produces: regression tests for EOL-only acceptance, non-EOL failure closure, and exact-file provenance mode.

- [ ] **Step 1: Write the EOL-only acceptance test.**

  ```python
  def test_crlf_declared_seed_hash_accepts_lf_worktree_file_and_records_audit_fields(self) -> None:
      seed = self.seed_archaeal
      raw = seed.read_bytes()
      crlf_sha = hashlib.sha256(raw.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n")).hexdigest()
      rows = read_tsv(self.seeds)
      rows[0]["sequence_sha256"] = crlf_sha
      write_tsv(self.seeds, tuple(rows[0].keys()), rows)
      references = read_tsv(self._prepare()["family_reference_manifest"])
      reference = next(row for row in references if row["record_id"] == rows[0]["seed_id"])
      self.assertEqual(reference["checksum_verification_mode"], "eol_normalized_file_sha256")
      self.assertEqual(reference["checksum_declared_representation"], "crlf")
      self.assertEqual(reference["observed_file_sha256"], hashlib.sha256(raw).hexdigest())
  ```

- [ ] **Step 2: Write the non-EOL mutation rejection test.**

  ```python
  def test_non_eol_seed_mutation_does_not_match_declared_crlf_hash(self) -> None:
      seed = self.seed_archaeal
      original = seed.read_bytes()
      declared = hashlib.sha256(original.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n")).hexdigest()
      seed.write_bytes(original.replace(b"MKKK", b"MKKR"))
      rows = read_tsv(self.seeds)
      rows[0]["sequence_sha256"] = declared
      write_tsv(self.seeds, tuple(rows[0].keys()), rows)
      with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
          self._prepare()
  ```

- [ ] **Step 3: Run the two new tests before implementation.**

  Run:

  ```powershell
  python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_crlf_declared_seed_hash_accepts_lf_worktree_file_and_records_audit_fields tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_non_eol_seed_mutation_does_not_match_declared_crlf_hash -v
  ```

  Expected: EOL-only test fails under existing raw-file-only verification; non-EOL test continues to fail closed.

### Task 2: 实现严格 EOL-only 验证与 manifest 审计

**Files:**
- Modify: `scripts/p08_prepare_phylogeny.py: checksum helper, P05 reference loading, reference and family-input manifest fields`
- Test: `tests/test_p08_prepare_phylogeny.py: Task 1 tests`

**Interfaces:**
- Produces: `_verify_file_sha256_with_eol_adapter(path: Path, declared_sha256: str) -> dict[str, str]`.
- Returns: `declared_sequence_sha256`, `observed_file_sha256`, `canonical_lf_sha256`, `checksum_verification_mode`, and `checksum_declared_representation`.

- [ ] **Step 1: Implement the verification helper.**

  ```python
  def _verify_file_sha256_with_eol_adapter(path: Path, declared_sha256: str) -> dict[str, str]:
      raw = path.read_bytes()
      observed = hashlib.sha256(raw).hexdigest()
      lf = raw.replace(b"\r\n", b"\n")
      canonical_lf = hashlib.sha256(lf).hexdigest()
      if observed == declared_sha256:
          return {"observed_file_sha256": observed, "canonical_lf_sha256": canonical_lf, "checksum_verification_mode": "exact_file_sha256", "checksum_declared_representation": "raw"}
      crlf = lf.replace(b"\n", b"\r\n")
      matched = "lf" if canonical_lf == declared_sha256 else "crlf" if hashlib.sha256(crlf).hexdigest() == declared_sha256 else ""
      if not matched:
          raise ValueError("SHA-256 mismatch")
      return {"observed_file_sha256": observed, "canonical_lf_sha256": canonical_lf, "checksum_verification_mode": "eol_normalized_file_sha256", "checksum_declared_representation": matched}
  ```

- [ ] **Step 2: Replace only P05 seed/control raw-file verification with the helper.**

  Preserve `residue_sha256` verification exactly as implemented. On mismatch call `_fail` with the existing reason and include declared/observed/canonical hashes in review notes.

- [ ] **Step 3: Add the five audit fields to reference and family-input manifests.**

  Keep `sequence_sha256` as the declared registry SHA and `verified_sha256` as the observed raw SHA for backwards readability. Add explicit fields for canonical LF SHA, verification mode, and declared representation.

- [ ] **Step 4: Run the focused tests and P08 preparation suite.**

  Run:

  ```powershell
  python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_crlf_declared_seed_hash_accepts_lf_worktree_file_and_records_audit_fields tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_non_eol_seed_mutation_does_not_match_declared_crlf_hash tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_seed_or_control_checksum_mismatch_blocks_and_writes_review -v
  python -m unittest tests/test_p08_prepare_phylogeny.py -v
  ```

  Expected: EOL-only fixture passes with explicit audit mode; both non-EOL and arbitrary checksum mismatch fixtures fail closed.

### Task 3: 记录跨平台 provenance 语义并验证集成

**Files:**
- Modify: `07_phylogeny/README.md: EOL-only checksum scope`
- Modify: `docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md: checksum admissibility and audit fields`
- Create: `docs/superpowers/plans/2026-07-30-p08-eol-checksum-adapter-implementation.md`

- [ ] **Step 1: Document that the adapter accepts only Git CRLF/LF worktree conversion and records both identities; it does not normalize biological content or weaken non-EOL integrity checks.**

- [ ] **Step 2: Run final validation.**

  ```powershell
  python -m unittest discover -s tests -v
  python scripts/validate_repository.py
  git diff --check
  ```

  Expected: all tests pass, repository validation passes, and `git diff --check` is silent.

- [ ] **Step 3: Commit only the five implementation/docs/test/plan files with a Chinese-first message.**

  ```powershell
  git add scripts/p08_prepare_phylogeny.py tests/test_p08_prepare_phylogeny.py 07_phylogeny/README.md docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md docs/superpowers/plans/2026-07-30-p08-eol-checksum-adapter-implementation.md
  git commit -m "P08：兼容 Git 换行转换的 seed checksum 溯源"
  ```
