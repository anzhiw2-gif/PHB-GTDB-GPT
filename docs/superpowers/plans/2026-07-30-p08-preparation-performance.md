# P08 Preparation Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate repeated large-file SHA-256 work and allow deterministic P08 FASTA preloading with at most 60 workers.

**Architecture:** An invocation-local digest cache holds a canonical path, digest, and immutable file identity. A bounded `ThreadPoolExecutor` preloads separate P07 FASTA shards, while ordered validation and every output write remain on the main thread.

**Tech Stack:** Python standard library (`dataclasses`, `concurrent.futures`, `pathlib`, `hashlib`, `unittest`).

## Global Constraints

- `workers` must be an integer from 1 through 60, and defaults to 1.
- A cached input whose device, inode, size, or modification time changes must fail closed.
- P05/P06/P07/P03/GTDB provenance checks and sequence-evidence-only language remain unchanged.
- Workers must never write shared P08 manifests or FASTA files.
- Preparation must not execute MAFFT, FastTree, IQ-TREE, rooting, or biological analysis.

---

### Task 1: Add a stable input digest cache

**Files:**

- Modify: `scripts/p08_prepare_phylogeny.py`
- Modify: `tests/test_p08_prepare_phylogeny.py`

**Interfaces:**

- Produces `_InputDigestCache.sha256(path: Path) -> str`.
- Raises `ValueError("input changed during P08 preparation: ...")` on a changed cached input.

- [ ] **Step 1: Write the failing test**

```python
def test_input_digest_cache_reuses_stable_hash_and_blocks_mutation(self) -> None:
    cache = preparer._InputDigestCache()
    path = self.root / "stable.txt"
    path.write_text("first\n", encoding="utf-8")
    self.assertEqual(cache.sha256(path), self._sha256(path))
    self.assertEqual(cache.sha256(path), self._sha256(path))
    path.write_text("second\n", encoding="utf-8")
    with self.assertRaisesRegex(ValueError, "input changed during P08 preparation"):
        cache.sha256(path)
```

- [ ] **Step 2: Run it and verify RED**

Run: `python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_input_digest_cache_reuses_stable_hash_and_blocks_mutation -v`

Expected: FAIL because `_InputDigestCache` is absent.

- [ ] **Step 3: Implement the smallest cache**

```python
@dataclass(frozen=True)
class _InputFingerprint:
    device: int
    inode: int
    size: int
    mtime_ns: int

class _InputDigestCache:
    def sha256(self, path: Path) -> str:
        canonical = Path(path).resolve(strict=True)
        fingerprint = _input_fingerprint(canonical)
        if canonical in self._entries and self._entries[canonical][0] != fingerprint:
            raise ValueError(f"input changed during P08 preparation: {canonical}")
        return self._entries.setdefault(canonical, (fingerprint, _sha256(canonical)))[1]
```

- [ ] **Step 4: Run it and verify GREEN**

Run: `python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_input_digest_cache_reuses_stable_hash_and_blocks_mutation -v`

Expected: PASS.

### Task 2: Add bounded deterministic FASTA preloading

**Files:**

- Modify: `scripts/p08_prepare_phylogeny.py`
- Modify: `tests/test_p08_prepare_phylogeny.py`

**Interfaces:**

- Extends `prepare_p08_inputs(..., workers: int = 1)`.
- Extends CLI with `--workers`.
- Produces equivalent candidate-manifest rows for `workers=1` and `workers=60`.

- [ ] **Step 1: Write failing tests**

```python
def test_prepare_accepts_sixty_workers_with_deterministic_candidate_rows(self) -> None:
    one = read_tsv(self._prepare(outdir=self.root / "one", workers=1)["candidate_manifest"])
    sixty = read_tsv(self._prepare(outdir=self.root / "sixty", workers=60)["candidate_manifest"])
    self.assertEqual(one, sixty)

def test_prepare_rejects_worker_counts_outside_one_through_sixty(self) -> None:
    for workers in (0, 61):
        with self.assertRaisesRegex(ValueError, "workers must be between 1 and 60"):
            self._prepare(workers=workers)
```

- [ ] **Step 2: Run them and verify RED**

Run: `python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_prepare_accepts_sixty_workers_with_deterministic_candidate_rows tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_prepare_rejects_worker_counts_outside_one_through_sixty -v`

Expected: FAIL because the API does not accept `workers`.

- [ ] **Step 3: Implement bounded preloading**

```python
def _preload_candidate_fastas(paths: Iterable[Path], workers: int, cache: _InputDigestCache) -> dict[str, tuple[str, dict[str, str]]]:
    ordered = sorted({Path(path).resolve(strict=True) for path in paths}, key=str)
    with ThreadPoolExecutor(max_workers=min(workers, len(ordered) or 1)) as executor:
        values = executor.map(lambda path: (cache.sha256(path), _read_fasta(path)), ordered)
        return {str(path): value for path, value in zip(ordered, values)}
```

Validate the worker range before loading inputs. Replace repeated candidate-row `_sha256(...)` calls with the cache. Use a lock inside the cache because preloading calls it from several threads.

- [ ] **Step 4: Run them and verify GREEN**

Run: `python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_prepare_accepts_sixty_workers_with_deterministic_candidate_rows tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_prepare_rejects_worker_counts_outside_one_through_sixty -v`

Expected: PASS.

### Task 3: Record worker provenance and perform regression validation

**Files:**

- Modify: `scripts/p08_prepare_phylogeny.py`
- Modify: `tests/test_p08_prepare_phylogeny.py`
- Modify: `07_phylogeny/README.md`
- Modify: `docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md`
- Modify: `AGENTS.md`

**Interfaces:**

- Produces `p08_requested_workers` and `p08_effective_fasta_preload_workers` provenance rows.

- [ ] **Step 1: Write the failing provenance test**

```python
def test_worker_provenance_records_requested_and_effective_counts(self) -> None:
    provenance = {row["input_role"]: row for row in read_tsv(self._prepare(workers=60)["input_provenance"])}
    self.assertEqual(provenance["p08_requested_workers"]["input_path"], "60")
    self.assertEqual(provenance["p08_effective_fasta_preload_workers"]["input_path"], "2")
```

- [ ] **Step 2: Run it and verify RED**

Run: `python -m unittest tests.test_p08_prepare_phylogeny.PrepareP08InputsTests.test_worker_provenance_records_requested_and_effective_counts -v`

Expected: FAIL because the provenance rows are absent.

- [ ] **Step 3: Implement and document the contract**

Append compact worker provenance rows, document the single-writer rule and the fresh T141 form:

```bash
python scripts/p08_prepare_phylogeny.py ... --workers 60 --outdir <new-unused-directory>
```

- [ ] **Step 4: Run all checks**

Run:

```powershell
python -m unittest tests/test_p08_prepare_phylogeny.py -v
python -m unittest tests/test_p08_run_phylogeny.py -v
python -m unittest discover -s tests -v
python scripts/validate_repository.py
git diff --check
```

Expected: all tests pass, the repository is valid, and the whitespace check is silent.

- [ ] **Step 5: Commit**

```powershell
git add AGENTS.md 07_phylogeny/README.md docs/P08_PHYLOGENY_TAXONOMY_PLAN_2026-07-29.md scripts/p08_prepare_phylogeny.py tests/test_p08_prepare_phylogeny.py
git commit -m "P08：加入确定性并行预检与哈希缓存"
```

