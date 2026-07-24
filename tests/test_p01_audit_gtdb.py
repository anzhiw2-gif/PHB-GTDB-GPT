from __future__ import annotations

import hashlib
import tempfile
import textwrap
from pathlib import Path
import unittest

from scripts.p01_audit_gtdb import (
    collect_tool_versions,
    copy_support_files,
    ensure_sufficient_free_space,
    find_ar53_tree,
    iter_file_manifest,
    load_paths_config,
    select_checksum_sample,
    summarize_tree,
    validate_copy_plan,
)


class P01AuditGTDBTest(unittest.TestCase):
    def test_load_paths_config_reads_nested_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            config_path = Path(tmp) / "paths.yaml"
            config_path.write_text(
                textwrap.dedent(
                    """
                    paths:
                      project_dir: /home/data/haoyu/PHB-GTDB-GPT
                      gtdb_root: /home/data/haoyu/GTDB
                      gtdb_genomes_source: /home/data/haoyu/GTDB/gtdb_genomes_reps_r232/database
                      old_project_readonly: /home/data/haoyu/PHB_gtdb
                      bac120_taxonomy_source: /home/data/haoyu/GTDB/taxonomy/bac120_taxonomy_r232.tsv
                      ar53_taxonomy_source: /home/data/haoyu/GTDB/taxonomy/ar53_taxonomy_r232.tsv
                      bac120_tree_source: /home/data/haoyu/GTDB/GTDB_tree/bac120_r232.tree
                    """
                ).strip()
                + "\n",
                encoding="utf-8",
            )

            paths = load_paths_config(config_path)

            self.assertEqual(paths["project_dir"], "/home/data/haoyu/PHB-GTDB-GPT")
            self.assertEqual(paths["gtdb_root"], "/home/data/haoyu/GTDB")

    def test_validate_copy_plan_rejects_missing_source_target_match_and_escape(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            project_dir = root / "project"
            project_dir.mkdir()
            source_dir = root / "source"
            source_dir.mkdir()
            target_inside = project_dir / "00_raw_gtdb_r232" / "genomes"

            missing_errors = validate_copy_plan(root / "missing", target_inside, project_dir)
            self.assertIn("missing source", " | ".join(missing_errors))

            equal_errors = validate_copy_plan(project_dir, project_dir, project_dir)
            self.assertIn("source and target must differ", " | ".join(equal_errors))

            escape_errors = validate_copy_plan(source_dir, root / "outside", project_dir)
            self.assertIn("outside project_dir", " | ".join(escape_errors))

    def test_iter_file_manifest_records_relative_paths_and_sha256(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "payload"
            root.mkdir()
            (root / "alpha.txt").write_bytes(b"alpha\n")
            nested = root / "nested"
            nested.mkdir()
            (nested / "beta.bin").write_bytes(b"\x00\x01")

            records = list(iter_file_manifest(root))

            self.assertEqual([record["relative_path"] for record in records], ["alpha.txt", "nested/beta.bin"])
            self.assertEqual(records[0]["bytes"], 6)
            self.assertEqual(records[0]["sha256"], hashlib.sha256(b"alpha\n").hexdigest())
            self.assertEqual(records[1]["sha256"], hashlib.sha256(b"\x00\x01").hexdigest())

    def test_summarize_tree_reports_top_level_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "payload"
            root.mkdir()
            (root / "root.txt").write_bytes(b"abc")
            nested = root / "nested"
            nested.mkdir()
            (nested / "leaf.txt").write_bytes(b"defg")

            summary = summarize_tree(root)

            self.assertEqual(summary["file_count"], 2)
            self.assertEqual(summary["byte_count"], 7)
            self.assertEqual(summary["top_level_bytes"], {"nested": 4, "root.txt": 3})

    def test_select_checksum_sample_respects_minimum_and_fraction(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "payload"
            root.mkdir()
            for index in range(5):
                (root / f"f{index}.txt").write_bytes(f"{index}".encode("ascii"))

            sample = select_checksum_sample(root, minimum=1000, fraction=0.01)

            self.assertEqual([path.name for path in sample], [f"f{index}.txt" for index in range(5)])

    def test_ensure_sufficient_free_space_raises_when_insufficient(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "insufficient free space"):
            ensure_sufficient_free_space(free_bytes=10, required_bytes=11)

    def test_collect_tool_versions_uses_runner(self) -> None:
        calls: list[list[str]] = []

        def fake_runner(command: list[str], **_: object) -> object:
            calls.append(command)

            class Result:
                returncode = 0
                stdout = "tool 1.2.3\n"
                stderr = ""

            return Result()

        versions = collect_tool_versions({"python": ["python", "--version"]}, runner=fake_runner)

        self.assertEqual(calls, [["python", "--version"]])
        self.assertEqual(versions["python"], "tool 1.2.3")

    def test_find_ar53_tree_picks_first_sorted_match(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            later = root / "z_ar53.tree"
            earlier = root / "a_ar53.tree"
            later.write_text("later", encoding="utf-8")
            earlier.write_text("earlier", encoding="utf-8")

            self.assertEqual(find_ar53_tree(root), earlier)

    def test_copy_support_files_copies_expected_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            project = root / "project"
            project.mkdir()
            source_root = root / "gtdb"
            source_root.mkdir()
            raw_target = project / "00_raw_gtdb_r232"
            raw_target.mkdir(parents=True)

            bac_tax = source_root / "bac120_taxonomy_r232.tsv"
            ar_tax = source_root / "ar53_taxonomy_r232.tsv"
            bac_tree = source_root / "bac120_r232.tree"
            ar_tree = source_root / "nested" / "ar53_r232.tree"
            ar_tree.parent.mkdir()
            bac_tax.write_text("bac", encoding="utf-8")
            ar_tax.write_text("ar", encoding="utf-8")
            bac_tree.write_text("tree", encoding="utf-8")
            ar_tree.write_text("ar-tree", encoding="utf-8")

            paths = {
                "gtdb_root": str(source_root),
                "bac120_taxonomy_source": str(bac_tax),
                "ar53_taxonomy_source": str(ar_tax),
                "bac120_tree_source": str(bac_tree),
            }

            copied = copy_support_files(paths, raw_target)

            self.assertEqual((raw_target / "bac120_taxonomy_r232.tsv").read_text(encoding="utf-8"), "bac")
            self.assertEqual((raw_target / "ar53_taxonomy_r232.tsv").read_text(encoding="utf-8"), "ar")
            self.assertEqual((raw_target / "bac120_r232.tree").read_text(encoding="utf-8"), "tree")
            self.assertEqual((raw_target / "ar53_r232.tree").read_text(encoding="utf-8"), "ar-tree")
            self.assertEqual(copied["ar53_tree"], raw_target / "ar53_r232.tree")


if __name__ == "__main__":
    unittest.main()
