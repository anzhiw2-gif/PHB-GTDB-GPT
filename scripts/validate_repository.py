"""Validate the tracked PHB-GTDB-GPT repository skeleton without extra packages."""

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
REQUIRED_FILES = (
    "README.md",
    ".gitignore",
    ".gitattributes",
    "main.nf",
    "nextflow.config",
    "config/project.yaml",
    "config/paths.example.yaml",
    "docs/PROJECT_SCOPE.md",
)
REQUIRED_DIRECTORIES = (
    "00_raw_gtdb_r232",
    "01_reference_library",
    "02_prediction_benchmark",
    "03_gtdb_proteomes",
    "04_family_profiles",
    "05_hmmer_scan",
    "06_domain_annotation",
    "07_phylogeny",
    "08_reports",
    "modules",
    "subworkflows",
    "tests",
)


def main() -> int:
    missing_files = [path for path in REQUIRED_FILES if not (ROOT / path).is_file()]
    missing_dirs = [path for path in REQUIRED_DIRECTORIES if not (ROOT / path).is_dir()]
    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8") if not missing_files else ""
    protects_raw_inputs = "00_raw_gtdb_r232/*" in gitignore

    if missing_files or missing_dirs or not protects_raw_inputs:
        if missing_files:
            print("Missing files: " + ", ".join(missing_files), file=sys.stderr)
        if missing_dirs:
            print("Missing directories: " + ", ".join(missing_dirs), file=sys.stderr)
        if not protects_raw_inputs:
            print("Raw GTDB inputs are not protected by .gitignore", file=sys.stderr)
        return 1

    print("Repository layout is valid; raw GTDB inputs are ignored by Git.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
