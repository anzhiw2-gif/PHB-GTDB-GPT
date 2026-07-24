"""P01 GTDB audit helpers.

This module validates the server-side path template and builds compact file
manifests for the raw GTDB copy stage.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Sequence

import yaml


REQUIRED_PATH_KEYS = (
    "project_dir",
    "gtdb_root",
    "gtdb_genomes_source",
    "old_project_readonly",
    "bac120_taxonomy_source",
    "ar53_taxonomy_source",
    "bac120_tree_source",
)

CommandRunner = Callable[..., object]


def load_paths_config(config_path: Path) -> dict[str, str]:
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{config_path} must contain a mapping at the top level")

    paths = data.get("paths")
    if not isinstance(paths, dict):
        raise ValueError(f"{config_path} must contain a 'paths' mapping")

    return {str(key): str(value) for key, value in paths.items()}


def validate_copy_plan(source: Path, target: Path, project_dir: Path) -> list[str]:
    errors: list[str] = []
    resolved_project_dir = project_dir.resolve(strict=False)
    resolved_source = source.resolve(strict=False)
    resolved_target = target.resolve(strict=False)

    if not source.exists():
        errors.append(f"missing source: {source}")
    if resolved_source == resolved_target:
        errors.append("source and target must differ")
    if not _is_within_project_dir(resolved_target, resolved_project_dir):
        errors.append("target is outside project_dir")

    return errors


def iter_file_manifest(root: Path) -> Iterable[dict[str, object]]:
    root = root.resolve(strict=False)
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        payload = path.read_bytes()
        stat = path.stat()
        yield {
            "relative_path": path.relative_to(root).as_posix(),
            "bytes": len(payload),
            "mtime_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat().replace("+00:00", "Z"),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }


def summarize_tree(root: Path) -> dict[str, object]:
    root = root.resolve(strict=False)
    file_count = 0
    byte_count = 0
    top_level_bytes: dict[str, int] = {}

    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        size = path.stat().st_size
        file_count += 1
        byte_count += size
        relative = path.relative_to(root)
        top_level = relative.parts[0]
        top_level_bytes[top_level] = top_level_bytes.get(top_level, 0) + size

    return {
        "file_count": file_count,
        "byte_count": byte_count,
        "top_level_bytes": dict(sorted(top_level_bytes.items())),
    }


def find_ar53_tree(gtdb_root: Path) -> Path | None:
    candidates = sorted(path for path in gtdb_root.resolve(strict=False).rglob("*ar53*.tree*") if path.is_file())
    return candidates[0] if candidates else None


def copy_support_files(paths: dict[str, str], raw_target: Path) -> dict[str, Path | None]:
    raw_target.mkdir(parents=True, exist_ok=True)
    copied: dict[str, Path | None] = {}

    for key in ("bac120_taxonomy_source", "ar53_taxonomy_source", "bac120_tree_source"):
        source = Path(paths[key])
        destination = raw_target / source.name
        shutil.copy2(source, destination)
        copied[key] = destination

    ar53_tree = find_ar53_tree(Path(paths["gtdb_root"]))
    copied["ar53_tree"] = None
    if ar53_tree is not None:
        destination = raw_target / ar53_tree.name
        shutil.copy2(ar53_tree, destination)
        copied["ar53_tree"] = destination

    return copied


def select_checksum_sample(root: Path, minimum: int = 1000, fraction: float = 0.01) -> list[Path]:
    files = sorted(p for p in root.resolve(strict=False).rglob("*") if p.is_file())
    if not files:
        return []

    sample_size = min(len(files), max(minimum, math.ceil(len(files) * fraction)))
    if sample_size >= len(files):
        return files

    step = len(files) / sample_size
    chosen: list[Path] = []
    index = 0.0
    while len(chosen) < sample_size:
        candidate = files[min(int(index), len(files) - 1)]
        if not chosen or candidate != chosen[-1]:
            chosen.append(candidate)
        index += step
    return chosen


def collect_tool_versions(command_map: dict[str, Sequence[str]], runner: CommandRunner = subprocess.run) -> dict[str, str]:
    versions: dict[str, str] = {}
    for tool_name, command in command_map.items():
        versions[tool_name] = _probe_command_version(command, runner=runner)
    return versions


def _probe_command_version(command: Sequence[str], runner: CommandRunner = subprocess.run) -> str:
    try:
        result = runner(
            list(command),
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        return "not found"

    stdout = str(getattr(result, "stdout", "") or "").strip()
    stderr = str(getattr(result, "stderr", "") or "").strip()
    if stdout and stderr:
        return f"{stdout} {stderr}".strip()
    return stdout or stderr or "no version output"


def available_disk_bytes(path: Path) -> int:
    return shutil.disk_usage(path).free


def ensure_sufficient_free_space(free_bytes: int, required_bytes: int) -> None:
    if free_bytes < required_bytes:
        raise RuntimeError(
            f"insufficient free space: need {required_bytes} bytes, found {free_bytes} bytes"
        )


def build_tool_command_map() -> dict[str, list[str]]:
    return {
        "python": ["python", "--version"],
        "nextflow": ["nextflow", "-version"],
        "java": ["java", "-version"],
        "slurm": ["sinfo", "--version"],
        "prodigal": ["prodigal", "-v"],
        "pyrodigal": ["python", "-c", "import pyrodigal; print(pyrodigal.__version__)"],
        "hmmer": ["hmmsearch", "-h"],
        "mafft": ["mafft", "--version"],
        "iqtree2": ["iqtree2", "-h"],
        "fasttree": ["FastTree", "-version"],
        "interproscan": ["interproscan.sh", "-version"],
        "signalp": ["signalp", "--version"],
        "phobius": ["phobius", "--version"],
    }


def _is_within_project_dir(path: Path, project_dir: Path) -> bool:
    try:
        path.relative_to(project_dir)
    except ValueError:
        return False
    return True


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Audit the P01 GTDB copy plan.")
    parser.add_argument("--paths", required=True, type=Path, help="Path to config/paths.yaml")
    parser.add_argument("--source", type=Path, help="Raw GTDB source directory")
    parser.add_argument("--target", type=Path, help="Copy target directory")
    parser.add_argument(
        "--copy-support-files",
        action="store_true",
        help="Copy taxonomy and tree support files into the raw target directory",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    paths = load_paths_config(args.paths)

    missing = [key for key in REQUIRED_PATH_KEYS if key not in paths]
    if missing:
        raise SystemExit(f"Missing required path keys: {', '.join(missing)}")

    project_dir = Path(paths["project_dir"])
    source = args.source or Path(paths["gtdb_genomes_source"])
    target = args.target or project_dir / "00_raw_gtdb_r232" / "genomes"
    errors = validate_copy_plan(source, target, project_dir)
    if errors:
        raise SystemExit("; ".join(errors))

    summary = summarize_tree(source)
    tool_versions = collect_tool_versions(build_tool_command_map())
    free_bytes = available_disk_bytes(target.parent)
    ensure_sufficient_free_space(free_bytes, int(summary["byte_count"]))
    copied_support = copy_support_files(paths, target.parent) if args.copy_support_files else {}

    print(f"Validated P01 copy plan from {source} to {target}.")
    print(f"Source file count: {summary['file_count']}")
    print(f"Source byte count: {summary['byte_count']}")
    print(f"Free bytes at target parent: {free_bytes}")
    print(f"Collected tool versions for {len(tool_versions)} tools.")
    if copied_support:
        print(f"Copied support files: {', '.join(sorted(key for key, value in copied_support.items() if value is not None))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
