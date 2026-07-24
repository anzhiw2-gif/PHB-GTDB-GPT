"""P01 GTDB audit helpers.

This module validates the server-side path template and builds compact file
manifests for the raw GTDB copy stage.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterable, Sequence


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
    return _parse_paths_file(config_path.read_text(encoding="utf-8"), config_path)


def _parse_paths_file(text: str, config_path: Path) -> dict[str, str]:
    paths: dict[str, str] = {}
    current_section: str | None = None

    for raw_line in text.splitlines():
        line = raw_line.rstrip()
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if not line.startswith(" "):
            if not stripped.endswith(":"):
                raise ValueError(f"{config_path} has invalid top-level line: {raw_line}")
            current_section = stripped[:-1]
            continue
        if current_section != "paths":
            continue
        if ":" not in stripped:
            raise ValueError(f"{config_path} has invalid path line: {raw_line}")
        key, value = stripped.split(":", 1)
        paths[key.strip()] = value.strip()

    if current_section is None:
        raise ValueError(f"{config_path} must contain a top-level section")
    if not paths:
        raise ValueError(f"{config_path} must contain a 'paths' mapping")
    return paths


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


def iter_file_manifest(root: Path, max_workers: int | None = None) -> Iterable[dict[str, object]]:
    root = root.resolve(strict=False)
    files = sorted(p for p in root.rglob("*") if p.is_file())
    worker_count = min(max_workers or 60, os.cpu_count() or 1, max(1, len(files)))
    if worker_count == 1:
        for path in files:
            yield _manifest_row(root, path)
        return

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        for row in executor.map(lambda path: _manifest_row(root, path), files):
            yield row


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


def copy_raw_genomes(source: Path, target: Path, runner: CommandRunner = subprocess.run) -> None:
    target.mkdir(parents=True, exist_ok=True)
    command = [
        "rsync",
        "-aHAX",
        "--info=progress2",
        f"{source.resolve(strict=False).as_posix().rstrip('/')}/",
        f"{target.resolve(strict=False).as_posix().rstrip('/')}/",
    ]
    result = runner(command, capture_output=True, text=True, check=False)
    if getattr(result, "returncode", 1) != 0:
        raise RuntimeError(
            "rsync failed: "
            f"{getattr(result, 'stdout', '')}{getattr(result, 'stderr', '')}".strip()
        )


def write_manifest(manifest_path: Path, records: Iterable[dict[str, object]]) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["relative_path", "bytes", "mtime_utc", "sha256"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def verify_sampled_hashes(source_root: Path, target_root: Path, sample_paths: Sequence[Path]) -> int:
    verified = 0
    for source_path in sample_paths:
        relative = source_path.resolve(strict=False).relative_to(source_root.resolve(strict=False))
        target_path = target_root.resolve(strict=False) / relative
        if not target_path.is_file():
            raise RuntimeError(f"missing copied file: {target_path}")
        if _sha256(source_path) != _sha256(target_path):
            raise RuntimeError(f"checksum mismatch: {relative.as_posix()}")
        verified += 1
    return verified


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


def _manifest_row(root: Path, path: Path) -> dict[str, object]:
    payload = path.read_bytes()
    stat = path.stat()
    return {
        "relative_path": path.relative_to(root).as_posix(),
        "bytes": len(payload),
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat().replace("+00:00", "Z"),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def run_p01_audit(paths: dict[str, str], source: Path, target: Path, threads: int) -> dict[str, object]:
    source_summary = summarize_tree(source)
    free_bytes = available_disk_bytes(target.parent)
    ensure_sufficient_free_space(free_bytes, int(source_summary["byte_count"]))

    copy_raw_genomes(source, target)
    support_files = copy_support_files(paths, target.parent)

    target_summary = summarize_tree(target)
    if source_summary["file_count"] != target_summary["file_count"]:
        raise RuntimeError(
            f"file count mismatch: source={source_summary['file_count']} target={target_summary['file_count']}"
        )
    if source_summary["byte_count"] != target_summary["byte_count"]:
        raise RuntimeError(
            f"byte count mismatch: source={source_summary['byte_count']} target={target_summary['byte_count']}"
        )
    if source_summary["top_level_bytes"] != target_summary["top_level_bytes"]:
        raise RuntimeError("top-level byte counts do not match between source and target")

    sample_paths = select_checksum_sample(source, minimum=1000, fraction=0.01)
    verified = verify_sampled_hashes(source, target, sample_paths)
    manifest_path = target.parent / "manifests" / "raw_genomes_manifest.tsv"
    write_manifest(manifest_path, iter_file_manifest(target, max_workers=threads))
    tool_versions = collect_tool_versions(build_tool_command_map())

    return {
        "source_summary": source_summary,
        "target_summary": target_summary,
        "sample_verified": verified,
        "manifest_path": manifest_path,
        "support_files": support_files,
        "tool_versions": tool_versions,
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
    parser.add_argument("--threads", type=int, default=60, help="Maximum manifest hashing threads, capped at 60")
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

    result = run_p01_audit(paths, source, target, threads=min(60, args.threads))

    print(f"Validated P01 copy plan from {source} to {target}.")
    print(f"Source file count: {result['source_summary']['file_count']}")
    print(f"Source byte count: {result['source_summary']['byte_count']}")
    print(f"Target file count: {result['target_summary']['file_count']}")
    print(f"Target byte count: {result['target_summary']['byte_count']}")
    print(f"Sampled checksum verifications: {result['sample_verified']}")
    print(f"Manifest written: {result['manifest_path']}")
    print(f"Collected tool versions for {len(result['tool_versions'])} tools.")
    copied_support = result["support_files"]
    if copied_support:
        print(
            "Copied support files: "
            + ", ".join(sorted(key for key, value in copied_support.items() if value is not None))
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
