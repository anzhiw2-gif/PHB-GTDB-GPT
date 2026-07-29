"""Preflight P08 phylogeny command plans without biological interpretation.

The Task 3 command manifest is the sole authoritative plan input.  This module
checks local task integrity and tool availability; its statuses are never
evidence that a family lacks PHB/PHA degradation activity or any subtype.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shlex
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path


STATUS_FILENAME = "p08_phylogeny_run_status.tsv"
STATUS_VALUES = (
    "preflight_ok",
    "missing_executable",
    "missing_input",
    "checksum_mismatch",
    "completed",
    "skipped_existing",
    "failed_exit_code",
    "failed_missing_output",
)
COMMAND_MANIFEST_FIELDS = (
    "family_category", "command_status", "input_fasta_path", "input_sha256",
    "candidate_input_record_count", "total_input_record_count", "route", "alignment_fasta_path",
    "representative_input_fasta_path", "fasttree_tree_path", "iqtree_prefix", "representative_plan",
    "mafft_template", "fasttree_template", "iqtree2_template", "iqtree2_annotation",
    "rooting_policy", "evidence_boundary",
)
STATUS_FIELDS = (
    "family_category", "route", "output_path", "input_fasta_path", "input_sha256",
    "command_template", "selected_command", "executable", "status", "exit_code",
    "started_at_utc", "finished_at_utc", "notes",
)
PHYLOGENY_EXECUTABLES = frozenset({"mafft", "fasttree", "fasttreemp", "iqtree", "iqtree2"})
NON_BIOLOGICAL_NOTE = (
    "Execution and integrity status is not biological negative evidence; "
    "P08 remains sequence/annotation/tree-planning evidence only."
)


def _utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_manifest(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"{path} is missing a tabular header")
            missing = [field for field in COMMAND_MANIFEST_FIELDS if field not in reader.fieldnames]
            if missing:
                raise ValueError(f"{path} is missing Task 3 command-manifest columns: {', '.join(missing)}")
            rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]
    except OSError as error:
        raise ValueError(f"could not read command manifest {path}: {error}") from error
    if not rows:
        raise ValueError(f"{path} has no command-plan rows")
    return rows


def _selected_task(job: dict[str, str]) -> dict[str, str]:
    route = job["route"]
    if route in {"mafft_linsi_then_review", "mafft_auto_then_review"}:
        input_path = job["input_fasta_path"]
        output_path = job["alignment_fasta_path"]
        template = job["mafft_template"]
    elif route == "deterministic_representative_plan_then_fasttree_exploratory":
        input_path = job["representative_input_fasta_path"]
        output_path = job["fasttree_tree_path"]
        template = job["fasttree_template"]
    else:
        raise ValueError(f"unsupported Task 3 P08 route {route!r} for family {job['family_category']!r}")
    required = {
        "family_category": job["family_category"],
        "input_fasta_path": input_path,
        "input_sha256": job["input_sha256"],
        "output_path": output_path,
    }
    missing = [field for field, value in required.items() if not value]
    if missing:
        raise ValueError(f"Task 3 command manifest has missing selected values for {job['family_category']!r}: {', '.join(missing)}")
    return {
        "family_category": job["family_category"],
        "route": route,
        "input_fasta_path": input_path,
        "input_sha256": job["input_sha256"],
        "output_path": output_path,
        "command_template": template,
        "selected_command": "",
        "executable": "",
        "command_parts": [],
    }


def _parse_selected_command(task: dict[str, str], workers: int) -> dict[str, str]:
    template = task["command_template"]
    if not template:
        raise ValueError(f"Task 3 command manifest has missing selected command template for {task['family_category']!r}")
    format_values = {
        "input_fasta": task["input_fasta_path"],
        "alignment_fasta": task["output_path"],
        "representative_input_fasta": task["input_fasta_path"],
        "representative_alignment_fasta": task["input_fasta_path"],
        "fasttree_tree": task["output_path"],
        "iqtree_prefix": "",
        "threads": str(workers),
    }
    try:
        command = template.format(**format_values)
    except (KeyError, ValueError) as error:
        raise ValueError(f"could not render selected command template for {task['family_category']!r}: {error}") from error
    try:
        command_parts = shlex.split(command, posix=True)
    except ValueError as error:
        raise ValueError(f"could not parse selected command for {task['family_category']!r}: {error}") from error
    if not command_parts:
        raise ValueError(f"could not parse selected command for {task['family_category']!r}: command is empty")
    return {**task, "selected_command": command, "executable": command_parts[0].strip("\"'"), "command_parts": command_parts}


def _executable_available(executable: str) -> bool:
    executable_path = Path(executable)
    if executable_path.parent != Path("."):
        return executable_path.is_file()
    return shutil.which(executable) is not None


def _output_is_complete(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _status_row(task: dict[str, str], status: str, *, exit_code: int | None = None, notes: str = "") -> dict[str, str]:
    timestamp = _utc_now()
    message = NON_BIOLOGICAL_NOTE if not notes else f"{notes} {NON_BIOLOGICAL_NOTE}"
    return {
        "family_category": task["family_category"],
        "route": task["route"],
        "output_path": task["output_path"],
        "input_fasta_path": task["input_fasta_path"],
        "input_sha256": task["input_sha256"],
        "command_template": task["command_template"],
        "selected_command": task["selected_command"],
        "executable": task["executable"],
        "status": status,
        "exit_code": "" if exit_code is None else str(exit_code),
        "started_at_utc": timestamp,
        "finished_at_utc": timestamp,
        "notes": message,
    }


def _write_status_atomic(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in STATUS_FIELDS})
    os.replace(temporary, path)


def _task_key(task: dict[str, str]) -> tuple[str, str, str]:
    return (task["family_category"], task["route"], task["output_path"])


def _check_task(task: dict[str, str], *, workers: int, preflight_only: bool) -> tuple[dict[str, str] | None, dict[str, str]]:
    input_path = Path(task["input_fasta_path"])
    if not input_path.is_file() or input_path.stat().st_size == 0:
        return _status_row(task, "missing_input", notes="Input FASTA is absent or empty."), task
    if _sha256(input_path) != task["input_sha256"]:
        return _status_row(task, "checksum_mismatch", notes="Input FASTA whole-file SHA-256 differs from the Task 3 manifest."), task
    task = _parse_selected_command(task, workers)
    executable_name = Path(task["executable"]).name.lower()
    if not preflight_only and executable_name in PHYLOGENY_EXECUTABLES:
        raise ValueError(f"non-preflight execution of phylogeny executable {task['executable']!r} is not authorized")
    if not _executable_available(task["executable"]):
        return _status_row(task, "missing_executable", notes="Selected command executable was not found."), task
    if _output_is_complete(Path(task["output_path"])):
        return _status_row(task, "skipped_existing", notes="Route-specific planned output already exists and is nonempty."), task
    if preflight_only:
        return _status_row(task, "preflight_ok", notes="Input, checksum, command parsing, and executable checks passed; no command was run."), task
    return None, task


def _execute_task(task: dict[str, str]) -> dict[str, str]:
    completed = subprocess.run(task["command_parts"], check=False)
    if completed.returncode != 0:
        return _status_row(task, "failed_exit_code", exit_code=completed.returncode, notes="Generic test-only command returned a nonzero exit code.")
    if not _output_is_complete(Path(task["output_path"])):
        return _status_row(task, "failed_missing_output", exit_code=completed.returncode, notes="Generic test-only command exited 0 but did not create its nonempty expected output.")
    return _status_row(task, "completed", exit_code=completed.returncode, notes="Generic non-phylogeny test command completed with a nonempty expected output.")


def run_manifest(
    manifest_path: Path,
    status_dir: Path,
    *,
    workers: int = 1,
    preflight_only: bool = True,
) -> dict[str, int]:
    """Preflight, or unit-test-only run, the immutable Task 3 command manifest."""

    if workers < 1:
        raise ValueError("workers must be at least 1")
    tasks = [_selected_task(job) for job in _load_manifest(Path(manifest_path))]
    keys = [_task_key(task) for task in tasks]
    if len(keys) != len(set(keys)):
        raise ValueError("Task 3 command manifest has duplicate family_category/route/output_path task keys")

    status_by_key: dict[tuple[str, str, str], dict[str, str]] = {}
    pending: list[dict[str, str]] = []
    for task in tasks:
        status, parsed_task = _check_task(task, workers=workers, preflight_only=preflight_only)
        if status is None:
            pending.append(parsed_task)
        else:
            status_by_key[_task_key(task)] = status

    if pending:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_execute_task, task): task for task in pending}
            for future in as_completed(futures):
                task = futures[future]
                status_by_key[_task_key(task)] = future.result()

    rows = [status_by_key[key] for key in sorted(status_by_key)]
    _write_status_atomic(Path(status_dir) / STATUS_FILENAME, rows)
    counts = Counter(row["status"] for row in rows)
    return {status: counts[status] for status in STATUS_VALUES}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Preflight immutable P08 phylogeny command manifests.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--status-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--preflight-only", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if not args.preflight_only:
        parser.error("P08 CLI only permits --preflight-only; execution needs separate authorization")
    try:
        summary = run_manifest(args.manifest, args.status_dir, workers=args.workers, preflight_only=True)
    except ValueError as error:
        parser.error(str(error))
    for status in STATUS_VALUES:
        print(f"{status}: {summary[status]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
