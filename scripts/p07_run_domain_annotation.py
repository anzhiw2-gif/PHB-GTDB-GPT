"""Run P07 annotation commands with preflight and resumable status tracking.

P07 annotates P06 sequence-evidence candidates with domain architecture and
localization signals. This runner records command outcomes only; downstream
review must keep these annotations separate from PHB/PHA degradation phenotype
claims.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path


STATUS_FILENAME = "p07_domain_annotation_run_status.tsv"
REQUIRED_MANIFEST_FIELDS = (
    "tool",
    "fasta_shard",
    "input_fasta",
    "output_path",
    "command",
)
STATUS_FIELDNAMES = (
    "tool",
    "fasta_shard",
    "input_fasta",
    "output_path",
    "executable",
    "status",
    "exit_code",
    "stdout_log_path",
    "stderr_log_path",
    "started_at_utc",
    "finished_at_utc",
    "notes",
)
SUMMARY_KEYS = (
    "preflight_ok",
    "missing_executable",
    "completed",
    "skipped_existing",
    "failed_exit_code",
    "failed_missing_output",
)


def _utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat()


def _load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing = [field for field in REQUIRED_MANIFEST_FIELDS if field not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")

        rows: list[dict[str, str]] = []
        seen: set[tuple[str, str, str]] = set()
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            for field in REQUIRED_MANIFEST_FIELDS:
                if not normalized[field]:
                    raise ValueError(f"{path}:{line_number} is missing required value for {field}")
            key = (normalized["tool"], normalized["fasta_shard"], normalized["output_path"])
            if key in seen:
                raise ValueError(f"{path}:{line_number} has duplicate tool/fasta_shard/output_path {key!r}")
            seen.add(key)
            rows.append(normalized)
    return rows


def _write_status_atomic(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDNAMES, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in STATUS_FIELDNAMES})
    temporary.replace(path)


def _status_row(
    job: dict[str, str],
    *,
    executable: str,
    status: str,
    exit_code: int | None,
    stdout_log_path: Path | None = None,
    stderr_log_path: Path | None = None,
    started_at_utc: str | None = None,
    finished_at_utc: str | None = None,
    notes: str = "",
) -> dict[str, str]:
    started = started_at_utc or _utc_now()
    finished = finished_at_utc or started
    return {
        "tool": job["tool"],
        "fasta_shard": job["fasta_shard"],
        "input_fasta": job["input_fasta"],
        "output_path": job["output_path"],
        "executable": executable,
        "status": status,
        "exit_code": "" if exit_code is None else str(exit_code),
        "stdout_log_path": "" if stdout_log_path is None else stdout_log_path.as_posix(),
        "stderr_log_path": "" if stderr_log_path is None else stderr_log_path.as_posix(),
        "started_at_utc": started,
        "finished_at_utc": finished,
        "notes": notes,
    }


def _completed_checkpoint_keys(status_path: Path) -> set[tuple[str, str, str]]:
    if not status_path.is_file():
        return set()
    with status_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"tool", "fasta_shard", "output_path", "status"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            return set()
        return {
            (row["tool"], row["fasta_shard"], row["output_path"])
            for row in reader
            if row["status"] in {"completed", "skipped_existing"}
        }


def _manifest_key(job: dict[str, str]) -> tuple[str, str, str]:
    return (job["tool"], job["fasta_shard"], job["output_path"])


def _command_executable(command: str) -> str:
    try:
        parts = shlex.split(command, posix=(os.name != "nt"))
    except ValueError as exc:
        raise ValueError(f"Could not parse command executable from {command!r}") from exc
    if not parts:
        raise ValueError("Command is empty")
    return parts[0].strip("\"'")


def _executable_available(executable: str) -> bool:
    path = Path(executable)
    if path.parent != Path("."):
        return path.is_file()
    return shutil.which(executable) is not None


def _output_is_complete(output_path: Path) -> bool:
    if output_path.is_file():
        return output_path.stat().st_size > 0
    if output_path.is_dir():
        return any(path.is_file() and path.stat().st_size > 0 for path in output_path.rglob("*"))
    parent = output_path.parent
    if not parent.is_dir():
        return False
    return any(
        candidate.is_file() and candidate.stat().st_size > 0
        for candidate in parent.glob(f"{output_path.name}*")
    )


def _safe_identifier(value: str) -> str:
    return "".join(char if char.isalnum() or char in {".", "_", "-"} else "_" for char in value) or "unnamed"


def _execute_job(job: dict[str, str], log_dir: Path) -> dict[str, str]:
    executable = _command_executable(job["command"])
    if not _executable_available(executable):
        return _status_row(
            job,
            executable=executable,
            status="missing_executable",
            exit_code=None,
            notes="Executable was not found before command execution.",
        )

    output_path = Path(job["output_path"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_stem = f"{_safe_identifier(job['tool'])}__{_safe_identifier(job['fasta_shard'])}"
    stdout_log_path = log_dir / f"{log_stem}.stdout.log"
    stderr_log_path = log_dir / f"{log_stem}.stderr.log"
    started_at_utc = _utc_now()
    run_kwargs: dict[str, object] = {"shell": True}
    if os.name == "posix":
        run_kwargs["executable"] = "/bin/bash"
    with stdout_log_path.open("w", encoding="utf-8") as stdout_handle:
        with stderr_log_path.open("w", encoding="utf-8") as stderr_handle:
            completed = subprocess.run(
                job["command"],
                stdout=stdout_handle,
                stderr=stderr_handle,
                **run_kwargs,
            )
    finished_at_utc = _utc_now()
    if completed.returncode != 0:
        return _status_row(
            job,
            executable=executable,
            status="failed_exit_code",
            exit_code=completed.returncode,
            stdout_log_path=stdout_log_path,
            stderr_log_path=stderr_log_path,
            started_at_utc=started_at_utc,
            finished_at_utc=finished_at_utc,
        )
    if not _output_is_complete(output_path):
        return _status_row(
            job,
            executable=executable,
            status="failed_missing_output",
            exit_code=completed.returncode,
            stdout_log_path=stdout_log_path,
            stderr_log_path=stderr_log_path,
            started_at_utc=started_at_utc,
            finished_at_utc=finished_at_utc,
            notes="Command exited 0, but expected output path or prefix has no nonempty files.",
        )
    return _status_row(
        job,
        executable=executable,
        status="completed",
        exit_code=completed.returncode,
        stdout_log_path=stdout_log_path,
        stderr_log_path=stderr_log_path,
        started_at_utc=started_at_utc,
        finished_at_utc=finished_at_utc,
    )


def run_manifest(
    manifest_path: Path,
    status_dir: Path,
    *,
    workers: int = 1,
    preflight_only: bool = False,
) -> dict[str, int]:
    """Run or preflight P07 annotation commands and write atomic status rows."""

    if workers < 1:
        raise ValueError("workers must be at least 1")
    jobs = _load_manifest(manifest_path)
    status_path = status_dir / STATUS_FILENAME
    log_dir = status_dir / "logs"
    completed_checkpoints = _completed_checkpoint_keys(status_path)
    status_by_index: dict[int, dict[str, str]] = {}
    pending: list[tuple[int, dict[str, str]]] = []

    for index, job in enumerate(jobs):
        executable = _command_executable(job["command"])
        if not _executable_available(executable):
            status_by_index[index] = _status_row(
                job,
                executable=executable,
                status="missing_executable",
                exit_code=None,
                notes="Executable was not found during preflight.",
            )
            continue
        if preflight_only:
            status_by_index[index] = _status_row(
                job,
                executable=executable,
                status="preflight_ok",
                exit_code=None,
            )
            continue
        if _manifest_key(job) in completed_checkpoints and _output_is_complete(Path(job["output_path"])):
            status_by_index[index] = _status_row(
                job,
                executable=executable,
                status="skipped_existing",
                exit_code=None,
            )
        else:
            pending.append((index, job))

    _write_status_atomic(status_path, [status_by_index[index] for index in sorted(status_by_index)])

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(_execute_job, job, log_dir): index
            for index, job in pending
        }
        for future in as_completed(futures):
            status_by_index[futures[future]] = future.result()
            _write_status_atomic(status_path, [status_by_index[index] for index in sorted(status_by_index)])

    counts = Counter(row["status"] for row in status_by_index.values())
    return {name: counts[name] for name in SUMMARY_KEYS}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run P07 InterPro/localization command manifest with checkpoints.")
    parser.add_argument("--manifest", type=Path, default=Path("06_domain_annotation/manifests/p07_domain_annotation_command_manifest.tsv"))
    parser.add_argument("--status-dir", type=Path, default=Path("06_domain_annotation/run_status"))
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--preflight-only", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = run_manifest(
        args.manifest,
        args.status_dir,
        workers=args.workers,
        preflight_only=args.preflight_only,
    )
    for name in SUMMARY_KEYS:
        print(f"{name}: {summary[name]}")
    return 1 if any(summary[name] for name in ("missing_executable", "failed_exit_code", "failed_missing_output")) else 0


if __name__ == "__main__":
    raise SystemExit(main())
