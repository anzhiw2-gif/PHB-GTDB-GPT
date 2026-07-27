"""Execute checksum-locked P06 HMMER jobs with resumable status tracking.

The P06 planner owns model selection, checksums, and calibrated thresholds.
This runner only executes its manifest commands, preserving completed nonempty
``domtblout`` files and recording each invocation outcome.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path


STATUS_FILENAME = "p06_hmmer_run_status.tsv"
REQUIRED_MANIFEST_FIELDS = ("family_category", "proteome_shard", "domtblout_path", "command")
STATUS_FIELDNAMES = (
    "family_category",
    "proteome_shard",
    "domtblout_path",
    "status",
    "exit_code",
    "stderr_log_path",
    "started_at_utc",
    "finished_at_utc",
    "command",
)


def _utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat()


def _safe_identifier(value: str) -> str:
    return "".join(char if char.isalnum() or char in {".", "_", "-"} else "_" for char in value) or "unnamed"


def _load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is missing a tabular header")
        missing = [field for field in REQUIRED_MANIFEST_FIELDS if field not in reader.fieldnames]
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")

        rows: list[dict[str, str]] = []
        seen_pairs: set[tuple[str, str]] = set()
        for line_number, row in enumerate(reader, start=2):
            normalized = {key: (value or "").strip() for key, value in row.items()}
            for field in REQUIRED_MANIFEST_FIELDS:
                if not normalized[field]:
                    raise ValueError(f"{path}:{line_number} is missing required value for {field}")
            pair = (normalized["family_category"], normalized["proteome_shard"])
            if pair in seen_pairs:
                raise ValueError(f"{path}:{line_number} has duplicate family_category/proteome_shard {pair!r}")
            seen_pairs.add(pair)
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
    status: str,
    exit_code: int | None,
    stderr_log_path: Path | None,
    started_at_utc: str,
    finished_at_utc: str,
) -> dict[str, str]:
    return {
        "family_category": job["family_category"],
        "proteome_shard": job["proteome_shard"],
        "domtblout_path": job["domtblout_path"],
        "status": status,
        "exit_code": "" if exit_code is None else str(exit_code),
        "stderr_log_path": "" if stderr_log_path is None else stderr_log_path.as_posix(),
        "started_at_utc": started_at_utc,
        "finished_at_utc": finished_at_utc,
        "command": job["command"],
    }


def _execute_job(job: dict[str, str], log_dir: Path) -> dict[str, str]:
    domtblout_path = Path(job["domtblout_path"])
    domtblout_path.parent.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_name = f"{_safe_identifier(job['family_category'])}__{_safe_identifier(job['proteome_shard'])}.stderr.log"
    stderr_log_path = log_dir / log_name
    started_at_utc = _utc_now()
    run_kwargs: dict[str, object] = {"shell": True}
    if os.name == "posix":
        run_kwargs["executable"] = "/bin/bash"
    with stderr_log_path.open("w", encoding="utf-8") as stderr_handle:
        completed = subprocess.run(job["command"], stdout=subprocess.DEVNULL, stderr=stderr_handle, **run_kwargs)
    finished_at_utc = _utc_now()
    if completed.returncode != 0:
        return _status_row(
            job,
            status="failed_exit_code",
            exit_code=completed.returncode,
            stderr_log_path=stderr_log_path,
            started_at_utc=started_at_utc,
            finished_at_utc=finished_at_utc,
        )
    if not domtblout_path.is_file() or domtblout_path.stat().st_size == 0:
        return _status_row(
            job,
            status="failed_empty_domtblout",
            exit_code=completed.returncode,
            stderr_log_path=stderr_log_path,
            started_at_utc=started_at_utc,
            finished_at_utc=finished_at_utc,
        )
    return _status_row(
        job,
        status="completed",
        exit_code=completed.returncode,
        stderr_log_path=stderr_log_path,
        started_at_utc=started_at_utc,
        finished_at_utc=finished_at_utc,
    )


def run_manifest(manifest_path: Path, status_dir: Path, *, workers: int = 1) -> dict[str, int]:
    """Run missing P06 jobs and write an atomic status snapshot after each result."""

    if workers < 1:
        raise ValueError("workers must be at least 1")
    jobs = _load_manifest(manifest_path)
    status_path = status_dir / STATUS_FILENAME
    stderr_log_dir = status_dir / "stderr_logs"
    status_by_index: dict[int, dict[str, str]] = {}
    pending: list[tuple[int, dict[str, str]]] = []
    for index, job in enumerate(jobs):
        domtblout_path = Path(job["domtblout_path"])
        if domtblout_path.is_file() and domtblout_path.stat().st_size > 0:
            timestamp = _utc_now()
            status_by_index[index] = _status_row(
                job,
                status="skipped_existing",
                exit_code=None,
                stderr_log_path=None,
                started_at_utc=timestamp,
                finished_at_utc=timestamp,
            )
        else:
            pending.append((index, job))
    _write_status_atomic(status_path, [status_by_index[index] for index in sorted(status_by_index)])

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(_execute_job, job, stderr_log_dir): index
            for index, job in pending
        }
        for future in as_completed(futures):
            status_by_index[futures[future]] = future.result()
            _write_status_atomic(status_path, [status_by_index[index] for index in sorted(status_by_index)])

    counts = Counter(row["status"] for row in status_by_index.values())
    return {
        "completed": counts["completed"],
        "skipped_existing": counts["skipped_existing"],
        "failed_exit_code": counts["failed_exit_code"],
        "failed_empty_domtblout": counts["failed_empty_domtblout"],
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run P06 HMMER manifest commands with resumable checkpoints.")
    parser.add_argument("--manifest", type=Path, default=Path("05_hmmer_scan/p06_hmmer_scan_manifest.tsv"))
    parser.add_argument("--status-dir", type=Path, default=Path("05_hmmer_scan/run_status"))
    parser.add_argument("--workers", type=int, default=1)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = run_manifest(args.manifest, args.status_dir, workers=args.workers)
    for name in ("completed", "skipped_existing", "failed_exit_code", "failed_empty_domtblout"):
        print(f"{name}: {summary[name]}")
    return 1 if summary["failed_exit_code"] or summary["failed_empty_domtblout"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
