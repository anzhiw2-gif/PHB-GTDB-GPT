"""Monitor the P03 translation-fix rerun on T141.

The original P03 run can look complete by file count while still containing
empty protein sequences. This monitor therefore tracks rewrite progress from
the fix-run pidfile and samples stable FASTA files for residue counts.
"""

from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path


DEFAULT_ROOT = Path("/home/data/haoyu/PHB-GTDB-GPT")
DEFAULT_TOTAL_GENOMES = 199_923
DEFAULT_PID_GLOB = "p03_translation_fix_*.pid"


def latest_pidfile(run_log_dir: Path, pattern: str = DEFAULT_PID_GLOB) -> Path | None:
    pidfiles = sorted(
        (path for path in run_log_dir.glob(pattern) if path.is_file()),
        key=lambda path: path.stat().st_mtime,
    )
    return pidfiles[-1] if pidfiles else None


def read_pid(pidfile: Path | None) -> int | None:
    if pidfile is None:
        return None
    try:
        text = pidfile.read_text(encoding="utf-8").strip()
        return int(text) if text else None
    except (OSError, ValueError):
        return None


def process_status(pid: int | None) -> dict[str, object]:
    if pid is None:
        return {"running": False, "pid": None}
    result = subprocess.run(
        ["ps", "-p", str(pid), "-o", "pid=,ppid=,stat=,etimes=,%cpu=,%mem=,cmd="],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return {"running": False, "pid": pid}
    parts = result.stdout.strip().split(maxsplit=6)
    return {
        "running": True,
        "pid": int(parts[0]),
        "ppid": int(parts[1]),
        "state": parts[2],
        "elapsed_seconds": int(parts[3]),
        "cpu_percent": parts[4],
        "mem_percent": parts[5],
        "command": parts[6] if len(parts) > 6 else "",
    }


def inventory_faa(faa_root: Path, run_start: float) -> dict[str, object]:
    total = older = newer = stable_new = 0
    newest_mtime = 0.0
    newest_path = ""
    now = time.time()
    for path in iter_faa_files(faa_root):
        total += 1
        mtime = path.stat().st_mtime
        if mtime < run_start:
            older += 1
        else:
            newer += 1
            if now - mtime >= 120:
                stable_new += 1
        if mtime > newest_mtime:
            newest_mtime = mtime
            newest_path = str(path)
    return {
        "total": total,
        "older": older,
        "newer": newer,
        "stable_new": stable_new,
        "newest_path": newest_path,
        "newest_mtime": newest_mtime,
    }


def iter_faa_files(faa_root: Path):
    if not faa_root.is_dir():
        return
    for current_root, _, filenames in os.walk(faa_root):
        for filename in filenames:
            if filename.endswith(".faa.gz"):
                yield Path(current_root) / filename


def fasta_stats(path: Path) -> tuple[int, int, str]:
    seqs = 0
    residues = 0
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith(">"):
                    seqs += 1
                else:
                    residues += len("".join(stripped.split()))
    except Exception as error:  # pragma: no cover - depends on live filesystem state
        return seqs, residues, f"{type(error).__name__}: {error}"
    return seqs, residues, ""


def sample_files(faa_root: Path, *, run_start: float, newer: bool, limit: int, stable_age: int) -> list[Path]:
    now = time.time()
    selected: list[Path] = []
    for path in iter_faa_files(faa_root):
        mtime = path.stat().st_mtime
        if newer:
            if mtime >= run_start and now - mtime >= stable_age:
                selected.append(path)
        elif mtime < run_start:
            selected.append(path)
        if len(selected) >= limit:
            break
    return selected


def tail_lines(path: Path, limit: int) -> list[str]:
    if not path.is_file():
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return [line.rstrip("\n") for line in handle.readlines()[-limit:]]
    except OSError:
        return []


def format_duration(seconds: float | int | None) -> str:
    if seconds is None:
        return "unknown"
    return str(timedelta(seconds=max(0, int(seconds))))


def format_time(timestamp: float | None) -> str:
    if not timestamp:
        return "unknown"
    return datetime.fromtimestamp(timestamp).strftime("%F %T")


def progress_bar(done: int, total: int, width: int = 42) -> str:
    if total <= 0:
        return "[" + "?" * width + "]"
    ratio = min(max(done / total, 0.0), 1.0)
    filled = int(width * ratio)
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def render(args: argparse.Namespace) -> None:
    root = args.root
    output_root = root / "03_gtdb_proteomes"
    faa_root = output_root / "faa"
    run_log_dir = output_root / "run_logs"
    pidfile = args.pid_file or latest_pidfile(run_log_dir)
    run_start = pidfile.stat().st_mtime if pidfile and pidfile.exists() else time.time()
    status = process_status(read_pid(pidfile))
    inv = inventory_faa(faa_root, run_start)
    rewritten = int(inv["newer"])
    total_expected = args.total
    elapsed = time.time() - run_start
    rate = rewritten / elapsed if elapsed > 0 and rewritten else 0.0
    remaining = max(total_expected - rewritten, 0)
    eta = remaining / rate if rate else None

    print("P03 translation-fix monitor")
    print("=" * 78)
    print(f"time                 : {datetime.now().strftime('%F %T')}")
    print(f"project_root         : {root}")
    print(f"pidfile              : {pidfile or 'not found'}")
    print(f"fix_start            : {format_time(run_start)}")
    print(f"process              : {'RUNNING' if status.get('running') else 'STOPPED'} pid={status.get('pid')}")
    print(
        "cpu/mem/state        : "
        f"{status.get('cpu_percent', '?')} / {status.get('mem_percent', '?')} / {status.get('state', '?')}"
    )
    print(f"process_elapsed      : {format_duration(status.get('elapsed_seconds') if status.get('running') else None)}")
    print()
    print(f"main FAA total       : {inv['total']:,}")
    print(f"rewritten since start: {rewritten:,} / {total_expected:,}")
    print(f"still old by mtime   : {inv['older']:,}")
    print(f"stable rewritten     : {inv['stable_new']:,}  (mtime older than 120s)")
    print(f"newest FAA           : {inv['newest_path'] or 'none'}")
    print(f"newest mtime         : {format_time(float(inv['newest_mtime']) if inv['newest_mtime'] else None)}")
    print(f"{progress_bar(rewritten, total_expected)} {rewritten / total_expected * 100 if total_expected else 0:6.2f}%")
    print(f"rewrite speed        : {rate:,.2f} genomes/sec ({rate * 3600:,.0f} genomes/hour)" if rate else "rewrite speed        : unknown")
    print(f"ETA                  : {format_duration(eta)}")
    print()
    print("old untouched FAA sample")
    for path in sample_files(faa_root, run_start=run_start, newer=False, limit=args.sample, stable_age=args.stable_age):
        seqs, residues, error = fasta_stats(path)
        print(f"  {path.relative_to(root)} seqs={seqs} residues={residues} error={error or '-'}")
    print("stable rewritten FAA sample")
    for path in sample_files(faa_root, run_start=run_start, newer=True, limit=args.sample, stable_age=args.stable_age):
        seqs, residues, error = fasta_stats(path)
        print(f"  {path.relative_to(root)} seqs={seqs} residues={residues} error={error or '-'}")
    print()
    print("stderr tail")
    stderr_path = run_log_dir / "p03_translation_fix_20260726.stderr"
    for line in tail_lines(stderr_path, args.tail):
        print(f"  {line}")
    print()
    print("Note: p03_prediction_qc.tsv is rewritten at the end of P03, so it may still show old mean_protein_length=0 while the rerun is active.")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Monitor P03 translation-fix progress and sample FAA residue counts.")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="PHB-GTDB-GPT root on T141")
    parser.add_argument("--pid-file", type=Path, default=None, help="Explicit P03 translation-fix pidfile")
    parser.add_argument("--total", type=int, default=DEFAULT_TOTAL_GENOMES, help="Expected total genome count")
    parser.add_argument("--interval", type=int, default=60, help="Refresh interval for watch mode")
    parser.add_argument("--once", action="store_true", help="Print one snapshot and exit")
    parser.add_argument("--sample", type=int, default=5, help="Number of old/new FAA files to sample")
    parser.add_argument("--stable-age", type=int, default=120, help="Only sample rewritten files older than this many seconds")
    parser.add_argument("--tail", type=int, default=12, help="stderr lines to print")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.once:
        render(args)
        return 0
    try:
        while True:
            print("\033[2J\033[H", end="")
            render(args)
            sys.stdout.flush()
            time.sleep(max(1, args.interval))
    except KeyboardInterrupt:
        print("\nMonitor stopped. The P03 job is not interrupted.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
