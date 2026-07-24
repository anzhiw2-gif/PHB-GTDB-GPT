"""Realtime terminal monitor for P03 GTDB proteome prediction."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path


DEFAULT_ROOT = Path("/home/data/haoyu/PHB-GTDB-GPT")
DEFAULT_TOTAL_GENOMES = 199_923


def count_outputs(output_root: Path) -> tuple[int, int]:
    faa_count = count_files(output_root / "faa", ".faa.gz")
    gff_count = count_files(output_root / "gff", ".gff.gz")
    return faa_count, gff_count


def count_files(root: Path, suffix: str) -> int:
    if not root.exists():
        return 0

    count = 0
    for _, _, filenames in os.walk(root):
        count += sum(1 for filename in filenames if filename.endswith(suffix))
    return count


def detect_p03_pid() -> int | None:
    try:
        result = subprocess.run(
            ["pgrep", "-f", "scripts/p03_predict_proteomes.py"],
            check=False,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        return None

    pids = [int(line) for line in result.stdout.splitlines() if line.strip().isdigit()]
    current_pid = os.getpid()
    pids = [pid for pid in pids if pid != current_pid]
    return pids[-1] if pids else None


def process_status(pid: int | None) -> dict[str, object]:
    if pid is None:
        return {"running": False, "pid": None}

    result = subprocess.run(
        ["ps", "-p", str(pid), "-o", "etimes=,%cpu=,%mem=,stat="],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return {"running": False, "pid": pid}

    parts = result.stdout.split()
    return {
        "running": True,
        "pid": pid,
        "elapsed_seconds": int(parts[0]),
        "cpu_percent": parts[1] if len(parts) > 1 else "?",
        "mem_percent": parts[2] if len(parts) > 2 else "?",
        "state": parts[3] if len(parts) > 3 else "?",
    }


def format_duration(seconds: float | int | None) -> str:
    if seconds is None:
        return "unknown"
    return str(timedelta(seconds=max(0, int(seconds))))


def progress_bar(done: int, total: int, width: int) -> str:
    if total <= 0:
        return "[" + "?" * width + "]"
    ratio = min(max(done / total, 0.0), 1.0)
    filled = int(width * ratio)
    return "[" + "#" * filled + "-" * (width - filled) + "]"


def latest_monitor_lines(project_root: Path, limit: int = 5) -> list[str]:
    monitor_path = project_root / "03_gtdb_proteomes" / "qc" / "p03_monitor.log"
    if not monitor_path.exists():
        return []
    try:
        with monitor_path.open("r", encoding="utf-8", errors="replace") as handle:
            return [line.rstrip("\n") for line in handle.readlines()[-limit:]]
    except OSError:
        return []


def clear_screen() -> None:
    sys.stdout.write("\033[2J\033[H")
    sys.stdout.flush()


def render(project_root: Path, total_genomes: int, pid: int | None, no_clear: bool = False) -> None:
    output_root = project_root / "03_gtdb_proteomes"
    faa_count, gff_count = count_outputs(output_root)
    completed = min(faa_count, gff_count)
    pair_gap = abs(faa_count - gff_count)

    status = process_status(pid)
    if not status.get("running"):
        auto_pid = detect_p03_pid()
        if auto_pid != pid:
            status = process_status(auto_pid)
            pid = auto_pid

    elapsed_seconds = status.get("elapsed_seconds")
    rate = None
    if isinstance(elapsed_seconds, int) and elapsed_seconds > 0:
        rate = completed / elapsed_seconds

    remaining = max(total_genomes - completed, 0)
    eta_seconds = remaining / rate if rate else None
    finish_at = datetime.now() + timedelta(seconds=eta_seconds) if eta_seconds else None
    percent = (completed / total_genomes * 100) if total_genomes else 0.0

    columns = shutil.get_terminal_size((100, 30)).columns
    bar_width = max(20, min(60, columns - 32))

    if not no_clear:
        clear_screen()

    print("P03 GTDB Proteome Prediction Monitor")
    print("=" * min(columns, 90))
    print(f"time             : {datetime.now().strftime('%F %T')}")
    print(f"project_root     : {project_root}")
    print(f"status           : {'RUNNING' if status.get('running') else 'STOPPED'}")
    print(f"pid              : {status.get('pid') or 'not found'}")
    print(
        "cpu / mem / state: "
        f"{status.get('cpu_percent', '?')}% / {status.get('mem_percent', '?')}% / {status.get('state', '?')}"
    )
    print(f"job elapsed      : {format_duration(elapsed_seconds if isinstance(elapsed_seconds, int) else None)}")
    print()
    print(f"FAA files        : {faa_count:,}")
    print(f"GFF files        : {gff_count:,}")
    print(f"completed genomes: {completed:,} / {total_genomes:,}")
    print(f"file pair gap    : {pair_gap:,}")
    print()
    print(f"{progress_bar(completed, total_genomes, bar_width)} {percent:6.2f}%")
    print()
    if rate:
        print(f"speed            : {rate:,.2f} genomes/sec  ({rate * 3600:,.0f} genomes/hour)")
    else:
        print("speed            : unknown")
    print(f"remaining        : {remaining:,}")
    print(f"ETA              : {format_duration(eta_seconds)}")
    print(f"estimated finish : {finish_at.strftime('%F %T') if finish_at else 'unknown'}")

    lines = latest_monitor_lines(project_root)
    if lines:
        print()
        print("latest monitor log")
        print("-" * min(columns, 90))
        for line in lines:
            print(line)

    print()
    print("Ctrl+C stops this monitor only. It does not stop the P03 job.")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Monitor the P03 prediction run in the terminal.")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="PHB-GTDB-GPT project root on T141")
    parser.add_argument("--total", type=int, default=DEFAULT_TOTAL_GENOMES, help="Total expected genome count")
    parser.add_argument("--pid", type=int, default=None, help="P03 process PID; auto-detected if omitted")
    parser.add_argument("--interval", type=int, default=30, help="Refresh interval in seconds")
    parser.add_argument("--once", action="store_true", help="Print one snapshot and exit")
    parser.add_argument("--no-clear", action="store_true", help="Do not clear the terminal before each refresh")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.once:
        render(args.root, args.total, args.pid, no_clear=True)
        return 0

    try:
        while True:
            render(args.root, args.total, args.pid, no_clear=args.no_clear)
            time.sleep(max(1, args.interval))
    except KeyboardInterrupt:
        print("\nMonitor stopped.")
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
