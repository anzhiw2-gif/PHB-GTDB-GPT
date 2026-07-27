"""Stream P06 proteomes while recording HMMER-length-limited targets.

HMMER 3.4 cannot compare protein targets longer than 100,000 residues with the
standard comparison pipeline. Those records are retained in an explicit audit
table rather than being interpreted as non-hits.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from pathlib import Path
from typing import Iterable, TextIO


EXCLUSION_FIELDNAMES = (
    "source_proteome_path",
    "target_id",
    "sequence_length_aa",
    "reason",
)
OVERLONG_REASON = "hmmsearch_target_length_gt_100000"


def _open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _iter_fasta_records(path: Path) -> Iterable[tuple[str, str]]:
    header: str | None = None
    sequence_lines: list[str] = []
    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_lines)
                header = line[1:]
                sequence_lines = []
            elif header is None:
                raise ValueError(f"{path} contains sequence data before its first FASTA header")
            else:
                sequence_lines.append(line)
    if header is not None:
        yield header, "".join(sequence_lines)


def stream_proteomes(
    proteome_paths: Iterable[Path],
    *,
    max_protein_length: int,
    exclusion_path: Path,
    output: TextIO,
) -> dict[str, int]:
    """Write eligible FASTA records to stdout and atomically audit overlong targets."""

    if max_protein_length < 1:
        raise ValueError("max_protein_length must be at least 1")
    exclusion_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = exclusion_path.with_name(f".{exclusion_path.name}.{os.getpid()}.tmp")
    emitted_sequences = 0
    excluded_overlong_sequences = 0
    try:
        with temporary.open("w", encoding="utf-8", newline="") as exclusion_handle:
            writer = csv.DictWriter(exclusion_handle, fieldnames=EXCLUSION_FIELDNAMES, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for path in proteome_paths:
                if not path.is_file():
                    raise FileNotFoundError(f"Proteome FASTA is missing: {path}")
                for header, sequence in _iter_fasta_records(path):
                    if len(sequence) > max_protein_length:
                        writer.writerow(
                            {
                                "source_proteome_path": path.as_posix(),
                                "target_id": header.split(maxsplit=1)[0],
                                "sequence_length_aa": str(len(sequence)),
                                "reason": OVERLONG_REASON,
                            }
                        )
                        excluded_overlong_sequences += 1
                        continue
                    output.write(f">{header}\n{sequence}\n")
                    emitted_sequences += 1
        temporary.replace(exclusion_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return {
        "emitted_sequences": emitted_sequences,
        "excluded_overlong_sequences": excluded_overlong_sequences,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Stream P06 proteomes after excluding HMMER-length-limited proteins.")
    parser.add_argument("proteome_paths", nargs="+", type=Path)
    parser.add_argument("--max-protein-length", type=int, default=100000)
    parser.add_argument("--exclusion-path", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    summary = stream_proteomes(
        args.proteome_paths,
        max_protein_length=args.max_protein_length,
        exclusion_path=args.exclusion_path,
        output=sys.stdout,
    )
    print(
        f"emitted_sequences={summary['emitted_sequences']} excluded_overlong_sequences={summary['excluded_overlong_sequences']}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
