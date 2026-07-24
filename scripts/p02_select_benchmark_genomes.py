"""Select deterministic P02 benchmark genomes from copied GTDB inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from typing import Iterable, Sequence


FASTA_SUFFIXES = (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")


def select_benchmark_genomes(
    genomes_root: Path,
    taxonomy_path: Path | Sequence[Path],
    sample_size: int = 240,
    seed: int = 20260724,
) -> list[dict[str, object]]:
    taxonomy = load_taxonomy(taxonomy_path)
    genomes = discover_genomes(genomes_root)
    rows: list[dict[str, object]] = []

    for accession, genome_path in genomes.items():
        lineage = taxonomy.get(accession)
        if lineage is None:
            continue
        rows.append(
            {
                "accession": accession,
                "domain": taxonomy_rank(lineage, "d"),
                "stratum": stratum_for_lineage(lineage),
                "genome_path": str(genome_path),
                "genome_bytes": genome_path.stat().st_size,
            }
        )

    selected = stratified_sample(rows, sample_size=sample_size, seed=seed)
    return sorted(selected, key=lambda row: str(row["accession"]))


def load_taxonomy(taxonomy_path: Path | Sequence[Path]) -> dict[str, str]:
    taxonomy: dict[str, str] = {}
    paths = [taxonomy_path] if isinstance(taxonomy_path, Path) else list(taxonomy_path)
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                accession, lineage = stripped.split("\t", 1)
                taxonomy[accession] = lineage
                taxonomy[strip_gtdb_prefix(accession)] = lineage
    return taxonomy


def discover_genomes(genomes_root: Path) -> dict[str, Path]:
    genomes: dict[str, Path] = {}
    for path in sorted(p for p in genomes_root.rglob("*") if p.is_file() and is_fasta(p)):
        accession = infer_accession(path, genomes_root)
        if accession not in genomes:
            genomes[accession] = path
    return genomes


def is_fasta(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def infer_accession(path: Path, genomes_root: Path) -> str:
    name = path.name
    for suffix in FASTA_SUFFIXES:
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_genomic"):
        return name[: -len("_genomic")]
    return name


def strip_gtdb_prefix(accession: str) -> str:
    if accession.startswith(("RS_", "GB_")):
        return accession[3:]
    return accession


def taxonomy_rank(lineage: str, prefix: str) -> str:
    wanted = prefix + "__"
    for part in lineage.split(";"):
        if part.startswith(wanted):
            return part
    return f"{prefix}__unknown"


def stratum_for_lineage(lineage: str) -> str:
    domain = taxonomy_rank(lineage, "d")
    if domain == "d__Archaea":
        return f"{domain}|{taxonomy_rank(lineage, 'c')}"
    return f"{domain}|{taxonomy_rank(lineage, 'p')}"


def stratified_sample(rows: list[dict[str, object]], sample_size: int, seed: int) -> list[dict[str, object]]:
    strata: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        strata.setdefault(str(row["stratum"]), []).append(row)

    selected: list[dict[str, object]] = []
    for stratum in sorted(strata):
        candidates = sorted(
            strata[stratum],
            key=lambda row: str(row["accession"]),
        )
        selected.append(candidates[0])
        if len(selected) == sample_size:
            return selected

    remaining = [row for row in rows if row not in selected]
    remaining.sort(key=lambda row: stable_score(str(row["accession"]), seed))
    selected.extend(remaining[: max(0, sample_size - len(selected))])
    return selected


def stable_score(value: str, seed: int) -> str:
    return hashlib.sha256(f"{seed}|{value}".encode("utf-8")).hexdigest()


def write_selection(path: Path, rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession", "domain", "stratum", "genome_path", "genome_bytes"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Select deterministic P02 benchmark genomes.")
    parser.add_argument("--genomes-root", required=True, type=Path)
    parser.add_argument("--taxonomy", required=True, type=Path, action="append")
    parser.add_argument("--sample-size", type=int, default=240)
    parser.add_argument("--seed", type=int, default=20260724)
    parser.add_argument("--out", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    rows = select_benchmark_genomes(
        genomes_root=args.genomes_root,
        taxonomy_path=args.taxonomy,
        sample_size=args.sample_size,
        seed=args.seed,
    )
    write_selection(args.out, rows)
    print(f"Selected {len(rows)} benchmark genomes: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
