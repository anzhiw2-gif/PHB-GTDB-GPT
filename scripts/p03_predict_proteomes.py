"""Predict complete GTDB proteomes with the locked P03 Pyrodigal route."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from statistics import mean
from typing import Callable, Iterable, Sequence


KNOWN_FASTA_SUFFIXES = (
    ".fna.gz",
    ".fa.gz",
    ".fasta.gz",
    ".fna",
    ".fa",
    ".fasta",
)
ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")
DEFAULT_MAX_THREADS = 60
DEFAULT_GZIP_LEVEL = 0

try:  # pragma: no cover - optional runtime accelerator on T141
    from isal import igzip as gzip_backend
except ImportError:  # pragma: no cover - fallback for local test environments
    gzip_backend = gzip

Predictor = Callable[[Path, str, str], dict[str, object]]


def infer_accession_from_genome_path(genome_path: Path) -> str:
    """Infer a GTDB accession from a canonical GTDB genome filename."""

    name = genome_path.name
    stripped = True
    while stripped:
        stripped = False
        for suffix in KNOWN_FASTA_SUFFIXES:
            if name.lower().endswith(suffix):
                name = name[: -len(suffix)]
                stripped = True
                break
    if name.endswith("_genomic"):
        name = name[: -len("_genomic")]
    return name


def discover_genome_files(genomes_root: Path) -> list[Path]:
    """Return a deterministic list of GTDB genome FASTA files under a root."""

    genomes_root = genomes_root.resolve(strict=False)
    discovered: dict[str, Path] = {}
    for path in sorted(
        (candidate for candidate in genomes_root.rglob("*") if candidate.is_file()),
        key=lambda candidate: candidate.as_posix(),
    ):
        if not _is_genome_fasta(path):
            continue
        accession = infer_accession_from_genome_path(path)
        discovered.setdefault(accession, path)
    return list(discovered.values())


def load_prediction_policy(policy_path: Path) -> dict[str, str]:
    """Load the locked P03 prediction policy from a small YAML-like file."""

    policy = _parse_prediction_policy(policy_path)
    if policy.get("selected_predictor") != "pyrodigal":
        raise ValueError("P03 policy must select pyrodigal")
    if policy.get("selected_mode") != "meta":
        raise ValueError("P03 policy must select meta mode")
    return policy


def output_paths_for_genome(genome_path: Path, genomes_root: Path, output_root: Path) -> tuple[Path, Path]:
    """Mirror the input genome tree under FAA and GFF output roots."""

    relative = _relative_path(genome_path, genomes_root)
    accession = infer_accession_from_genome_path(genome_path)
    faa_path = output_root / "faa" / relative.parent / f"{accession}.faa.gz"
    gff_path = output_root / "gff" / relative.parent / f"{accession}.gff.gz"
    return faa_path, gff_path


def predict_genome(
    genome_path: Path,
    policy: dict[str, str],
    predictor: Predictor | None = None,
) -> dict[str, object]:
    """Predict one genome and return the in-memory genome result."""

    accession = infer_accession_from_genome_path(genome_path)
    mode = policy.get("selected_mode")
    if policy.get("selected_predictor") != "pyrodigal":
        raise ValueError("P03 only supports pyrodigal")
    if mode != "meta":
        raise ValueError("P03 only supports meta mode")

    predictor = predictor or run_pyrodigal_prediction
    result = predictor(genome_path, accession, mode)
    result.setdefault("accession", accession)
    result.setdefault("genome_path", str(genome_path))
    result.setdefault("status", "ok")
    if "summary" not in result:
        result["summary"] = summarize_genome_result(result)
    return result


def run_p03_prediction(
    genome_paths: Iterable[Path],
    policy: dict[str, str],
    output_root: Path,
    threads: int = DEFAULT_MAX_THREADS,
    predictor: Predictor | None = None,
    genomes_root: Path | None = None,
) -> dict[str, object]:
    """Predict a batch of genomes, write outputs, and return a run summary."""

    predictor = predictor or run_pyrodigal_prediction
    paths = list(genome_paths)
    if not paths:
        return {
            "policy": policy,
            "genome_count": 0,
            "genomes": [],
            "summary": {
                "ok_genomes": 0,
                "failed_genomes": 0,
                "predicted_genes_total": 0,
                "mean_coding_density": 0.0,
            },
        }

    genomes_root = genomes_root or _infer_genomes_root(paths)
    worker_count = min(max(1, threads), DEFAULT_MAX_THREADS, os.cpu_count() or 1, len(paths))

    def process_one(genome_path: Path) -> dict[str, object]:
        try:
            genome_result = predict_genome(genome_path, policy, predictor=predictor)
            faa_path, gff_path = output_paths_for_genome(genome_path, genomes_root, output_root)
            write_prediction_outputs(genome_result, faa_path, gff_path)
            qc_row = _qc_row_from_result(genome_result, faa_path, gff_path, status="ok")
            qc_row["error"] = ""
            return qc_row
        except Exception as error:  # pragma: no cover - exercised in smoke/robustness runs
            faa_path, gff_path = output_paths_for_genome(genome_path, genomes_root, output_root)
            return _failed_qc_row(genome_path, faa_path, gff_path, error)

    if worker_count == 1:
        genomes = [process_one(path) for path in paths]
    else:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            genomes = list(executor.map(process_one, paths))

    qc_path = output_root / "qc" / "p03_prediction_qc.tsv"
    manifest_path = output_root / "manifests" / "p03_prediction_manifest.tsv"
    write_tabular_rows(qc_path, genomes, _qc_fieldnames())
    write_tabular_rows(manifest_path, genomes, _manifest_fieldnames())

    return {
        "policy": policy,
        "genome_count": len(genomes),
        "genomes": genomes,
        "summary": summarize_run(genomes),
        "qc_path": qc_path,
        "manifest_path": manifest_path,
    }


def write_prediction_outputs(genome_result: dict[str, object], faa_path: Path, gff_path: Path) -> None:
    """Write gzipped FAA and GFF3 outputs for one predicted genome."""

    accession = str(genome_result["accession"])
    contigs = list(genome_result.get("contigs", []))
    faa_path.parent.mkdir(parents=True, exist_ok=True)
    gff_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip_backend.open(faa_path, "wt", encoding="utf-8", newline="\n", compresslevel=DEFAULT_GZIP_LEVEL) as faa_handle:
        for contig in contigs:
            contig_id = str(contig["contig_id"])
            for gene in contig.get("genes", []):
                protein_id = _protein_id(accession, contig_id, int(gene["orf_index"]))
                faa_handle.write(f">{protein_id}\n")
                faa_handle.write(f"{gene['translation']}\n")

    with gzip_backend.open(gff_path, "wt", encoding="utf-8", newline="\n", compresslevel=DEFAULT_GZIP_LEVEL) as gff_handle:
        gff_handle.write("##gff-version 3\n")
        for contig in contigs:
            contig_id = str(contig["contig_id"])
            sequence = str(contig["sequence"])
            gff_handle.write(f"##sequence-region {contig_id} 1 {len(sequence)}\n")
            for gene in contig.get("genes", []):
                protein_id = _protein_id(accession, contig_id, int(gene["orf_index"]))
                start = int(gene["begin"])
                end = int(gene["end"])
                strand = "+" if int(gene["strand"]) >= 0 else "-"
                partial_begin = "true" if bool(gene.get("partial_begin")) else "false"
                partial_end = "true" if bool(gene.get("partial_end")) else "false"
                attrs = [
                    f"ID={protein_id}",
                    f"Name={protein_id}",
                    f"protein_id={protein_id}",
                    f"partial_begin={partial_begin}",
                    f"partial_end={partial_end}",
                ]
                gff_handle.write(
                    "\t".join(
                        [
                            contig_id,
                            "Pyrodigal",
                            "gene",
                            str(start),
                            str(end),
                            ".",
                            strand,
                            ".",
                            f"ID={protein_id};Name={protein_id}",
                        ]
                    )
                    + "\n"
                )
                gff_handle.write(
                    "\t".join(
                        [
                            contig_id,
                            "Pyrodigal",
                            "CDS",
                            str(start),
                            str(end),
                            ".",
                            strand,
                            "0",
                            ";".join(attrs),
                        ]
                    )
                    + "\n"
                )


def run_pyrodigal_prediction(genome_path: Path, accession: str, mode: str) -> dict[str, object]:
    """Predict coding sequences with Pyrodigal in metagenomic mode."""

    if mode != "meta":
        raise ValueError(f"unsupported P03 mode: {mode}")

    try:
        import pyrodigal
    except ModuleNotFoundError as error:  # pragma: no cover - depends on host environment
        raise RuntimeError("pyrodigal is required for P03 prediction") from error

    contigs = read_fasta_records(genome_path)
    finder = pyrodigal.GeneFinder(meta=True)
    genome_contigs: list[dict[str, object]] = []
    total_bases = 0
    coding_bases = 0
    protein_lengths: list[int] = []
    internal_stops = 0
    illegal_amino_acids = 0
    short_orfs = 0
    overlaps = 0

    for contig_id, sequence in contigs:
        total_bases += len(sequence)
        raw_prediction = finder.find_genes(sequence)
        genes = _normalize_gene_objects(raw_prediction)
        contig_genes: list[dict[str, object]] = []
        previous_end = 0

        for orf_index, gene in enumerate(genes, start=1):
            begin = _gene_int(gene, ("begin", "start"))
            end = _gene_int(gene, ("end", "stop"))
            strand = _gene_strand(gene)
            translation = _gene_translation(gene)
            partial_begin = _gene_bool(gene, ("partial_begin", "start_partial"))
            partial_end = _gene_bool(gene, ("partial_end", "end_partial"))

            coding_bases += abs(end - begin) + 1
            protein_lengths.append(len(translation))
            internal_stops += translation.count("*")
            illegal_amino_acids += sum(1 for aa in translation if aa not in ALLOWED_AA)
            short_orfs += 1 if len(translation) < 30 else 0
            if previous_end and begin <= previous_end:
                overlaps += 1
            previous_end = max(previous_end, end)

            contig_genes.append(
                {
                    "begin": begin,
                    "end": end,
                    "strand": strand,
                    "orf_index": orf_index,
                    "translation": translation,
                    "partial_begin": partial_begin,
                    "partial_end": partial_end,
                }
            )

        genome_contigs.append(
            {
                "contig_id": contig_id,
                "sequence": sequence,
                "genes": contig_genes,
            }
        )

    summary = {
        "contig_count": len(genome_contigs),
        "predicted_genes": sum(len(contig["genes"]) for contig in genome_contigs),
        "total_bases": total_bases,
        "coding_bases": coding_bases,
        "coding_density": round(coding_bases / total_bases, 6) if total_bases else 0.0,
        "mean_protein_length": round(mean(protein_lengths), 3) if protein_lengths else 0.0,
        "internal_stops": internal_stops,
        "illegal_amino_acids": illegal_amino_acids,
        "short_orfs": short_orfs,
        "overlaps": overlaps,
        "control_profiles_recovered": 0,
    }
    return {
        "accession": accession,
        "genome_path": str(genome_path),
        "contigs": genome_contigs,
        "summary": summary,
    }


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    """Read a FASTA file and return (contig_id, sequence) pairs."""

    opener = gzip.open if path.name.lower().endswith(".gz") else open
    records: list[tuple[str, str]] = []
    current_id: str | None = None
    chunks: list[str] = []

    with opener(path, "rt", encoding="utf-8") as handle:  # type: ignore[arg-type]
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(chunks).upper()))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if current_id is not None:
        records.append((current_id, "".join(chunks).upper()))
    return records


def summarize_genome_result(genome_result: dict[str, object]) -> dict[str, object]:
    """Derive a summary from a predicted genome result."""

    contigs = list(genome_result.get("contigs", []))
    summary = {
        "contig_count": len(contigs),
        "predicted_genes": 0,
        "total_bases": 0,
        "coding_bases": 0,
        "coding_density": 0.0,
        "mean_protein_length": 0.0,
        "internal_stops": 0,
        "illegal_amino_acids": 0,
        "short_orfs": 0,
        "overlaps": 0,
        "control_profiles_recovered": 0,
    }
    protein_lengths: list[int] = []
    for contig in contigs:
        sequence = str(contig["sequence"])
        genes = list(contig.get("genes", []))
        summary["total_bases"] += len(sequence)
        summary["predicted_genes"] += len(genes)
        for gene in genes:
            begin = int(gene["begin"])
            end = int(gene["end"])
            translation = str(gene["translation"])
            summary["coding_bases"] += abs(end - begin) + 1
            protein_lengths.append(len(translation))
            summary["internal_stops"] += translation.count("*")
            summary["illegal_amino_acids"] += sum(1 for aa in translation if aa not in ALLOWED_AA)
            summary["short_orfs"] += 1 if len(translation) < 30 else 0

    if summary["total_bases"]:
        summary["coding_density"] = round(summary["coding_bases"] / summary["total_bases"], 6)
    if protein_lengths:
        summary["mean_protein_length"] = round(mean(protein_lengths), 3)
    return summary


def summarize_run(genomes: list[dict[str, object]]) -> dict[str, object]:
    """Summarize a batch prediction run."""

    ok = [row for row in genomes if row.get("status") == "ok"]
    failed = [row for row in genomes if row.get("status") != "ok"]
    return {
        "ok_genomes": len(ok),
        "failed_genomes": len(failed),
        "predicted_genes_total": sum(int(row["predicted_genes"]) for row in ok),
        "mean_coding_density": round(mean(float(row["coding_density"]) for row in ok), 6) if ok else 0.0,
    }


def write_tabular_rows(path: Path, rows: Iterable[dict[str, object]], fieldnames: Sequence[str]) -> None:
    """Write a TSV file with stable row ordering."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Predict GTDB proteomes with the locked P03 route.")
    parser.add_argument("--genomes-root", required=True, type=Path, help="Root of the copied GTDB genome tree")
    parser.add_argument("--policy", required=True, type=Path, help="Prediction policy YAML")
    parser.add_argument("--out", required=True, type=Path, help="P03 output root")
    parser.add_argument("--threads", type=int, default=DEFAULT_MAX_THREADS, help="Worker threads, capped at 60")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    policy = load_prediction_policy(args.policy)
    genomes = discover_genome_files(args.genomes_root)
    result = run_p03_prediction(
        genomes,
        policy,
        args.out,
        threads=args.threads,
        genomes_root=args.genomes_root,
    )
    print(
        "P03 finished: "
        f"{result['genome_count']} genomes, "
        f"{result['summary']['predicted_genes_total']} predicted genes, "
        f"QC: {result['qc_path']}"
    )
    return 0


def _is_genome_fasta(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in KNOWN_FASTA_SUFFIXES)


def _parse_prediction_policy(policy_path: Path) -> dict[str, str]:
    policy: dict[str, str] = {}
    in_prediction = False
    with policy_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if not line.startswith(" "):
                in_prediction = stripped == "prediction:"
                continue
            if in_prediction and ":" in stripped:
                key, value = stripped.split(":", 1)
                policy[key.strip()] = value.strip().strip('"').strip("'")

    if not policy:
        raise ValueError(f"{policy_path} does not contain a prediction policy")
    return policy


def _relative_path(path: Path, root: Path) -> Path:
    if path.is_absolute() == root.is_absolute():
        return path.relative_to(root)
    return path.resolve(strict=False).relative_to(root.resolve(strict=False))


def _infer_genomes_root(paths: Sequence[Path]) -> Path:
    parents = [str(path.resolve(strict=False).parent) for path in paths]
    common = Path(os.path.commonpath(parents))
    return common


def _normalize_gene_objects(raw_prediction: object) -> list[object]:
    if hasattr(raw_prediction, "genes"):
        genes = getattr(raw_prediction, "genes")
        if callable(genes):  # pragma: no cover - defensive
            genes = genes()
        return list(genes)
    return list(raw_prediction)  # type: ignore[arg-type]


def _gene_translation(gene: object) -> str:
    for attr in ("translation", "protein_sequence", "protein", "aa_sequence"):
        value = getattr(gene, attr, None)
        if value is not None:
            return str(value)
    return ""


def _gene_int(gene: object, attributes: Sequence[str]) -> int:
    for attr in attributes:
        value = getattr(gene, attr, None)
        if value is not None:
            return int(value)
    raise AttributeError(f"gene object missing coordinate attributes: {attributes}")


def _gene_bool(gene: object, attributes: Sequence[str]) -> bool:
    for attr in attributes:
        value = getattr(gene, attr, None)
        if value is not None:
            return bool(value)
    return False


def _gene_strand(gene: object) -> int:
    strand = getattr(gene, "strand", None)
    if strand in (1, "+", b"+"):
        return 1
    if strand in (-1, "-", b"-"):
        return -1
    if strand is None:
        return 1
    return 1 if int(strand) >= 0 else -1


def _protein_id(accession: str, contig_id: str, orf_index: int) -> str:
    return f"{accession}|{contig_id}|{orf_index}"


def _qc_row_from_result(
    genome_result: dict[str, object],
    faa_path: Path,
    gff_path: Path,
    status: str,
) -> dict[str, object]:
    summary = dict(genome_result.get("summary", {}))
    return {
        "accession": genome_result.get("accession", ""),
        "genome_path": genome_result.get("genome_path", ""),
        "faa_path": str(faa_path),
        "gff_path": str(gff_path),
        "status": status,
        "error": "",
        "contig_count": summary.get("contig_count", 0),
        "predicted_genes": summary.get("predicted_genes", 0),
        "total_bases": summary.get("total_bases", 0),
        "coding_bases": summary.get("coding_bases", 0),
        "coding_density": summary.get("coding_density", 0.0),
        "mean_protein_length": summary.get("mean_protein_length", 0.0),
        "internal_stops": summary.get("internal_stops", 0),
        "illegal_amino_acids": summary.get("illegal_amino_acids", 0),
        "short_orfs": summary.get("short_orfs", 0),
        "overlaps": summary.get("overlaps", 0),
        "control_profiles_recovered": summary.get("control_profiles_recovered", 0),
    }


def _failed_qc_row(genome_path: Path, faa_path: Path, gff_path: Path, error: Exception) -> dict[str, object]:
    accession = infer_accession_from_genome_path(genome_path)
    return {
        "accession": accession,
        "genome_path": str(genome_path),
        "faa_path": str(faa_path),
        "gff_path": str(gff_path),
        "status": "failed",
        "error": str(error),
        "contig_count": 0,
        "predicted_genes": 0,
        "total_bases": 0,
        "coding_bases": 0,
        "coding_density": 0.0,
        "mean_protein_length": 0.0,
        "internal_stops": 0,
        "illegal_amino_acids": 0,
        "short_orfs": 0,
        "overlaps": 0,
        "control_profiles_recovered": 0,
    }


def _qc_fieldnames() -> list[str]:
    return [
        "accession",
        "genome_path",
        "faa_path",
        "gff_path",
        "status",
        "error",
        "contig_count",
        "predicted_genes",
        "total_bases",
        "coding_bases",
        "coding_density",
        "mean_protein_length",
        "internal_stops",
        "illegal_amino_acids",
        "short_orfs",
        "overlaps",
        "control_profiles_recovered",
    ]


def _manifest_fieldnames() -> list[str]:
    return [
        "accession",
        "genome_path",
        "faa_path",
        "gff_path",
        "status",
    ]


if __name__ == "__main__":
    raise SystemExit(main())
