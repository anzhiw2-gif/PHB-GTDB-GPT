from __future__ import annotations

import csv
from pathlib import Path
import unittest


ROOT = Path(__file__).parents[1]
MODEL_REGISTRY = ROOT / "04_family_profiles" / "manifests" / "p05_hmm_model_registry.tsv"
CORE_SEEDS = ROOT / "04_family_profiles" / "manifests" / "p05_extracellular_core_seed_registry.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class P05ExtracellularProvenanceTest(unittest.TestCase):
    def test_core_seed_registry_is_locked_to_approved_core_model(self) -> None:
        models = {row["family_category"]: row for row in read_tsv(MODEL_REGISTRY)}
        core = models["extracellular_pha_depolymerase_core"]
        core_seeds = read_tsv(CORE_SEEDS)

        self.assertEqual(core["approved_for_p06"], "yes")
        self.assertEqual(core["scan_permission"], "approved")
        self.assertEqual(len(core_seeds), 17)
        self.assertEqual(
            {row["model_sha256"] for row in core_seeds},
            {core["model_sha256"]},
        )
        self.assertEqual(len({row["source_accession"] for row in core_seeds}), 17)

    def test_archived_subtypes_remain_blocked_from_direct_p06_scanning(self) -> None:
        models = {row["family_category"]: row for row in read_tsv(MODEL_REGISTRY)}
        subtypes = (
            "extracellular_mcl_pha_dep",
            "extracellular_scl_pha_dep_type_I",
            "extracellular_scl_pha_dep_type_II",
        )
        for family in subtypes:
            self.assertEqual(models[family]["approved_for_p06"], "no")
            self.assertEqual(models[family]["scan_permission"], "blocked")
            self.assertEqual(models[family]["model_status"], "superseded_by_extracellular_core")


if __name__ == "__main__":
    unittest.main()
