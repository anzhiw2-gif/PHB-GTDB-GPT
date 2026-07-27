from __future__ import annotations

import csv
import gzip
import io
import tempfile
from pathlib import Path
import unittest

from scripts import p06_stream_proteomes as streamer


class P06StreamProteomesTest(unittest.TestCase):
    def test_streams_eligible_proteins_and_records_overlong_targets(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "proteome.faa.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write(">normal protein\nMSTN\n>overlong protein\n" + "A" * 11 + "\n")
            exclusions = root / "overlong.tsv"
            output = io.StringIO()

            summary = streamer.stream_proteomes([source], max_protein_length=10, exclusion_path=exclusions, output=output)

            self.assertEqual(summary, {"emitted_sequences": 1, "excluded_overlong_sequences": 1})
            self.assertEqual(output.getvalue(), ">normal protein\nMSTN\n")
            with exclusions.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows, [{"source_proteome_path": source.as_posix(), "target_id": "overlong", "sequence_length_aa": "11", "reason": "hmmsearch_target_length_gt_100000"}])
