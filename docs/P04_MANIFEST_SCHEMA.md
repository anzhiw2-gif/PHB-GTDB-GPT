# P04 Manifest Schema

P04 stores seed-sequence provenance in a normalized TSV manifest.

## File

`01_reference_library/reference_library.normalized.tsv`

## Required Columns

- `seed_id`
- `reference_library`
- `taxonomic_domain`
- `family_category`
- `seed_name`
- `evidence_level`
- `source_database`
- `source_accession`
- `organism`
- `taxon_id`
- `retrieval_date`
- `sequence_format`
- `sequence_length_aa`
- `sequence_path`

## Optional Columns

- `profile_seed_status`
- `family_label`
- `source_release`
- `source_version`
- `source_url`
- `retrieval_method`
- `retrieval_query`
- `retrieval_endpoint`
- `retrieval_batch_id`
- `retrieval_log_path`
- `accession_version`
- `taxon`
- `database_version`
- `doi`
- `pmid`
- `pmcid`
- `reference_title`
- `reference_year`
- `exclusion_reason`
- `record_kind`
- `literature_support_scope`
- `supporting_sources`
- `supporting_accessions`
- `supporting_notes`
- `notes`

## Field Rules

### `seed_id`

Stable project-local identifier. Keep it short, unique, and reproducible.

### `reference_library`

Must be one of:

- `bacteria_high_confidence`
- `archaea_literature_supported`

This field controls which normalized split library receives the row.

### `taxonomic_domain`

Must be `Bacteria` or `Archaea`.

Bacterial and archaeal references are built as separate libraries. Do not mix
archaeal candidate rows into the bacterial high-confidence seed set.

### `family_category`

Must be one of the canonical P04 categories listed in
[docs/P04_REFERENCE_LIBRARY.md](P04_REFERENCE_LIBRARY.md).

### `evidence_level`

Must be one of `E1`, `E2`, `E3`, or `Excluded`.

### `source_database`

Primary database for the seed sequence. Use one of:

- `UniProtKB`
- `NCBI Protein`
- `NCBI Nucleotide`
- `PDB`
- `Legacy verified`

### `source_accession`

Primary accession or accession.version from the source database. Keep the
version suffix when the source provides one.

### `organism`

Free-text organism name copied from the source record or verified literature.

### `taxon_id`

Numeric NCBI taxon ID. Do not leave it blank for a seed row.

### `retrieval_date`

ISO date in `YYYY-MM-DD` format.

### `sequence_format`

Must be `fasta` or `faa` for seed sequences.

### `sequence_length_aa`

Positive integer amino-acid length of the stored seed sequence.

### `sequence_path`

Repository-relative path to the stored seed FASTA file.

### `profile_seed_status`

Optional P05 HMM-admission state. Leave blank only for legacy rows that have
not yet received an accession-level review. Accepted non-blank values are:

- `approved`: the row may contribute to a P05 seed bundle when its evidence
  level also satisfies the bacterial or archaeal library rule.
- `boundary_candidate`: retain the accession and provenance for architecture
  review or calibration controls, but exclude it from every P05 HMM seed
  bundle.
- `pending_evidence_audit`: retain the row while experimental support is
  checked; exclude it from every P05 HMM seed bundle until it is approved.

For bacterial rows, `approved` additionally requires experimental support
tied to the accession. For archaeal E3 coverage rows, it records only
architecture-supported coverage and must not be interpreted as phenotype
evidence.

### Conditional fields

- `exclusion_reason` is required when `evidence_level` is `Excluded`
- `reference_library` must be `bacteria_high_confidence` when
  `taxonomic_domain` is `Bacteria`
- `reference_library` must be `archaea_literature_supported` when
  `taxonomic_domain` is `Archaea`
- non-excluded `bacteria_high_confidence` rows must use `E1` or `E2`; `E3`
  annotation-only rows are not allowed in the bacterial seed library
- non-excluded `archaea_literature_supported` rows may use `E1`, `E2`, or
  `E3`; `E1`/`E2` rows must include at least one of `pmid`, `doi`, or `pmcid`,
  and all archaeal rows must describe the support scope in
  `literature_support_scope`
- `source_url` should be filled when the source provides a stable landing page
- `doi`, `pmid`, and `pmcid` should be filled whenever a seed is anchored to a
  paper
- `retrieval_query`, `retrieval_endpoint`, and `retrieval_method` should be
  filled when the row was produced from a database search instead of a single
  accession fetch
- `retrieval_batch_id` and `retrieval_log_path` should point to the log record
  that reconstructed the row

## Normalization Rules

- Sort rows by `family_category`, evidence level, and `seed_id`
- The normalizer writes a combined manifest plus separate bacterial and
  archaeal normalized manifests
- Keep empty optional fields as empty strings
- Do not duplicate the same seed accession in multiple rows unless the rows
  differ by evidence status or source role and that distinction is documented
- Keep one primary source per row; secondary support belongs in the supporting
  provenance fields or the notes column

## Validation Rules

- reject unknown family categories
- reject unknown reference libraries
- reject unknown taxonomic domains
- reject mismatches between `taxonomic_domain` and `reference_library`
- reject unknown evidence levels
- reject unknown source databases
- reject non-ISO retrieval dates
- reject non-numeric taxon IDs
- reject non-positive amino-acid lengths
- reject unknown sequence formats
- reject unknown non-blank `profile_seed_status` values
- reject `Excluded` rows that do not explain why they were excluded
- reject bacterial high-confidence rows with annotation-only `E3` evidence
- reject archaeal rows that omit `literature_support_scope`; reject `E1`/`E2`
  archaeal rows without a literature identifier
