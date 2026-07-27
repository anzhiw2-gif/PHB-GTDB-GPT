# P04 Retrieval Contract

P04 uses a bounded, auditable retrieval process to collect reference sequences
and literature evidence for the seed library.

## Scope

The contract covers:

- literature-backed seed discovery
- sequence retrieval from primary databases
- accession and taxon reconciliation
- retrieval logging for later audit

The contract does not mix family-definition anchors with seed sequences.
DED and TIGRFAM family-definition rows are tracked separately from seed rows.

P04 also keeps bacterial and archaeal seed libraries separate:

- bacterial rows belong to `bacteria_high_confidence`
- archaeal rows belong to `archaea_literature_supported`

This separation prevents sparse archaeal evidence or automatic annotations from
changing the bacterial high-confidence seed set. Archaeal rows may be
literature-backed or annotation-supported, but they stay in the archaeal
library and never bleed into the bacterial seed set.

## Primary Source Rules

Use one primary source per seed row.

Priority order:

1. `UniProtKB` for curated protein sequence records and names
2. `NCBI Protein` for archival RefSeq or GenBank protein records
3. `NCBI Nucleotide` when the authoritative record is the coding nucleotide
   sequence rather than the translated protein
4. `PDB` only when the seed is explicitly structure-derived and the sequence is
   backed by the entry
5. `Legacy verified` only for historical records that were independently checked
   against the original source

Literature identifiers such as `PMID`, `DOI`, and `PMCID` are supporting
provenance, not replacements for a sequence accession.

## Taxonomic Library Rules

### Bacterial high-confidence library

Use this library for bacterial PHB/PHA depolymerase seeds with `E1` or `E2`
evidence only. `E1` requires direct experimental activity. `E2` can include
strong functional evidence with explicit limitations, such as curated protein
records backed by literature or characterized gene products where the sequence
identity is clear.

Keep the PHB-specific nucleus strict. Admit a small number of literature-
supported PHA bridge rows only when they are needed to keep the major DED
family branches represented; mark them clearly in notes so they do not blur the
core PHB set.

Do not admit annotation-only `E3` bacterial rows as seeds. Store them as
tracked candidates or excluded rows when they are useful for boundary review.

For HMM construction, "experimental support" means that the exact bacterial
accession is tied to direct biochemical activity, purification/assay of the
gene product, or an unambiguous knockout/complementation or physiological
experiment. A cloned locus, product name, or homology annotation alone is not
enough for a bacterial profile seed.

### Archaeal literature-supported library

Use this library for archaeal PHB/PHA depolymerase-related seeds when the
sequence is traceable to a source accession and the support scope is explicit.
Direct literature support is preferred for the anchor rows, but annotation-only
archaeal E3 records may be admitted when they are the best available evidence
for family coverage. Every non-excluded archaeal row must fill
`literature_support_scope` with the exact support type.

Acceptable archaeal support scopes include:

- direct depolymerase assay on native PHA, nPHB, or nPHBV
- granule association plus mutational or biochemical evidence
- literature-curated gene or protein association with clear accession mapping

Automatic names such as "PHB depolymerase family esterase", InterPro/Pfam
membership, or NCBI/UniProt TrEMBL annotation alone are not sufficient for an
E1/E2 archaeal anchor row. They can support an E3 archaeal seed row only when
the accession, organism, and support scope are explicitly recorded.

The first archaeal literature anchor is *Haloferax mediterranei* `PhaZh1` /
`HFX_6464`, supported by PMID `25710370` and DOI `10.1128/AEM.04269-14`.
Use NCBI Protein accession `AFK21580.1` as the archival protein accession and
UniProtKB `I3RBH0` as a secondary cross-check. Its patatin-like architecture
(`PF01734`, `IPR002641`, `IPR016035`) should not be merged with the bacterial
`PF10503`/`TIGR01840`, `TIGR01849`, or `TIGR02240` families.

Before rebuilding the archaeal patatin-like model, require a coherent
patatin-like architecture and conserved N-terminal `GxSxG` evidence for the
coverage rows. Automatic archaeal records with PHB-synthase-like, HHH, PKD,
or AxeA-related architecture are boundary candidates rather than seeds for
this model, even when their product name says "PHB depolymerase".

The `intracellular_mcl_pha_dep` family is expected to be Pseudomonas-enriched
because intracellular mcl-PHA mobilization has been experimentally described
primarily in *Pseudomonas*. This is an intrinsic family property, not a reason
to force artificial cross-genus seed balance. The result must still be labelled
as sequence/family evidence, not phenotype proof.

## Retrieval Interfaces

### PubMed literature lookup

Use NCBI E-utilities for literature discovery and citation verification.

- official help: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- quick start: https://www.ncbi.nlm.nih.gov/books/NBK25500/

Use `esearch` to find candidate papers, then `esummary` or `efetch` to record
the publication metadata that supports the seed decision. Record the exact query
terms, publication identifiers, and access date.

### NCBI Protein retrieval

Use NCBI E-utilities against `db=protein` for accessioned protein records.

- official help: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- base URL: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/

Preserve accession versions where they exist. Prefer `esearch` for lookup,
`esummary` for metadata, and `efetch` for the actual record or FASTA sequence.

### UniProtKB retrieval

Use the UniProt REST API for curated protein records.

- programmatic access: https://www.uniprot.org/help/api_queries
- query fields: https://www.uniprot.org/help/query-fields
- API base URL: https://rest.uniprot.org

Prefer accession-based retrieval. For search-based discovery, record the query,
field filters, page size, cursor, and returned accession list. If the row is
built from a UniProtKB record, preserve the primary accession and the reviewed
status.

### InterPro and family-anchor retrieval

Use InterPro only for family and domain support. InterPro, Pfam, and NCBIfam
membership can support a family label or boundary check, but it is not by
itself a seed-level phenotype claim.

## Retrieval Log Rules

Every retrieval batch should record:

- `retrieval_id`
- `source_database`
- `retrieval_endpoint`
- `retrieval_query`
- `retrieval_date`
- `database_release` or `record_version`
- `expected_count`
- `retrieved_count`
- `page_size`
- `page_cursor` or `page_offset`
- `local_filters`
- `supporting_accessions`
- `taxonomic_domain`
- `reference_library`
- `literature_support_scope`
- `notes`

For exhaustive retrievals, count first when the API exposes a count. Paginate
in a stable order and stop visibly if counts disagree, pagination ends early,
or the retrieval would exceed 10,000 records or 100 API calls without explicit
confirmation.

## Retrieval Boundaries

- Do not infer evidence from search-engine snippets.
- Do not use unverified secondary summaries as the only proof of a seed.
- Keep literature evidence separate from sequence evidence.
- Keep family-definition anchors separate from seed rows.
- Keep bacterial and archaeal seed libraries separate.
- Do not promote archaeal automatic annotation-only records to seed status.
- Record any accession mapping explicitly when a source accession must be
  converted before retrieval.
