# P05 HMM Seed Selection Decision

**Decision date:** 2026-07-27
**Status:** user-approved on 2026-07-27; the three affected models were rebuilt
and checksum-locked. This historical gate was satisfied after all six profiles
received calibration decisions: exactly four models were approved for P06, and
their scan/parse acceptance is recorded in `P06_COMPLETION_2026-07-28.md`.

## Locked Screening Rules

These rules are part of the HMM construction contract.

1. **Bacterial profile seeds require experimental support.** Acceptable support
   includes direct biochemical activity, purification/assay of the exact gene
   product, or a gene knockout/complementation or physiological experiment
   that maps unambiguously to the accession. Annotation-only E3 rows are not
   profile seeds.
2. **Archaeal coverage seeds do not all require direct enzyme experiments.**
   Accessioned E3 rows may be used when experimental coverage is sparse, but
   the organism, accession, architecture evidence, and exact support scope
   must be recorded. E3 sequence evidence must not be described as phenotype
   evidence.
3. A family HMM must represent a coherent sequence architecture. A product
   name such as "PHB depolymerase" is not sufficient when domain annotations,
   catalytic motifs, or sequence alignment place the record in another branch.
4. A family does not need artificial taxonomic balance. Taxonomic enrichment
   is retained when it is supported by the biology and the literature; it is
   not corrected by adding weak or unrelated seeds.

## `intracellular_mcl_pha_dep`

The Pseudomonas enrichment is treated as an intrinsic family characteristic,
not as a seed-selection artifact. De Eugenio et al. described the known
intracellular mcl-PHA mobilization systems as occurring in *Pseudomonas*
species and biochemically characterized the *P. putida* KT2442 PhaZ as a
granule-localized intracellular mcl-PHA depolymerase. The family is therefore
allowed to remain Pseudomonas-enriched. The result label must state that this
is family coverage and sequence evidence, not proof that every hit has a PHB
or PHA degradation phenotype.

Primary anchor: *Pseudomonas putida* KT2442, UniProtKB `Q5Y152`, PMID
`17170116`, DOI `10.1074/jbc.M608119200`.

The user-approved profile bundle contains `Q5Y152`, `B7UCC9`, `Q88D24`, and
`Q9R9W3`. The remaining audited records `Q5Q135`, `Q8VV57`, `Q9AGB5`, and
`Q9Z3Y0` remain in the reference library as `boundary_candidate` controls;
they do not enter the next HMM alignment. This preserves the genuine
*Pseudomonas* enrichment without using weak non-*Pseudomonas* sequences to
manufacture taxonomic balance.

## `archaeal_patatin_like_pha_dep`

### Direct anchor

| Accession | Organism | Evidence | Sequence evidence |
|---|---|---|---|
| `AFK21580.1` | *Haloferax mediterranei* ATCC 33500 | E1; PMID `25710370`; DOI `10.1128/AEM.04269-14` | PhaZh1/HFX_6464; patatin `PF01734`, `IPR002641`, `IPR016035`; `GTSGG` at residue 45 |

`AFK21580.1` is the archival NCBI Protein accession. UniProtKB `I3RBH0`
is a secondary cross-check and maps to KEGG `HFX_6464` and RefSeq
`WP_004060664.1`.

### Proposed coherent coverage set

The following accessioned candidates have the patatin-like RssA annotation,
the `IPR002641/IPR016035` cross-references where available, and a conserved
`GxSxG` motif near the N terminus. They are approved E3 coverage candidates;
their NCBI FASTA retrieval and SHA256 verification are recorded in the P05
retrieval log and seed provenance:

| Accession | Organism | Length aa | Motif | Role |
|---|---|---:|---|---|
| `CCQ36014.1` | *Natronomonas moolapensis* | 323 | `GTSGG` | retain as current E3 cross-genus coverage |
| `CCQ32286.1` | *Halorhabdus tiamatea* | 322 | `GSSGG` | add as independent E3 coverage |
| `AGN01047.1` | *Salinarchaeum* sp. Harcht-Bsk1 | 321 | `GTSGG` | add as independent E3 coverage |
| `KYH27761.1` | *Halalkalicoccus paucihalophilus* | 329 | `GTSGG` | add as independent E3 coverage |
| `ELY43313.1` | *Natronorubrum tibetense* | 269 | `GTSGG` | add as shorter E3 architecture coverage |

The observed lengths and motif positions are compatible with the 321-aa
PhaZh1-like core and are more coherent than the current long records.

### Demote from this model

The following current rows remain useful for boundary review, but should not
be profile seeds for the patatin-like HMM without new evidence:

`AHB64615.1`, `AHZ23723.1`, `AJF25805.1`, `EFW93255.1`, `KOX95185.1`, and
`KZX50211.1` are 442--490 aa records whose NCBI annotations include
PHB-synthase-like, HHH, PKD, or AxeA-related architectures and which lack the
PhaZh1-like N-terminal `GxSxG` motif. Their automatic product name is not
enough to merge them into the patatin-like branch. Keep them as explicitly
labelled non-patatin archaeal boundary candidates or a future separate family.

## `intracellular_phaZ_no_lipase_box`

### Proposed profile set

| Accession | Organism | Evidence | Role |
|---|---|---|---|
| `O87189` | *Cupriavidus necator* H16 PhaZ1/PhaZre | E2; PMID `11114905`; DOI `10.1128/JB.183.1.94-100.2001` | direct intracellular PHB assay and granule localization |
| `Q0K7T2` | *Cupriavidus necator* H16 PhaZ2 | E1; PMID `12775684`; DOI `10.1128/JB.185.12.3485-3490.2003` | purified enzyme and PHB-oligomer/PHB activity |
| `Q71KW6` | *Azospirillum brasilense* PhaZ | E2; PMID `12898135`; DOI `10.1007/s00203-003-0590-z` | cloned gene and phaZ mutant unable to degrade PHB |
| `Q92TD3` | *Sinorhizobium meliloti* PhaZ | E2; PMID `20346169`; DOI `10.1186/1471-2180-10-92` | cloned gene, mutant phenotype, and complementation |
| `Q9WX79` | *Paracoccus denitrificans* PhaZ | E2; PMID `11267773`; DOI `10.1111/j.1574-6968.2001.tb10558.x` | functional expression and PHB-granule degradation assay |

This set represents four genera while retaining the two experimentally
distinct H16 core paralogs. `Q9WX79` is a better cross-genus profile seed
than an uncharacterized H16 paralog because it has direct experimental
support and maps to the same `PF06850/TIGR01849` family space.

### Boundary control

`Q0K4D5` is current *Cupriavidus* PhaZ3. The 2003 deletion study identified
and tested its in vivo role, but concluded that its depolymerase role remained
to be established; later work associates it more strongly with granule
architecture than with primary PHB mobilization. It is now a labelled
`boundary_candidate` control, not an equal-weight profile seed.

## Model and Scan Gate

The existing six HMM files are provisional artifacts. No new P06 scan may use
them until the revised seed manifest is approved, the profiles are rebuilt or
explicitly retained, and the exact HMM paths plus SHA256 values are written to
the scan manifest. HMMER hits remain sequence evidence only.

## Sources

- De Eugenio et al. 2007, DOI `10.1074/jbc.M608119200`, PMID `17170116`:
  https://digital.csic.es/handle/10261/160503
- Liu et al. 2015, DOI `10.1128/AEM.04269-14`, PMID `25710370`:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4393451/
- Saegusa et al. 2001, DOI `10.1128/JB.183.1.94-100.2001`:
  https://journals.asm.org/doi/10.1128/jb.183.1.94-100.2001
- Kobayashi et al. 2003, DOI `10.1128/JB.185.12.3485-3490.2003`:
  https://journals.asm.org/doi/10.1128/jb.185.12.3485-3490.2003
- York et al. 2003, DOI `10.1128/JB.185.13.3788-3794.2003`, PMCID `PMC161563`:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC161563/
- Trainer et al. 2010, DOI `10.1186/1471-2180-10-92`, PMCID `PMC2867953`:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2867953/
- Gao et al. 2001, DOI `10.1111/j.1574-6968.2001.tb10558.x`, PMID `11267773`:
  https://doi.org/10.1111/j.1574-6968.2001.tb10558.x
- NCBI Protein E-utilities:
  https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
