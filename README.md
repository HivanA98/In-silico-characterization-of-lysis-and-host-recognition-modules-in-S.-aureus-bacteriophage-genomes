# Supplementary Code and Data — *In silico* characterization of lysis and host-recognition modules in *Staphylococcus aureus* bacteriophage genomes

**Phage Characterization Toolkit v3.2 — *Staphylococcus aureus* release**

**Associated manuscript:**
> Hasugian IA, Alifiyah NI. *In silico* characterization of lysis and host-recognition modules in *Staphylococcus aureus* bacteriophage genomes. *Memórias do Instituto Oswaldo Cruz* (under review; MIOC-2026-0245).

This archive contains the code, host and phage genome records, profile, alignments and
intermediate results required to reproduce every table and figure in the manuscript.

---

## 1. Scope of this release

The toolkit is built as **one shared engine (`phagecore`) plus a host profile**: the
analysis code is host-neutral, and everything specific to a bacterial host lives in a
YAML profile that supplies plausibility bounds, annotation vocabulary and curated
evidence. Every stage that depends on host biology **requires** `--profile` and exits
with an error if it is not given, so no host is ever a silent default.

**This deposit contains the *Staphylococcus aureus* profile only.** Profiles for other
bacterial hosts exist in the development version and are not part of this release;
they are unrelated to the manuscript and are withheld pending separate publication.
Stages S8 and S9 of the development pipeline are likewise outside the scope of the
manuscript and are not deposited.

| Stage | Script | Manuscript output |
|-------|--------|-------------------|
| S1 | `S1_genome_statistics.py` | **Table I** — size, GC, CDS, tRNA, taxonomy, QC |
| S2 | `S2_holin_rbp_annotation.py` | **Table II** (holin, RBP columns), **Table III** |
| S3 | `S3_terl_extractor.py` | **Figure 1** — TerL multi-FASTA for alignment |
| S4 | `S4_endolysin_extractor.py` | **Table II** (endolysin product, length, domains) |
| S4b | `S4b_endolysin_domain_identity.py` | **Supplementary Table S4** — per-domain pairwise identity |
| S5 | `S5_rbp_extractor.py` | not used in the manuscript; included for completeness |
| S6 | `S6_tRNA_analyzer.py` | de novo tRNA counts (Results) |
| S7 | `S7_codon_trna_coverage.py` | host RSCU, phage tRNA codon coverage |
| S7b | `S7b_phage_codon_usage.py` | **Supplementary Table S2** — phage-vs-host codon usage |

`S4_endolysin_extractor_for_interpro.py` is the script used for the **original
submission** and is retained so that the first version of Table II can be reproduced.
`S4_endolysin_extractor.py` is the current v3.2 script and is the one used for the
revision. They are not interchangeable.

---

## 2. Requirements

| Dependency | Version used | Notes |
|---|---|---|
| Python | 3.12.10 | ≥ 3.10 required |
| biopython | 1.87 | all stages |
| pandas | 3.0.3 | CSV output |
| PyYAML | any | loading the profile |
| BLAST+ (`tblastn`) | 2.16.0 | optional, `--run-tblastn` |
| MAFFT | web server, run 16 August 2026 | https://mafft.cbrc.jp/alignment/server/ |
| MEGA | 12.1.2 (build 12251216-x86_64) | https://www.megasoftware.net/ |
| InterPro | 108.0 | https://www.ebi.ac.uk/interpro/ |
| tRNAscan-SE | 2.0 | de novo tRNA prediction |
| VIRIDIC | web server | intergenomic similarity |
| Pharokka | 1.10.1 | re-annotation arm |
| Phold | 1.3.0 | re-annotation arm |
| Phynteny_transformer | 0.1.3 | residual hypothetical proteins |

```
pip install biopython==1.87 pandas==3.0.3 pyyaml
```

Re-annotation with Pharokka, Phold and Phynteny was run in Docker; PyTorch 2.8.0+cu128
on an NVIDIA RTX 4060, with `setuptools<82` pinned (Phynteny 0.1.3 imports
`pkg_resources`, which setuptools removed in 82.0.0). All other stages run on Windows
with the Python packages above and require no compilation.

---

## 3. Contents

```
GenBank/                          22 deposited GenBank records (one genome per file)
GenBank/Pharokka_Phold/           the same 22 genomes re-annotated (Pharokka -> Phold)
GenBank/Phynteny/phynteny.gbk     Phynteny output: ALL 22 genomes in ONE file
host/Staphylococcus_aureus/       host genomes for S7 codon reference
phage_characterization/profiles/  AureusPhage.yaml — the host profile
phagecore/                        shared engine
S1..S7, S4b, S7b                  stage scripts
Mega Phylogeny Result/            Newick export, MEGA session, Figure 1 SVG
alignments/                       TerL and endolysin-domain alignments
interpro/                         InterPro TSV results
results/                          per-stage CSV outputs for the three annotation arms
viridic/                          intergenomic similarity matrix and heatmap
```

**`GenBank/Phynteny/phynteny.gbk` holds all 22 genomes in a single file.** The stage
scripts take one genome per file (`phagecore/genbank_io.py` reads the first record of
each file), so this file must be split before it can be used as input. It is deposited
in the form Phynteny writes it.

---

## 4. Reproducing the manuscript

Input convention: phage genomes in `GenBank/`, host genomes in
`host/Staphylococcus_aureus/`, outputs in `results/`.

### Table I

```
python S1_genome_statistics.py -i GenBank -o results\result_aureus_04 ^
       --profile phage_characterization\profiles\AureusPhage.yaml
```

Class, family and subfamily are read from `record.annotations["taxonomy"]` by ICTV
rank suffix (`-viricetes`, `-viridae`, `-virinae`). Subfamilies that NCBI places in no
family — Azeredovirinae, i.e. EW (NC_007056) and SA13 (NC_021863) — are recorded as
family *Unassigned* and flagged `family_unresolved_verify_ICTV` rather than left blank.

### Table II — holin and RBP columns

```
python S2_holin_rbp_annotation.py -i GenBank -o results\result_aureus_04 ^
       --profile phage_characterization\profiles\AureusPhage.yaml
```

Holin and RBP presence are called from the wording of the CDS `product` qualifier
using the keyword sets in the profile (`holin_keywords`, 5 terms; `rbp_keywords`,
12 terms). Every positive call is written with the matched string in the
`Holin_Evidence` / `RBP_Evidence` columns so it can be checked by hand. The structural
tail-tube protein is deliberately excluded from the RBP vocabulary, since it is present
in essentially every tailed phage and would make the column uninformative.

### Table II — endolysin columns

```
python S4_endolysin_extractor.py -i GenBank -o results\result_aureus_04 ^
       --profile phage_characterization\profiles\AureusPhage.yaml --run-tblastn
```

The script collects **every** lysis-keyword CDS and ranks them, rather than taking the
first keyword match; the ranked audit is written to `endolysin_audit_*.csv` with the
evidence for each row. Unique candidate sequences are written to
`endolysin_unique_*.faa` for InterPro submission.

`--run-tblastn` is applied **conditionally, not to all 22 genomes**: it runs only where
the keyword scan produced a weak pick (inferred domain containing `nlpc`, or lacking
CHAP/amidase/LysK evidence) or produced no free endolysin at all. Command and
parameters:

```
tblastn -query <MN336261 Sb1_8383 LysK, 495 aa> -subject <genome.fna>
        -outfmt "6 qseqid sseqid pident length sstart send evalue bitscore"
        -max_target_seqs 5
```

All other parameters are BLAST+ defaults (e-value 10, word size 3, BLOSUM62, gap open
11, gap extend 1, SEG filtering on). Acceptance threshold `--identity-threshold`,
default **90.0 %**: at or above it the hit is recorded as `tblastn-recovered`, below it
as `divergent-uncertain`. Three genomes met the trigger — Maine (MN045228), JD007
(NC_019726) and Twort (NC_007021).

### Supplementary Table S4 — per-domain identity

```
python S4b_endolysin_domain_identity.py ^
       --faa results\result_aureus_04\endolysin_unique_AureusPhage_confirmed.faa ^
       --interpro results\result_aureus_04\interpro_endolysin\endolysin_unique_AureusPhage.tsv ^
       -o results\result_aureus_04\endolysin_identity ^
       --profile phage_characterization\profiles\AureusPhage.yaml ^
       --status free-endolysin
```

Domain boundaries are **read from the InterPro output already used for Table II**, not
re-predicted, so the domain definition cannot drift between the two tables. Where more
than one member database annotates the same domain the widest span is taken and its
source recorded: CHAP `PS50911`, amidase `G3DSA:3.40.80.10`, SH3b `G3DSA:2.30.30.40`.
The Pfam amidase signature `PF01510` was returned for too few sequences to define the
set. Each partition is aligned separately with MAFFT L-INS-i and percent identity is
computed over columns in which **both** sequences carry a residue; columns with a gap
in either sequence are excluded from numerator and denominator alike.

If `mafft` is not on `PATH` the script stops after writing the per-domain FASTA files
and prints the submission instructions; re-run with `--align read` after saving each
web-server result as `<domain>.aln.faa`.

### Figure 1 — TerL phylogeny

```
python S3_terl_extractor.py -i GenBank -o results\result_aureus_04\TerL_combined.faa
```

S3 takes no profile: the terminase large subunit is universal.

Two detection mechanisms are applied to each CDS `product` qualifier:

1. **keyword match** — any of 8 substrings (`terminase large subunit`,
   `large terminase`, `terl`, …) in the lower-cased product;
2. **exact product match** — the product is exactly `Ter` or `ter`. Seven Kayvirus
   genomes (EU418428, NC_047722–NC_047727) annotate their 605-aa TerL with this
   three-letter abbreviation and are missed by keyword search alone.

**20 of 22 genomes yielded a TerL.** Portland (MT926124) and vB_SauP-436A1 (MN150710)
carry no CDS annotated as a terminase; the nearest candidates are a 415-aa
*"putative encapsidation protein"* and a 415-aa *"DNA packaging protein"* respectively.
Both are micro-class genomes below 20 kb. They are excluded from the phylogeny, and
this is stated in the manuscript Methods.

Alignment and tree:

```
MAFFT web server, 16 August 2026
  strategy : L-INS-i   (all other settings left at default)
  output   : TerL_aligned.faa  (20 sequences, 519 aligned sites)

MEGA 12.1.2 -> Phylogeny -> Construct/Test Maximum Likelihood Tree
  Statistical method        : Maximum Likelihood
  Test of phylogeny         : Standard Bootstrap, 1000 replicates
  Substitution model        : LG
  Rates among sites         : Gamma distributed with invariant sites (G+I)
  Discrete gamma categories : 5
  Gaps/missing data         : Partial deletion
  Site coverage cutoff      : 80 %
  ML heuristic method       : Nearest-Neighbour-Interchange (NNI)
  Initial tree for ML       : automatic (NJ/MP)
  Branch swap filter        : None
  Threads                   : 12
```

Run statistics from the MEGA session (`Mega Phylogeny Result/`): 20 taxa, 519 sites,
516 common sites, 118 invariant sites, lnL −3926.043, BIC 8210.524, AICc 7930.405,
sum of branch lengths 4.662.

### Supplementary Table S2 — codon usage

```
python S6_tRNA_analyzer.py  -i GenBank -o results\result_aureus_04 ^
       --profile phage_characterization\profiles\AureusPhage.yaml

python S7_codon_trna_coverage.py --host host\Staphylococcus_aureus ^
       --trna results\result_aureus_04\tRNA_detailed_Aureus.csv ^
       -o results\result_aureus_04\codon_analysis ^
       --profile phage_characterization\profiles\AureusPhage.yaml

python S7b_phage_codon_usage.py --phage-dir GenBank ^
       --host host\Staphylococcus_aureus ^
       -o results\result_aureus_04\codon_analysis ^
       --profile phage_characterization\profiles\AureusPhage.yaml
```

S7 computes host RSCU and maps each de novo tRNA anticodon to the codons it decodes.
S7b adds the missing half — per-genome phage RSCU and its similarity to the host —
reusing `count_codons_from_cds()` and `calculate_rscu()` from `phagecore.codon` so the
RSCU definition cannot differ between host and phage. Similarity is computed over the
59 codons of degenerate synonymous families; methionine and tryptophan are excluded,
being invariant at RSCU = 1.

**Host set.** Codon usage is a species-level property, so counts are summed across the
files in `host/Staphylococcus_aureus/`. That directory contains six files representing
**five** distinct strains: `GCA_002310435.gbff` (CP023390.1) and `GCF_002310435.gbff`
(NZ_CP023390.1) are the GenBank and RefSeq copies of the same 2,878,897-bp assembly,
so that strain contributes twice. The redundancy was present in the run reported in the
manuscript and is left in place here so the deposit reproduces it exactly. Its effect
is negligible and was measured: removing the duplicate changes host GC3 from 22.73 % to
22.72 %, gives an RSCU correlation of r = 0.999999 between the two host tables, and a
maximum per-codon |ΔRSCU| of 0.0042.

### Intergenomic similarity

VIRIDIC web server, default BLASTN parameter set
`-word_size 7 -reward 2 -penalty -3 -gapopen 5 -gapextend 2`, species demarcation 95 %
and genus demarcation 70 % as recommended by the ICTV Bacterial and Archaeal Viruses
Subcommittee. Input: the 22 retained genomes plus the two excluded during curation
(vB_SauP_EBHT NC_055906 and MarsHill MW248466), 24 sequences in total.

### Re-annotation arms (Table III)

```
Pharokka 1.10.1  ->  Phold 1.3.0  ->  Phynteny_transformer 0.1.3
```

Phold was run on GPU at half precision (ProstT5 encoder + Foldseek structural search);
Phynteny assigns a PHROG **category** to proteins that remain of unknown function, not
a product name, and was therefore reported separately rather than used to make module
calls. The three detection arms differ **only** in the source of the CDS product
names: keyword sets, InterPro protocol and curation funnel are identical.

---

## 5. The profile

`phage_characterization/profiles/AureusPhage.yaml` carries every host-specific value
used by the analysis, and it is deposited in full so that each call is auditable:

- QC plausibility envelope — genome size 15–300 kb, GC 24–38 %, minimum genome size
  for the zero-CDS check. **These bounds were derived empirically from a wider
  calibration set of *S. aureus* phage genomes, not from the 22 analysed here**; the
  observed range in this dataset is 17.5–148.6 kb and 29.34–35.99 % GC.
- Annotation vocabulary — `holin_keywords` (5), `rbp_keywords` (12), `lysis_keywords`
  (16 terms, a deliberately broad collection net), `vaph_markers` (9),
  `non_endolysin_markers` (18, a precision denylist of broad-keyword false positives
  confirmed non-endolysin by InterPro).
- Endolysin windows — canonical 440–520 aa, plausible 200–600 aa, VAPH length
  threshold 700 aa.
- `known_cases` — eight manually validated exceptions with status, identity and a
  literature note. These are curated evidence, not hardcoded shortcuts: each records a
  case where the deposited annotation alone gives the wrong answer (intron-split
  lysK.1, HNH-disrupted ORFs, tBLASTn-recovered endolysins, a tail-anchored VAPH, and
  one divergent endolysin).
- S6/S7 anchors — tRNA ground-truth accession, canonical isotypes, host codon
  reference NC_007795.1.

The profile also contains fields consumed only by stage S5, which is included in this
archive for completeness but produces no manuscript output.

---

## 6. Known limitations of this release

- **Annotation dependence is the principal limitation.** Stages S1–S4 extract and
  curate what the deposited annotation contains; they are not independent gene callers.
  The magnitude of that dependence is quantified in Table III of the manuscript.
- **Detection vocabulary is a second dependence.** A call can be missed because the
  submitter used wording outside the keyword set. Two micro-class genomes lose their
  RBP call under re-annotation for exactly this reason, and this is reported in the
  manuscript rather than corrected silently.
- **The Twort TerL entry is a fragment.** Twort's terminase is split across three CDS
  by HNH endonuclease insertions; S3 extracts the longest fragment (286 aa against a
  605-aa reference), and InterPro returns no domain for that fragment. Its position in
  Figure 1 rests on partial data.
- **Internal branching within the core clade is unresolved.** Ninety-two of the 190
  TerL pairwise distances are exactly zero and bootstrap support for nodes inside the
  clade is 27–64 %. Only the separation of the core clade from Twort and from EW/SA13
  is supported at 100 %. Figure 1 should be read as a statement about lineage
  separation, not about branching order within the clade.
- **Genotype is not phenotype.** The pipeline establishes genomic capacity. Whether a
  module produces a given lytic or biofilm phenotype is a wet-lab readout outside its
  scope.

---

## 7. Licence and citation

Code released under the GNU General Public License v3.0 (see `LICENSE`). Genome records
are reproduced from NCBI GenBank and remain subject to their original terms.

Please cite the manuscript, and this archive by its Zenodo DOI.

---

## 8. References for the external tools

Bouras G, Nepal R, Houtak G, Psaltis AJ, Wormald PJ, Vreugde S. Pharokka: a fast scalable bacteriophage annotation tool. *Bioinformatics*. 2023;39(1):btac776.

Bouras G, Grigson SR, Mirdita M, Heinzinger M, Papudeshi B, Mallawaarachchi V, et al. Protein structure informed bacteriophage genome annotation with Phold. *Nucleic Acids Res*. 2026;54(1):gkaf1448.

Chan PP, Lin BY, Mak AJ, Lowe TM. tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Res*. 2021;49(16):9077–96.

Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422–3.

Katoh K, Rozewicki J, Yamada KD. MAFFT online service: multiple sequence alignment, interactive sequence choice and visualization. *Brief Bioinform*. 2019;20(4):1160–6.

Kumar S, Stecher G, Suleski M, Sanderford M, Sharma S, Tamura K. Molecular Evolutionary Genetics Analysis Version 12 for adaptive and green computing. *Mol Biol Evol*. 2024;41:1–9.

Moraru C, Varsani A, Kropinski AM. VIRIDIC — a novel tool to calculate the intergenomic similarities of prokaryote-infecting viruses. *Viruses*. 2020;12(11):1268.

Paysan-Lafosse T, Blum M, Chuguransky S, Grego T, Pinto BL, Salazar GA, et al. InterPro in 2022. *Nucleic Acids Res*. 2023;51(D1):D418–27.

Sharp PM, Tuohy TMF, Mosurski KR. Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. *Nucleic Acids Res*. 1986;14(13):5125–43.
