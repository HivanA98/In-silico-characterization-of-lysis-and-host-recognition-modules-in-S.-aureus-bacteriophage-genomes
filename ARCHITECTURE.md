# ARCHITECTURE — Phage Characterization Toolkit v3.2

How the toolkit is layered, and the isolation contract that keeps a host's
calibration auditable. This document describes the design of the code deposited in
this archive; it is not a development roadmap.

---

## 1. One engine, a profile per host

The code analyses the **phage** genome; host biology supplies only configuration and
curation. Parsing, QC, deduplication, taxonomy, manifests and every detection
mechanism are identical regardless of which bacterial host is being studied. The
toolkit is therefore **one shared engine (`phagecore`) plus a host profile** — never
forked codebases.

A profile carries two kinds of information:

- **CONFIG** — simple values: size, GC and CDS bounds, endolysin length windows,
  taxonomy overrides, tRNA anchors, host codon reference, biofilm matrix.
- **CONTENT** — authored curation: matrix-depolymerase terms, host-specific
  receptor-binding vocabulary, precision denylists, per-host annotation-gap lists, and
  structurally confirmed `depolymerase_known_cases` with source tags.

The **universal invariants live once in the engine** so they cannot drift between
hosts: the Ile2-CAT→ATA lysidine rule and the RSCU family-sum arithmetic (S7), the CAT
Met/Ile2 ambiguity and pseudogene handling (S6), the InterPro reconciliation logic and
the chaperone precision guard (S5).

**This archive contains the *Staphylococcus aureus* profile only.** The
engine/profile split is described here because it determines how the code is
organised and how `--profile` behaves, not because other profiles are included.

---

## 2. The tiers

| Tier | Scripts | Engine module | To change host |
|---|---|---|---|
| **1 — engine + profile** | S1, S2, S4 | `phagecore` (qc, lysis, …) | `--profile <host>` |
| **2 — host-agnostic** | S3 (TerL) | none (the terminase is universal) | nothing |
| **3 — engine + profile** | S5, S6, S7 | `phagecore.rbp / .trna / .codon` | `--profile <host>` |

`--profile` is **required** on S1, S2, S4, S5, S6 and S7; `load_profile(None)` raises.
S3 takes no profile because the terminase large subunit is universal.

S4b and S7b are companion consumers rather than stages in their own right: S4b reads
the InterPro output that S4 already produced, and S7b reuses the RSCU functions from
`phagecore.codon` that S7 uses. Neither reimplements anything, so their definitions
cannot diverge from the stage they extend.

---

## 3. Why vocabulary belongs in the profile, not the engine

Any word list baked into the analysis code silently encodes an assumption about how a
particular group of submitters names its proteins. Three separate defects in
`phagecore/rbp.py` shared that single root cause:

1. **`profile.rbp_keywords` was never read.** S5 used only its module-level
   `PRIMARY_RBP_KEYWORDS` / `EXTENDED_RBP_KEYWORDS`. The profile field existed and was
   ignored.
2. **The InterPro submission filter** admitted only `confirmed`/`high` confidence
   candidates. Anything found through a profile vocabulary arrives as `putative` by
   construction, so it was filtered out before reaching InterPro.
3. **The large-carrier fold rescue** re-checked product names against another
   hardcoded list (`tail fiber`, `tailspike`, `spike`, `receptor`).

The fixes are all **additive and host-neutral**:

- `configure()` loads `profile.rbp_keywords` into `_HOST_RBP_TERMS`, filtered against
  the universal lists, and the extended scan iterates both. Verified empty for the
  previously validated configuration, so it cannot regress.
- The submission filter always includes carriers ≥ `LARGE_CARRIER_SUBMIT_AA` (700 aa)
  regardless of confidence tier.
- The rescue filters on length only — membership in `all_candidates` already means
  "this is a receptor-binding candidate".

**Design rule:** any vocabulary that describes *how proteins are named* belongs in the
profile. Only mechanisms that are true for all phages belong in the engine.

The same principle is what makes Table III of the manuscript interpretable: because the
detection vocabulary and the funnel are held fixed in the profile and the engine, the
only variable between the three annotation arms is the source of the CDS product names.

---

## 4. Three screens, cheapest resource first

Structure prediction is the scarce resource, so the pipeline spends CPU before GPU.

```
S4 lysis-axis triage ──► re-annotation ──► re-run S1–S4
S5 RBP-axis triage   ──► re-annotation ──► re-run S5
S5 fold targets      ──► ColabFold/AlphaFold2 ──► Dali
```

- **Lysis-axis triage** (S4): zero CDS on a genome ≥10 kb; a structural module with no
  lysis gene; accessory-only annotation; no features at all.
- **RBP-axis triage** (S5, after InterPro): a large carrier with neither an
  enzyme-class product nor a depolymerase verdict; or every RBP candidate returning no
  domain.
- **Fold targets** (S5): deduplicated, ≤1200 aa, clean headers, and proteins already in
  `depolymerase_known_cases` excluded so confirmed work is never repeated. This is what
  makes structural evidence **cumulative across versions** instead of being re-derived
  each release.

Only the lysis-axis triage bears on the manuscript; the S5 screens are included for
completeness and produce no manuscript output.

---

## 5. Evidence tiers and the curated override

`depolymerase_known_cases` is a curated registry, not a hardcode: each entry carries
`accession`, `protein_id`, `verdict`, `evidence` (Dali Z, % identity, PDB) and
`source`. `configure()` loads it and the reconciliation applies it as an override with
provenance, because a depolymerase **domain** inside a 2700–3000 aa carrier is
invisible to both a keyword scan and a whole-protein InterPro verdict.

Two verdict tiers:

- `matrix_depolymerase` — the enzyme class matches the host's own matrix polymer.
- `depolymerase` — the fold is Dali-confirmed but the substrate is not
  host-matrix-specific.

Candidates below the confirmation floor are written into the profile as **excluded,
with the reason**. An absence that comes from an incomplete fold run is recorded as an
**incomplete negative** and must not be read as "no depolymerase is present".

The same override mechanism carries the lysis-axis `known_cases`, which record the
eight genomes where the deposited annotation alone gives the wrong endolysin answer:
intron-split `lysK.1`, HNH-disrupted ORFs, tBLASTn-recovered endolysins, one
tail-anchored VAPH, and one divergent endolysin. Because these are applied as an
override, the corresponding entries in the audit carry the curated note rather than a
fresh tBLASTn result; the note names its source.

---

## 6. Isolation contract

- **Engine shared, stateless per run.** Host content enters through the profile at load
  time; each CLI invocation is one process with one profile.
- **No silent default.** A missing `--profile` is an error, never a fallback.
- **Profiles are independent.** Keyword sets, matrix terms, vocabulary and
  `known_cases` are never inherited.
- **`validate_outputs.py` reads outputs only** and never touches the scripts.
- **Each stage owns its outputs**; no stage writes into another's files.
- **Two InterPro result folders** (`interpro_endolysin/`, `interpro_rbp/`) keep the
  lysis and host-recognition axes from cross-contaminating.
- **One genome per file.** `phagecore/genbank_io.py` reads the first record of each
  file. A multi-record file is truncated to its first genome without error, which is
  why the Phynteny output in this archive must be split before use.

---

## 7. Calibration workflow

1. Derive size, GC, CDS and endolysin bounds **empirically**; re-derive on every change
   to the genome set. A stale bound either floods the QC report or silently stops
   flagging anything.
2. Check the naming conventions of the genome set and extend `rbp_keywords` if they
   differ from the universal vocabulary. Profile keywords are additive, so they can be
   added without disturbing a validated configuration.
3. Run both triage loops and re-annotate what they flag before trusting any extraction.
4. Set the S6 tRNA anchor from a well-annotated phage and the S7 codon reference from
   the host genome. Never assume a signature.
5. Fold only the selected targets, compare with Dali, and record confirmed cases with
   the correct tier and a source tag.

Step 1 makes S1–S4 usable immediately. Steps 4 and 5 are what make S5–S7 trustworthy,
and they are calibration work that cannot be inherited.

---

## 8. What the architecture does not do

The pipeline establishes **genomic capacity**. Whether an endolysin lyses a given
strain, whether a depolymerase reduces or reinforces biofilm, and an esterase's exact
substrate are wet-lab readouts outside its scope. A domain call reports fold or
enzyme-class capacity confirmed by InterPro, or by InterPro→Dali concordance — never a
predicted phenotype.

Equally, S1–S4 are extraction and curation layers over an existing annotation, not
independent gene callers. Their output represents what the annotation contains, which
is why the manuscript reports the same funnel run under three annotation sources rather
than a single set of counts.
