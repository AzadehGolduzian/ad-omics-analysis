# Genomic Annotation

Functional annotation pipeline for large-scale variant classification in population genomics studies (UK Biobank, All of Us).

## Files

### `cisencode_transeqtl.qmd`
**Goal:** Classify array-genotyped variants into functional regulatory categories for downstream GWAS annotation.

**Three output categories:**
- **CIS** — SNPs overlapping ENCODE candidate cis-Regulatory Elements (cCREs, GRCh38 Registry V3)
- **TRANS** — SNPs within ±0.5 Mb windows around GTEx trans-eQTL loci
- **REG** — Union of CIS and TRANS (used as the full regulatory variant set)

**Pipeline steps:**
1. Download ENCODE cCREs (GRCh38) from the Weng Lab registry
2. Load and clean cCRE BED file — filter to autosomes (chr1–22)
3. Load array variants from PLINK `.bim` file
4. Load pre-built trans-eQTL windows (±0.5 Mb)
5. Interval overlap: array SNPs ∩ ENCODE cCREs → CIS set
6. Interval overlap: array SNPs ∩ trans-eQTL windows → TRANS set
7. Union → REG set; summary counts
8. Save CIS, TRANS, REG ID lists and full annotation table (`.csv.gz`)

**Tools:** `data.table::foverlaps()` for fast genomic interval overlaps

**Environment:** Runs on HPC cluster (Jupyter-based) with PLINK `.bim` files. Not executable locally without input data.

---

> Part of a larger genomic annotation project for Chapter 2 of my dissertation on regulatory variant classification across ancestries.
