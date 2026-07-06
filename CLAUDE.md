# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository purpose

This repo holds analysis code for Oxford Nanopore (ONT) sequencing data, currently focused on
multi-sample VCF parsing and clonality/heterogeneity analysis (e.g. tumor sample comparisons).

- `vcf_utils/` — Python utilities for parsing and analyzing multi-sample VCFs annotated with VEP.
  See `vcf_utils/vcf_utils_README.md` for the full API reference (in Italian) and
  `vcf_utils/vcf_utils.py` for the implementation. `analyze_vcfs.ipynb` is the notebook that
  drives this code interactively; `vcf_utils.py` was extracted from an earlier
  `parse_merged_vcfs.ipynb` notebook.
- `Lumos/` — currently empty, reserved for future work.

There is no build system, package manifest, test suite, or linter configured in this repo yet —
work here is exploratory/analysis-driven, primarily through Jupyter notebooks calling into
`vcf_utils.py`.

## Dependencies

`vcf_utils.py` requires: `cyvcf2`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`,
`scikit-learn`. No requirements file exists yet — install these manually if working in this area.

## Architecture: vcf_utils.py

The module centers on a single container class, `VariantsSummary`, produced by
`get_VariantsSummary(vcf_path)`:

- Parses a multi-sample VCF (via `cyvcf2`) whose INFO field contains a VEP `CSQ` annotation
  string (pipe-delimited; indices 1/2/3/4/8 are read as consequence/impact/symbol/gene/exon —
  this positional indexing assumes a specific VEP annotation field order and will break silently
  if the CSQ format differs).
- Builds three aligned `variants × samples` matrices: `X` (presence/absence), `AF` (allele
  frequency), `depth_matrix` (read depth) — plus a `metadata` DataFrame (one row per variant,
  same row order/count as the matrices).
- All matrix-shaped attributes (`X`, `AF`, `depth_matrix`) and `metadata` must stay row-aligned;
  `filter()` mutates all of them together in place using a shared boolean mask.

Downstream analysis functions operate on `VariantsSummary` objects rather than raw VCFs:

- `find_clonal_threshold(af_matrix, ...)` — fits a Bayesian Gaussian Mixture Model on non-zero
  AFs to separate clonal from subclonal variants, returning a threshold AF value.
- `plot_clonality_summary(summary, threshold, ...)` — visualizes clonal/subclonal counts and AF
  distribution per sample using a previously computed threshold.
- `plot_sample_clustering(summary, ...)` — hierarchical clustering of samples by variant
  presence/absence (Jaccard distance by default).
- `VariantsSummary.high_mod_heatmap()` — deep-copies the summary, filters to
  HIGH/MODERATE-impact PASS variants, and renders a custom RGBA heatmap (red=HIGH, yellow=MODERATE)
  of max AF per gene per sample.

Typical workflow: `get_VariantsSummary` → `find_clonal_threshold` → `plot_clonality_summary` →
`VariantsSummary.filter(...)` → `plot_sample_clustering` / `high_mod_heatmap`.
