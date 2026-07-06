"""
Utilities for parsing and analysing Wakhan CNA VCF files and checking
gene-level copy-number status against a gene coordinate BED.
"""

from cyvcf2 import VCF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

DEFAULT_GENE_BED = "/CTGlab/projects/GBM/CEINGE/old_beds/genes_coordinates_short_robin.bed"


def load_gene_regions(bed_path: str = DEFAULT_GENE_BED) -> pd.DataFrame:
    """Load a BED of gene coordinates (chr, start, end, gene, ...).

    Only the first four columns are used; any additional BED columns
    (score, strand, etc.) are read but discarded.
    """
    bed = pd.read_csv(bed_path, sep="\t", header=None, comment="#")
    bed = bed.iloc[:, :4]
    bed.columns = ["chr", "start", "end", "gene"]
    bed = bed[bed["chr"].astype(str).str.lower() != "track"]
    bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
    bed["end"] = pd.to_numeric(bed["end"], errors="coerce")
    bed = bed.dropna(subset=["start", "end"]).astype({"start": "int64", "end": "int64"})
    return bed.reset_index(drop=True)


def parse_wakhan_cna_vcf(vcf_path: str, sample: str | None = None) -> pd.DataFrame:
    """Parse a Wakhan CNA VCF into one row per segment.

    The `ID` field is expected in the form `wakhan:<TYPE>:<chr>:<start>-<end>`,
    where `<TYPE>` is one of GAIN, LOSS, CNLOH.

    If `sample` is given, a `sample` column with that name is added as the
    first column of the result.
    """
    vcf = VCF(vcf_path)
    rows = []
    for v in vcf:
        seg_type = v.ID.split(":")[1] if v.ID is not None else None
        fmt = v.FORMAT

        def get(field):
            return v.format(field)[0][0] if field in fmt else None

        rows.append({
            "chr": v.CHROM,
            "start": v.start,
            "end": v.INFO.get("END"),
            "seg_id": v.ID,
            "seg_type": seg_type,
            "alt": ",".join(str(a) for a in v.ALT),
            "gt": "/".join(str(a) for a in v.genotypes[0][:-1]) if v.genotypes else None,
            "tcn": get("TCN"),
            "cn1": get("CN1"),
            "cn2": get("CN2"),
            "cnq1": get("CNQ1"),
            "cnq2": get("CNQ2"),
            "cov1": get("COV1"),
            "cov2": get("COV2"),
        })

    df = pd.DataFrame(rows)
    if sample is not None:
        df.insert(0, "sample", sample)
    return df


def annotate_genes_cn(
    vcf_path: str,
    gene_bed: str = DEFAULT_GENE_BED,
    sample: str | None = None,
    output_path: str | None = None,
    missing_genes_output_path: str | None = None,
) -> pd.DataFrame:
    """For each gene in `gene_bed`, return the overlapping CNA segment(s).

    A gene spanning more than one segment yields one row per overlapping
    segment (no information is collapsed away). The result includes a `cna`
    column summarising the call (e.g. `GAIN (TCN=3, CN1=2, CN2=1)`).

    If `sample` is given (e.g. `"GSC11"`), a `sample` column with that name
    is added as the first column, so it carries through to `output_path` and
    to downstream plots (`plot_gene_cn_heatmap`).

    If `output_path` is given, the result is also saved there (TSV if the
    extension is `.tsv`/`.txt`, CSV otherwise).

    Genes with no overlapping CNA segment are dropped from the result (no
    row is produced for them); their names are printed as a warning, and,
    if `missing_genes_output_path` is given, also saved there as a one-column
    (`gene`, plus `sample` if given) table.
    """
    segments = parse_wakhan_cna_vcf(vcf_path)
    genes = load_gene_regions(gene_bed)

    rows = []
    missing_genes = []
    for gene in genes.itertuples(index=False):
        overlaps = segments[
            (segments["chr"] == gene.chr)
            & (segments["start"] <= gene.end)
            & (segments["end"] >= gene.start)
        ]
        if overlaps.empty:
            missing_genes.append(gene.gene)
            continue
        for seg in overlaps.itertuples(index=False):
            overlap_bp = min(gene.end, seg.end) - max(gene.start, seg.start)
            seg_dict = seg._asdict()
            rows.append({
                "gene": gene.gene,
                "gene_chr": gene.chr,
                "gene_start": gene.start,
                "gene_end": gene.end,
                "overlap_bp": overlap_bp,
                **seg_dict,
                "cna": (
                    f"{seg_dict['seg_type']} "
                    f"(TCN={seg_dict['tcn']}, CN1={seg_dict['cn1']}, CN2={seg_dict['cn2']})"
                ),
            })

    if missing_genes:
        label = f" ({sample})" if sample is not None else ""
        print(f"No CNA segment overlaps{label}: {', '.join(missing_genes)}")
        if missing_genes_output_path is not None:
            missing_df = pd.DataFrame({"gene": missing_genes})
            if sample is not None:
                missing_df.insert(0, "sample", sample)
            save_df(missing_df, missing_genes_output_path)

    df = pd.DataFrame(rows)
    if sample is not None:
        df.insert(0, "sample", sample)
    if output_path is not None:
        save_df(df, output_path)
    return df


def annotate_genes_cn_multi(
    vcf_paths: dict[str, str],
    gene_bed: str = DEFAULT_GENE_BED,
    output_path: str | None = None,
    missing_genes_output_path: str | None = None,
) -> pd.DataFrame:
    """Run `annotate_genes_cn` across multiple samples and concatenate.

    `vcf_paths` maps sample name -> VCF path. If `output_path` is given, the
    concatenated result is also saved there (TSV if the extension is
    `.tsv`/`.txt`, CSV otherwise).

    If `missing_genes_output_path` is given, genes with no overlapping CNA
    segment are saved there as one row per (sample, gene), across all
    samples in `vcf_paths`.
    """
    dfs = [
        annotate_genes_cn(path, gene_bed=gene_bed, sample=sample)
        for sample, path in vcf_paths.items()
    ]

    result = pd.concat(dfs, ignore_index=True)
    if output_path is not None:
        save_df(result, output_path)

    if missing_genes_output_path is not None:
        all_genes = load_gene_regions(gene_bed)["gene"]
        missing_rows = [
            {"sample": sample, "gene": gene}
            for sample, df in zip(vcf_paths, dfs)
            for gene in all_genes[~all_genes.isin(df["gene"])]
        ]
        save_df(pd.DataFrame(missing_rows), missing_genes_output_path)

    return result


def save_df(df: pd.DataFrame, output_path: str) -> None:
    """Save `df` to `output_path`, using tab separation for `.tsv`/`.txt`."""
    sep = "\t" if output_path.lower().endswith((".tsv", ".txt")) else ","
    df.to_csv(output_path, sep=sep, index=False)


def plot_gene_cn_heatmap(df: pd.DataFrame, value: str = "tcn", output_path: str | None = None):
    """Heatmap of gene x sample colored by `value` (tcn, cn1, cn2, or seg_type).

    When a gene has multiple overlapping segments for a sample, the max
    value is shown (or, for seg_type, an arbitrary but deterministic pick).

    `df` must have a `sample` column (added automatically by
    `annotate_genes_cn`/`annotate_genes_cn_multi` when a `sample` name is
    passed to them); it is used both as the x-axis and in the plot title.

    For `"tcn"`/`"cn1"`/`"cn2"`, no colorbar is drawn (each cell is already
    annotated with its value); for `"seg_type"` a legend-style colorbar with
    `LOSS`/`CNLOH`/`GAIN` labels is kept, since the cells carry no numbers.

    If `output_path` is given, the figure is also saved there (format
    inferred from the extension, e.g. `.png`/`.pdf`/`.svg`).
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm

    if value == "seg_type":
        categories = {"LOSS": 0, "CNLOH": 1, "GAIN": 2}
        labels = list(categories.keys())
        pivot = df.pivot_table(
            index="gene", columns="sample",
            values="seg_type", aggfunc=lambda s: categories[s.iloc[0]],
        )
        cmap = ListedColormap(["#3b4cc0", "#dddddd", "#b40426"])
        norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    else:
        pivot = df.pivot_table(index="gene", columns="sample", values=value, aggfunc="max")
        cmap = "viridis"
        norm = None

    plt.figure(figsize=(max(6, pivot.shape[1] * 1.2), max(6, pivot.shape[0] * 0.4)))
    ax = sns.heatmap(
        pivot, cmap=cmap, norm=norm,
        annot=value != "seg_type", fmt=".2g",
        cbar=value == "seg_type",
    )
    if value == "seg_type":
        cbar = ax.collections[0].colorbar
        cbar.set_ticks(list(categories.values()))
        cbar.set_ticklabels(labels)
    plt.xticks(rotation=60, fontsize=9)
    plt.yticks(fontsize=9)
    samples = df["sample"].unique()
    title = f"{value} — {samples[0]}" if len(samples) == 1 else value
    plt.title(title)
    plt.tight_layout()
    if output_path is not None:
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
