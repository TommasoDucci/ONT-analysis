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


def parse_wakhan_cna_vcf(vcf_path: str) -> pd.DataFrame:
    """Parse a Wakhan CNA VCF into one row per segment.

    The `ID` field is expected in the form `wakhan:<TYPE>:<chr>:<start>-<end>`,
    where `<TYPE>` is one of GAIN, LOSS, CNLOH.
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

    return pd.DataFrame(rows)


def annotate_genes_cn(
    vcf_path: str, gene_bed: str = DEFAULT_GENE_BED, output_path: str | None = None
) -> pd.DataFrame:
    """For each gene in `gene_bed`, return the overlapping CNA segment(s).

    A gene spanning more than one segment yields one row per overlapping
    segment (no information is collapsed away). The result includes a `cna`
    column summarising the call (e.g. `GAIN (TCN=3, CN1=2, CN2=1)`).

    If `output_path` is given, the result is also saved there (TSV if the
    extension is `.tsv`/`.txt`, CSV otherwise).
    """
    segments = parse_wakhan_cna_vcf(vcf_path)
    genes = load_gene_regions(gene_bed)

    rows = []
    for gene in genes.itertuples(index=False):
        overlaps = segments[
            (segments["chr"] == gene.chr)
            & (segments["start"] <= gene.end)
            & (segments["end"] >= gene.start)
        ]
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

    df = pd.DataFrame(rows)
    if output_path is not None:
        save_df(df, output_path)
    return df


def annotate_genes_cn_multi(
    vcf_paths: dict[str, str],
    gene_bed: str = DEFAULT_GENE_BED,
    output_path: str | None = None,
) -> pd.DataFrame:
    """Run `annotate_genes_cn` across multiple samples and concatenate.

    `vcf_paths` maps sample name -> VCF path. If `output_path` is given, the
    concatenated result is also saved there (TSV if the extension is
    `.tsv`/`.txt`, CSV otherwise).
    """
    dfs = []
    for sample, path in vcf_paths.items():
        df = annotate_genes_cn(path, gene_bed=gene_bed)
        df.insert(0, "sample", sample)
        dfs.append(df)

    result = pd.concat(dfs, ignore_index=True)
    if output_path is not None:
        save_df(result, output_path)
    return result


def save_df(df: pd.DataFrame, output_path: str) -> None:
    """Save `df` to `output_path`, using tab separation for `.tsv`/`.txt`."""
    sep = "\t" if output_path.lower().endswith((".tsv", ".txt")) else ","
    df.to_csv(output_path, sep=sep, index=False)


def plot_gene_cn_heatmap(df: pd.DataFrame, value: str = "tcn"):
    """Heatmap of gene x sample colored by `value` (tcn, cn1, cn2, or seg_type).

    When a gene has multiple overlapping segments for a sample, the max
    value is shown (or, for seg_type, an arbitrary but deterministic pick).
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
    )
    if value == "seg_type":
        cbar = ax.collections[0].colorbar
        cbar.set_ticks(list(categories.values()))
        cbar.set_ticklabels(labels)
    plt.xticks(rotation=60, fontsize=9)
    plt.yticks(fontsize=9)
    plt.title(value)
    plt.tight_layout()
    plt.show()
