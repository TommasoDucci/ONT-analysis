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
    """Load a 4-column BED (chr, start, end, gene) of gene coordinates."""
    return pd.read_csv(
        bed_path, sep="\t", header=None, names=["chr", "start", "end", "gene"]
    )


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


def annotate_genes_cn(vcf_path: str, gene_bed: str = DEFAULT_GENE_BED) -> pd.DataFrame:
    """For each gene in `gene_bed`, return the overlapping CNA segment(s).

    A gene spanning more than one segment yields one row per overlapping
    segment (no information is collapsed away).
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
            rows.append({
                "gene": gene.gene,
                "gene_chr": gene.chr,
                "gene_start": gene.start,
                "gene_end": gene.end,
                "overlap_bp": overlap_bp,
                **seg._asdict(),
            })

    return pd.DataFrame(rows)


def annotate_genes_cn_multi(
    vcf_paths: dict[str, str], gene_bed: str = DEFAULT_GENE_BED
) -> pd.DataFrame:
    """Run `annotate_genes_cn` across multiple samples and concatenate.

    `vcf_paths` maps sample name -> VCF path.
    """
    dfs = []
    for sample, path in vcf_paths.items():
        df = annotate_genes_cn(path, gene_bed=gene_bed)
        df.insert(0, "sample", sample)
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


def plot_gene_cn_heatmap(df: pd.DataFrame, value: str = "tcn"):
    """Heatmap of gene x sample colored by `value` (tcn, cn1, cn2, or seg_type).

    When a gene has multiple overlapping segments for a sample, the max
    value is shown (or, for seg_type, an arbitrary but deterministic pick).
    """
    if value == "seg_type":
        categories = {"LOSS": 0, "CNLOH": 1, "GAIN": 2}
        pivot = df.pivot_table(
            index="gene", columns="sample",
            values="seg_type", aggfunc=lambda s: categories[s.iloc[0]],
        )
        cmap = sns.color_palette(["#3b4cc0", "#dddddd", "#b40426"], as_cmap=False)
    else:
        pivot = df.pivot_table(index="gene", columns="sample", values=value, aggfunc="max")
        cmap = "viridis"

    plt.figure(figsize=(max(6, pivot.shape[1] * 1.2), max(6, pivot.shape[0] * 0.4)))
    sns.heatmap(pivot, cmap=cmap, annot=True, fmt=".2g" if value != "seg_type" else "")
    plt.xticks(rotation=60, fontsize=9)
    plt.yticks(fontsize=9)
    plt.title(value)
    plt.tight_layout()
    plt.show()
