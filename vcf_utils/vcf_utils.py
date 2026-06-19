"""
Generic utilities for parsing and analysing multi-sample VCF files.
Extracted from parse_merged_vcfs.ipynb.
"""

from cyvcf2 import VCF
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from sklearn.mixture import BayesianGaussianMixture
import copy


class VariantsSummary:
    def __init__(
        self,
        path: str,
        X: np.ndarray,
        metadata: pd.DataFrame,
        samples: list,
        AF: np.ndarray | None = None,
        depth_matrix: np.ndarray | None = None,
    ):
        self.path = path
        self.X = X
        self.metadata = metadata
        self.samples = samples
        self.AF = AF
        self.depth_matrix = depth_matrix

        if self.metadata is not None and self.X.shape[0] != len(self.metadata):
            raise ValueError("X rows and metadata rows must match")

        if self.depth_matrix is not None and self.depth_matrix.shape[0] != self.X.shape[0]:
            raise ValueError("depth_matrix rows must match X rows")

    def filter(self, filters: dict[str, list]):
        mask = np.ones(len(self.metadata), dtype=bool)
        for k, allowed in filters.items():
            mask &= self.metadata[k].isin(allowed)
        mask = mask.to_numpy().astype(bool)

        self.X = self.X[mask, :]
        if self.depth_matrix is not None:
            self.depth_matrix = self.depth_matrix[mask, :]
        if self.AF is not None:
            self.AF = self.AF[mask, :]
        self.metadata = self.metadata.loc[mask].reset_index(drop=True)

    def plot_variants(self, attribute: str = "presence"):
        if attribute == "presence":
            data = self.X
            plot_value = False
        elif attribute == "frequency":
            data = self.AF
            plot_value = True
        elif attribute == "depth":
            data = self.depth_matrix
            plot_value = True
        else:
            raise ValueError(f"Unknown attribute: {attribute}")

        labels = [
            f"{self.metadata['symbol'][i]}_{self.metadata['consequence'][i]}"
            for i in range(self.X.shape[0])
        ]
        df = pd.DataFrame(data, index=labels, columns=self.samples)

        plt.figure(figsize=(10, 8))
        sns.heatmap(df, cmap="viridis", annot=plot_value)
        plt.xticks(rotation=60, fontsize=9)
        plt.yticks(fontsize=9)
        plt.title(attribute)
        plt.tight_layout()
        plt.show()
    
    def high_mod_heatmap(self, order=None, gene_order=None):
        data = copy.deepcopy(self)
        data.filter({"impact": ["HIGH", "MODERATE"], "filter": ["PASS"]})

        genes = gene_order if gene_order is not None else list(data.metadata.groupby("symbol").groups.keys())
        cols = order if order is not None else data.samples
        all_genes = list(data.metadata.groupby("symbol").groups.keys())
        df_high = pd.DataFrame(np.nan, index=all_genes, columns=data.samples)
        df_mod  = pd.DataFrame(np.nan, index=all_genes, columns=data.samples)

        for gene, df in data.metadata.groupby("symbol"):
            by_impact = df.groupby("impact")
            try:
                high = by_impact.get_group("HIGH")
                df_high.loc[gene] = data.AF[high.index, :].max(axis=0)
            except Exception:
                pass
            try:
                mod = by_impact.get_group("MODERATE")
                df_mod.loc[gene] = data.AF[mod.index, :].max(axis=0)
            except Exception:
                pass

        # reindex rows/cols to requested order, filling missing entries with 0
        df_high = df_high.reindex(index=genes, columns=cols, fill_value=np.nan)
        df_mod  = df_mod.reindex(index=genes, columns=cols, fill_value=np.nan)

        h = df_high.fillna(0).values
        m = df_mod.fillna(0).values
        af = np.where(h > 0, h, m)

        # raised baseline: HIGH → (1,0.7,0.7) light pink → red; MOD → (1,1,0.7) light lemon → yellow
        base = 0.7
        zero_mask = (af == 0)
        g = np.where(h > 0, base * (1 - af), 1.0)          # HIGH drops g; MOD keeps g=1
        b = np.where(zero_mask, 1.0, base * (1 - af))       # zero stays white; others drop b
        rgba = np.dstack([np.ones_like(af), g, b, np.ones_like(af)])

        fig, ax = plt.subplots(figsize=(max(6, len(cols) * 1.5), max(6, len(genes) * 0.5)))
        ax.imshow(rgba, aspect="auto")
        ax.set(xticks=range(len(cols)), xticklabels=cols,
               yticks=range(len(genes)), yticklabels=genes)
        ax.tick_params(axis="x", rotation=60)
        for (i, j), v in np.ndenumerate(af):
            if v > 0:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8)
        plt.tight_layout()
        plt.show()

        return df_high, df_mod



def get_VariantsSummary(vcf_path: str) -> VariantsSummary:
    """Parse a multi-sample VCF (with VEP CSQ annotation) into a VariantsSummary."""
    vcf = VCF(vcf_path)
    samples = vcf.samples

    meta = pd.DataFrame(
        columns=["chr", "pos", "end", "alt", "ref", "filter",
                 "consequence", "impact", "symbol", "gene", "exon"]
    )
    depth, X, AF = [], [], []

    for v in vcf:
        csq = v.INFO.get("CSQ")
        if csq is not None:
            csq = csq.split("|")
            consequence = csq[1].strip().upper()
            impact      = csq[2].strip().upper()
            symbol      = csq[3].strip().upper()
            gene        = csq[4].strip().upper()
            exon        = csq[8].strip().upper()
        else:
            consequence = impact = symbol = gene = exon = None

        meta.loc[len(meta)] = {
            "chr": v.CHROM, "pos": v.POS, "end": v.end,
            "alt": v.ALT,   "ref": v.REF, "filter": v.FILTERS[0],
            "consequence": consequence, "impact": impact,
            "symbol": symbol, "gene": gene, "exon": exon,
        }

        depth.append(v.gt_depths)
        AF.append(np.where(v.gt_alt_freqs < 0, 0, v.gt_alt_freqs))
        X.append(np.where(v.gt_depths < 0, 0, 1))

    return VariantsSummary(
        path=vcf_path,
        X=np.stack(X, axis=0),
        samples=samples,
        metadata=meta,
        depth_matrix=np.stack(depth, axis=0),
        AF=np.stack(AF, axis=0),
    )


def find_clonal_threshold(
    af_matrix: np.ndarray,
    n_components: int = 3,
    bins: int = 70,
    plot: bool = True,
    title: str = "",
) -> float:
    """
    Fit a Bayesian GMM on non-zero AFs and return the lower boundary of the
    second cluster (clonal/subclonal threshold).

    Parameters
    ----------
    af_matrix : np.ndarray
        2-D AF matrix (variants × samples); zeros are excluded before fitting.
    n_components : int
        Maximum number of GMM components (Dirichlet process prior).
    bins : int
        Number of histogram bins for visualization.
    plot : bool
        Whether to show the annotated histogram.
    title : str
        Plot title suffix.

    Returns
    -------
    float
        Threshold value separating the first cluster from the second.
    """
    X_flat = af_matrix.flatten()
    X_nz = X_flat[X_flat != 0].reshape(-1, 1)

    bgm = BayesianGaussianMixture(
        n_components=n_components,
        weight_concentration_prior_type="dirichlet_process",
        covariance_type="full",
        max_iter=1000,
    )
    bgm.fit(X_nz)
    probs = bgm.predict_proba(X_nz)

    counts, edges = np.histogram(X_nz, bins=bins)
    intervals = list(zip(edges[:-1], edges[1:]))
    X_ravel = X_nz.ravel()
    clusters = []
    for lo, hi in intervals:
        mask = (lo <= X_ravel) & (X_ravel < hi)
        p = probs[mask, :]
        clusters.append(int(np.argmax(np.sum(p, axis=0))) if p.shape[0] > 0 else -1)

    clusters = np.array(clusters)
    breaks = sorted([
        edges[np.min(np.where(clusters == s)[0])]
        for s in set(clusters)
        if s >= 0 and np.any(clusters == s)
    ])
    threshold = breaks[1] if len(breaks) > 1 else breaks[0]

    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))
        cmap = plt.get_cmap("tab10")
        _, _, patches = ax.hist(X_ravel, bins=edges)
        for j, patch in enumerate(patches):
            k = clusters[j]
            patch.set_facecolor(cmap(k % 10) if k >= 0 else "lightgrey")
        ax.vlines(threshold, linestyles=":", colors="black",
                  ymin=0, ymax=ax.get_ylim()[1])
        ax.set_title(title)
        plt.tight_layout()
        plt.show()

    return threshold


def plot_clonality_summary(
    summary: VariantsSummary,
    threshold: float,
    name: str = "",
):
    """
    Barplot (clonal vs subclonal counts per sample) and AF boxplot
    for a VariantsSummary object given a pre-computed threshold.
    """
    df_af = pd.DataFrame(summary.AF, columns=summary.samples).melt(
        var_name="sample", value_name="allele frequency"
    )
    df_af = df_af[df_af["allele frequency"] > 0].copy()
    df_af["clonality"] = np.where(
        df_af["allele frequency"] > threshold, "clonal_count", "subclonal_count"
    )

    df_summary = (
        df_af.groupby(["sample", "clonality"])["allele frequency"]
        .size()
        .reset_index(name="count")
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    sns.barplot(data=df_summary, x="sample", y="count", hue="clonality", ax=axes[0])
    axes[0].set_title(f"Per-spot clonal vs subclonal counts {name}")
    axes[0].tick_params(axis="x", rotation=90)

    sns.boxplot(data=df_af, x="sample", y="allele frequency", ax=axes[1], showfliers=False)
    axes[1].set_title(f"Per-spot AF boxplot {name}")
    axes[1].tick_params(axis="x", rotation=90)
    axes[1].hlines(threshold, xmin=axes[1].get_xlim()[0],
                   xmax=axes[1].get_xlim()[1], colors="red", linestyles=":")

    plt.tight_layout()
    plt.show()

    return df_af, df_summary


def plot_sample_clustering(
    summary: VariantsSummary,
    metric: str = "jaccard",
    method: str = "average",
    title: str = "",
):
    """
    Hierarchical clustering of samples based on variant presence/absence.
    Draws a clustermap and an optional dendrogram.
    """
    X = summary.X.T
    samples = np.array(summary.samples)

    D = pdist(X, metric=metric)
    Z = linkage(D, method=method)
    order = leaves_list(Z)

    df = pd.DataFrame(data=X[order], index=samples[order])

    sns.clustermap(df, cmap="viridis", row_cluster=True, col_cluster=False,
                   figsize=(10, 4))
    plt.title(title)
    plt.show()

    plt.figure(figsize=(10, 3))
    dendrogram(Z, labels=samples[order])
    plt.xticks(fontsize=7)
    plt.title(f"Dendrogram {title}")
    plt.show()

