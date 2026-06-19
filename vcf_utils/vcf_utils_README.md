# vcf_utils.py

Utilities per il parsing e l'analisi di VCF multi-campione annotati con VEP.
Estratte da `parse_merged_vcfs.ipynb`.

## Dipendenze

```
cyvcf2, pandas, numpy, matplotlib, seaborn, scipy, scikit-learn
```

## Utilizzo rapido

```python
from vcf_utils import get_VariantsSummary, find_clonal_threshold, plot_clonality_summary, plot_sample_clustering
```

---

## API

### `get_VariantsSummary(vcf_path)`

Parsifica un VCF multi-campione (con annotazione VEP `CSQ`) e restituisce un oggetto `VariantsSummary`.

```python
summary = get_VariantsSummary("sample_merged.vcf.gz")
```

---

### `VariantsSummary`

Oggetto container con i seguenti attributi:

| Attributo | Tipo | Descrizione |
|---|---|---|
| `path` | `str` | Path del VCF sorgente |
| `X` | `ndarray (variants × samples)` | Matrice di presenza/assenza (0/1) |
| `AF` | `ndarray (variants × samples)` | Allele frequencies |
| `depth_matrix` | `ndarray (variants × samples)` | Read depth per variante/campione |
| `metadata` | `DataFrame` | chr, pos, gene, consequence, impact, exon, ecc. |
| `samples` | `list[str]` | Nomi dei campioni |

**Metodi:**

`filter(filters: dict[str, list])` — filtra le varianti in-place in base ai valori dei metadati.
```python
summary.filter({"impact": ["HIGH", "MODERATE"]})
summary.filter({"symbol": ["TP53", "EGFR", "IDH1"]})
```

`plot_variants(attribute)` — heatmap varianti × campioni. `attribute` può essere `"presence"` (default), `"frequency"`, o `"depth"`.
```python
summary.plot_variants("frequency")
```

---

### `find_clonal_threshold(af_matrix, n_components=3, bins=70, plot=True, title="")`

Fitta un Bayesian Gaussian Mixture Model sulle AF non-zero per separare varianti clonali da subclonali. Restituisce il valore soglia (lower boundary del secondo cluster).

```python
thr = find_clonal_threshold(summary.AF, title="GB01")
```

---

### `plot_clonality_summary(summary, threshold, name="")`

Barplot (clonal vs subclonal per campione) e boxplot delle AF con soglia visualizzata. Restituisce `(df_af, df_summary)`.

```python
df_af, df_counts = plot_clonality_summary(summary, threshold=thr, name="GB01")
```

---

### `plot_sample_clustering(summary, metric="jaccard", method="average", title="")`

Clustering gerarchico dei campioni basato sulla matrice di presenza/assenza. Produce una clustermap e un dendrogramma.

```python
plot_sample_clustering(summary, title="GB01")
```

---

## Workflow tipico

```python
from vcf_utils import *

summary = get_VariantsSummary("GB01_merged.vcf.gz")

# analisi clonalità sull'intero set
thr = find_clonal_threshold(summary.AF, title="GB01")
plot_clonality_summary(summary, thr, name="GB01")

# filtraggio per geni/impatto
summary.filter({"impact": ["HIGH", "MODERATE"]})

# clustering campioni
plot_sample_clustering(summary, title="GB01 filtered")
```
