# wakhan_utils.py

Utilities per il parsing dei VCF di copy number di Wakhan e il controllo dello
stato del copy number a livello di singolo gene.

## Dipendenze

```
cyvcf2, pandas, numpy, matplotlib, seaborn
```

## Utilizzo rapido

```python
from wakhan_utils import annotate_genes_cn, annotate_genes_cn_multi, plot_gene_cn_heatmap
```

---

## API

### `load_gene_regions(bed_path=DEFAULT_GENE_BED)`

Carica un BED a 4 colonne (`chr, start, end, gene`) di coordinate geniche. Di default usa
`DEFAULT_GENE_BED` (`/CTGlab/projects/GBM/CEINGE/old_beds/genes_coordinates_short_robin.bed`,
pannello di 13 geni GBM-rilevanti: CDK4, CDK6, CDKN2A, CDKN2B, EGFR, FGFR1-3, MET, NTRK1-2,
PDGFRA, PTEN). Passare un path diverso per controllare altri geni.

```python
genes = load_gene_regions("my_genes.bed")
```

---

### `parse_wakhan_cna_vcf(vcf_path)`

Parsifica un VCF di segmenti CNA prodotto da Wakhan (es. `*_wakhan_cna_integers.vcf`) e
restituisce una `DataFrame` con una riga per segmento: `chr, start, end, seg_id, seg_type`
(`GAIN`/`LOSS`/`CNLOH`, estratto dal campo `ID`), `alt` (campo ALT grezzo), e i campi FORMAT
grezzi `gt, tcn, cn1, cn2, cnq1, cnq2, cov1, cov2`.

```python
segments = parse_wakhan_cna_vcf("Sample_wakhan_cna_integers.vcf")
```

---

### `annotate_genes_cn(vcf_path, gene_bed=DEFAULT_GENE_BED)`

Per ogni gene in `gene_bed`, trova il/i segmento/i CNA che lo sovrappongono e restituisce una
`DataFrame` con una riga per ogni coppia (gene, segmento sovrapposto) — se un gene ricade a
cavallo di più segmenti, vengono restituite più righe invece di collassare l'informazione.
Include tutte le colonne di `parse_wakhan_cna_vcf` più `gene, gene_chr, gene_start, gene_end,
overlap_bp`.

```python
ann = annotate_genes_cn("Sample_wakhan_cna_integers.vcf")
ann[["gene", "seg_type", "tcn", "cn1", "cn2"]]
```

---

### `annotate_genes_cn_multi(vcf_paths, gene_bed=DEFAULT_GENE_BED)`

Esegue `annotate_genes_cn` su più campioni (`vcf_paths`: `dict[nome_campione, path]`, es. più
soluzioni Wakhan o più campioni spaziali dello stesso paziente) e concatena i risultati con una
colonna `sample`.

```python
df = annotate_genes_cn_multi({
    "solution_1": "solution_1/vcf_output/Sample_..._integers.vcf",
    "solution_2": "solution_2/vcf_output/Sample_..._integers.vcf",
})
```

---

### `plot_gene_cn_heatmap(df, value="tcn")`

Heatmap gene × campione colorata per `value` (`"tcn"`, `"cn1"`, `"cn2"`, oppure `"seg_type"` come
categoria GAIN/CNLOH/LOSS). Se un gene ha più segmenti sovrapposti per lo stesso campione, viene
mostrato il valore massimo.

```python
plot_gene_cn_heatmap(df, value="tcn")
plot_gene_cn_heatmap(df, value="seg_type")
```

---

## Workflow tipico

```python
from wakhan_utils import annotate_genes_cn, annotate_genes_cn_multi, plot_gene_cn_heatmap

# singolo campione
ann = annotate_genes_cn("solution_1/vcf_output/Sample_3.44_0.87_0.91_wakhan_cna_integers.vcf")

# confronto tra più soluzioni/campioni
df = annotate_genes_cn_multi({
    "solution_1": "solution_1/vcf_output/Sample_3.44_0.87_0.91_wakhan_cna_integers.vcf",
    "solution_2": "solution_2/vcf_output/Sample_1.43_0.44_0.91_wakhan_cna_integers.vcf",
})
plot_gene_cn_heatmap(df, value="tcn")
```
