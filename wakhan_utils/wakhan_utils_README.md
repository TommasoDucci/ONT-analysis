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

Carica un BED di coordinate geniche (`chr, start, end, gene, ...`). Vengono usate solo le prime 4
colonne: eventuali colonne aggiuntive (score, strand, ecc.) vengono lette ma scartate. Righe di
header/track non numeriche (es. una riga `track ...` o valori non numerici in `start`/`end`)
vengono rimosse automaticamente. Di default usa `DEFAULT_GENE_BED`
(`/CTGlab/projects/GBM/CEINGE/old_beds/genes_coordinates_short_robin.bed`, pannello di 13 geni
GBM-rilevanti: CDK4, CDK6, CDKN2A, CDKN2B, EGFR, FGFR1-3, MET, NTRK1-2, PDGFRA, PTEN). Passare un
path diverso per controllare altri geni.

```python
genes = load_gene_regions("my_genes.bed")
```

---

### `parse_wakhan_cna_vcf(vcf_path, sample=None)`

Parsifica un VCF di segmenti CNA prodotto da Wakhan (es. `*_wakhan_cna_integers.vcf`) e
restituisce una `DataFrame` con una riga per segmento: `chr, start, end, seg_id, seg_type`
(`GAIN`/`LOSS`/`CNLOH`, estratto dal campo `ID`), `alt` (campo ALT grezzo), e i campi FORMAT
grezzi `gt, tcn, cn1, cn2, cnq1, cnq2, cov1, cov2`. Se si passa `sample` (es. `"GSC11"`), viene
aggiunta come prima colonna una colonna `sample` con quel nome.

```python
segments = parse_wakhan_cna_vcf("Sample_wakhan_cna_integers.vcf", sample="GSC11")
```

---

### `annotate_genes_cn(vcf_path, gene_bed=DEFAULT_GENE_BED, sample=None, output_path=None, missing_genes_output_path=None)`

Per ogni gene in `gene_bed`, trova il/i segmento/i CNA che lo sovrappongono e restituisce una
`DataFrame` con una riga per ogni coppia (gene, segmento sovrapposto) — se un gene ricade a
cavallo di più segmenti, vengono restituite più righe invece di collassare l'informazione.
Include tutte le colonne di `parse_wakhan_cna_vcf` più `gene, gene_chr, gene_start, gene_end,
overlap_bp` e `cna`, una stringa riassuntiva della chiamata (es. `GAIN (TCN=3, CN1=2, CN2=1)`).

Se si passa `sample` (es. `"GSC11"`), viene aggiunta come prima colonna una colonna `sample` con
quel nome: rimane così nel file salvato tramite `output_path` e viene usata come etichetta da
`plot_gene_cn_heatmap`.

Se si passa `output_path`, il risultato completo viene anche salvato su file (TSV se l'estensione
è `.tsv`/`.txt`, altrimenti CSV) tramite `save_df`, in aggiunta al valore restituito — si può
quindi visualizzare l'output nel notebook e salvarlo su file nella stessa chiamata.

I geni senza alcun segmento CNA sovrapposto non generano righe nel risultato; i loro nomi vengono
stampati come avviso (es. `No CNA segment overlaps (GSC11): GENE1, GENE2`) e, se si passa
`missing_genes_output_path`, salvati anche su file come tabella a una colonna (`gene`, più
`sample` se passato).

```python
ann = annotate_genes_cn(
    "Sample_wakhan_cna_integers.vcf",
    sample="GSC11",
    output_path="GSC11_sample_cna.tsv",
    missing_genes_output_path="GSC11_missing_genes.tsv",
)
ann[["sample", "gene", "seg_type", "cna", "tcn", "cn1", "cn2"]]
```

---

### `annotate_genes_cn_multi(vcf_paths, gene_bed=DEFAULT_GENE_BED, output_path=None, missing_genes_output_path=None)`

Esegue `annotate_genes_cn` su più campioni (`vcf_paths`: `dict[nome_campione, path]`, es. più
soluzioni Wakhan o più campioni spaziali dello stesso paziente) e concatena i risultati con una
colonna `sample`. Come per `annotate_genes_cn`, passare `output_path` salva anche il risultato
concatenato su file, senza impedirne la visualizzazione nel notebook.

Se si passa `missing_genes_output_path`, i geni senza segmento CNA sovrapposto vengono salvati lì,
una riga per coppia (`sample`, `gene`), aggregando su tutti i campioni di `vcf_paths` (in aggiunta
all'avviso stampato da `annotate_genes_cn` per ciascun campione).

```python
df = annotate_genes_cn_multi({
    "solution_1": "solution_1/vcf_output/Sample_..._integers.vcf",
    "solution_2": "solution_2/vcf_output/Sample_..._integers.vcf",
}, output_path="multi_sample_cna.tsv", missing_genes_output_path="multi_sample_missing_genes.tsv")
```

---

### `save_df(df, output_path)`

Salva una `DataFrame` su `output_path`, usando la tabulazione come separatore se l'estensione è
`.tsv`/`.txt`, altrimenti la virgola (CSV). Usata internamente da `annotate_genes_cn` e
`annotate_genes_cn_multi` quando si passa `output_path`, ma può essere chiamata anche
direttamente su qualsiasi `DataFrame`.

```python
save_df(df, "multi_sample_cna.csv")
```

---

### `plot_gene_cn_heatmap(df, value="tcn", output_path=None)`

Heatmap gene × campione colorata per `value` (`"tcn"`, `"cn1"`, `"cn2"`, oppure `"seg_type"` come
categoria GAIN/CNLOH/LOSS). Se un gene ha più segmenti sovrapposti per lo stesso campione, viene
mostrato il valore massimo. Per `"tcn"`/`"cn1"`/`"cn2"` non viene disegnata la colorbar (ogni
cella riporta già il proprio valore); per `"seg_type"` la colorbar viene mantenuta, con le
etichette testuali `LOSS`/`CNLOH`/`GAIN` invece dei codici numerici, perché le celle non
riportano numeri.

Richiede che `df` abbia una colonna `sample` (aggiunta automaticamente da `annotate_genes_cn`/
`annotate_genes_cn_multi` quando gli si passa un nome campione): viene usata sia sull'asse x sia
nel titolo del grafico (per un singolo campione, il titolo diventa `"<value> — <sample>"`).

Se si passa `output_path`, la figura viene anche salvata su file (formato dedotto
dall'estensione, es. `.png`/`.pdf`/`.svg`), oltre a essere mostrata nel notebook.

```python
plot_gene_cn_heatmap(df, value="tcn", output_path="tcn_heatmap.png")
plot_gene_cn_heatmap(df, value="seg_type", output_path="seg_type_heatmap.png")
```

---

## Workflow tipico

```python
from wakhan_utils import annotate_genes_cn, annotate_genes_cn_multi, plot_gene_cn_heatmap

# singolo campione: visualizzato nel notebook e salvato su file nella stessa chiamata
sample = "GSC11"
ann = annotate_genes_cn(
    "solution_1/vcf_output/Sample_3.44_0.87_0.91_wakhan_cna_integers.vcf",
    sample=sample,
    output_path=f"{sample}_single_sample_cna.tsv",
    missing_genes_output_path=f"{sample}_missing_genes.tsv",
)
plot_gene_cn_heatmap(ann, value="tcn", output_path=f"{sample}_tcn_heatmap.png")

# confronto tra più soluzioni per lo stesso campione
df = annotate_genes_cn_multi({
    "solution_1.integers": "solution_1/vcf_output/Sample_3.44_0.87_0.91_wakhan_cna_integers.vcf",
    "solution_1.subclonals": "solution_1/vcf_output/Sample_3.44_0.87_0.91_wakhan_cna_subclonals.vcf",
}, output_path="multi_sample_cna.tsv", missing_genes_output_path="multi_sample_missing_genes.tsv")
plot_gene_cn_heatmap(df, value="tcn", output_path="multi_sample_tcn_heatmap.png")
```

## Notebook (`analyze_wakhan_cna.ipynb`)

Il notebook è parametrizzato per essere rieseguito su linee cellulari diverse senza modificare le
celle di analisi: nella cella "Parametri" si cambia solo `sample` (es. `"GSC11"`), e il path del
VCF di `solution_1` — il cui nome file include stime di purezza/ploidia variabili per campione,
es. `Sample_3.44_0.87_0.91_..._integers.vcf` — viene trovato automaticamente con `glob` dentro
`{base_dir}/{sample}/wakhan_cna/solution_1/vcf_output/`. Richiede quindi che ogni cartella
campione contenga esattamente un file `*_wakhan_cna_integers.vcf` e uno
`*_wakhan_cna_subclonals.vcf`.
