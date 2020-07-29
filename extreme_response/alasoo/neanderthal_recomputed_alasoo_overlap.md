# Explore overlap between Neanderthal-introgressed SNPs and recomputed Alasoo eQTLs from EMBL-EBI eQTL catalogue

The Alasoo *et al.* (2018) study contains macrophage eQTL for the 4 conditions: IFNg, Salmonella, IFNg + Salmonella, and Naive. This document contains code to explore if any Neanderthal-introgressed SNPs are present in the list of recomputed eQTLs for Alasoo *et al.* (2018) from the EMBL-EBI eQTL catalogue. Overlapping recomputed eQTLs were filtered for p-value < 10<sup>-8</sup> and saved in CSV format. Gene symbols of genes associated with these eQTLs were obtained and saved to text files.

First, recomputed eQTLs for Alasoo *et al.* (2018) were downloaded from the EMBL-EBI eQTL catalogue [ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Alasoo_2018/ge/] and saved under: `/well/jknight/shiyao/data/alasoo/`.  

The recomputed Fairfax eQTLs and Neanderthal-introgressed SNPs were compared to look for overlap. Overlapping eQTLs were filtered for p-value < 10<sup>-8</sup> and saved in CSV format under: `/well/jknight/shiyao/data/alasoo/overlap_all_filtered.csv`. The number of eQTLs expressed per condition was computed and shown in Table 1.

```python
# Import modules
import pandas as pd

# Load Neanderthal SNPs with hg38 coordinates
neanderthal = pd.read_csv('/well/jknight/shiyao/data/alasoo/comparison_df_hg38.csv')

# Overlap between recomputed Alasoo eQTLs and Neanderthal SNPs
# IFNg+Salmonella
IS = pd.read_csv('/well/jknight/shiyao/data/alasoo/Alasoo_2018_ge_macrophage_IFNg+Salmonella.all.tsv.gz', sep='\t', compression='gzip')
IS.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_IS = neanderthal.merge(IS, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_IS.to_csv('/well/jknight/shiyao/data/alasoo/overlap_IfnSal.csv', index=False)

# IFNg
IFN = pd.read_csv('/well/jknight/shiyao/data/alasoo/Alasoo_2018_ge_macrophage_IFNg.all.tsv.gz', sep='\t', compression='gzip')
IFN.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_IFN = neanderthal.merge(IFN, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_IFN.to_csv('/well/jknight/shiyao/data/alasoo/overlap_IFN.csv', index=False)

# Salmonella
SAL = pd.read_csv('/well/jknight/shiyao/data/alasoo/Alasoo_2018_ge_macrophage_Salmonella.all.tsv.gz', sep='\t', compression='gzip')
SAL.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_SAL = neanderthal.merge(SAL, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_SAL.to_csv('/well/jknight/shiyao/data/alasoo/overlap_salmonella.csv', index=False)

# Naive
Naive = pd.read_csv('/well/jknight/shiyao/data/alasoo/Alasoo_2018_ge_macrophage_naive.all.tsv.gz', sep='\t', compression='gzip')
Naive.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_Naive = neanderthal.merge(Naive, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_Naive.to_csv('/well/jknight/shiyao/data/alasoo/overlap_naive.csv', index=False)

# Concatenate dfs
IS = pd.read_csv('overlap_IfnSal.csv')
IS['Condition'] = 'IFNg+Salmonella'
IFN = pd.read_csv('overlap_IFN.csv')
IFN['Condition'] = 'IFNg'
SAL = pd.read_csv('overlap_salmonella.csv')
SAL['Condition'] = 'Salmonella'
Naive = pd.read_csv('overlap_naive.csv')
Naive['Condition'] = 'Naive'
combined = pd.concat([IS, IFN, SAL, Naive])
combined.to_csv('/well/jknight/shiyao/data/alasoo/overlap_all.csv', index=False)
combined['rsid'].nunique()  # 170225

# Keep only eQTLs with p < 10e-8
filtered = combined.loc[combined['pvalue'] < 10e-8].copy()
filtered.to_csv('/well/jknight/shiyao/data/alasoo/overlap_all_filtered.csv', index=False)
filtered['rsid'].nunique()  # 2041

# Get number of eQTLs per condition
condition = filtered[['rsid', 'Condition']]
condition = condition.groupby(['rsid']).agg(lambda col: ','.join(col))
for i, row in condition.iterrows():
    condition.at[i, 'Condition'] = ','.join(set(condition.at[i, 'Condition'].split(',')))
condition.Condition.str.split(',', expand=True).stack().value_counts()
```

2,041 recomputed eQTLs (p < 10<sup>-8</sup>) were also identified as Neanderthal-introgressed SNPs.

**Table 1:** Number of condition-specific eQTLs per condition
| Condition | Number |
| --------- | ------ |
| Naive | 1,679 |
| IFNg | 1,463 |
| Salmonella | 1,202 | 
| IFNg + Salmonella | 860 | 

---

The gene symbols of genes associated with the eQTLs were than obtained based on the Ensembl gene IDs provided in the eQTL catalogue data and save as text files under: `/well/jknight/shiyao/data/alasoo/`.

```python
# Import modules
import pandas as pd

# Load recomputed eQTLs
filtered = pd.read_csv('/well/jknight/shiyao/data/alasoo/overlap_all_filtered.csv')
filtered.sort_values('pvalue', inplace=True)

# Save list of gene symbols to text file
filtered = pd.read_csv('/well/jknight/shiyao/data/alasoo/overlap_all_filtered.csv')
filtered.sort_values('pvalue', inplace=True)
ensembl_info = pd.read_csv('/well/jknight/shiyao/data/alasoo/hg38.ensGene.txt', usecols=['gene_id', 'symbol'])

# All conditions
genes = filtered[['gene_id', 'pvalue']].copy()
genes.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol = genes.merge(ensembl_info, how='left', on=['gene_id'])
symbol[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/alasoo/genes.txt', index=False, header=False, sep='\t')

# IFNg+Salmonella
ifnsal = filtered.loc[filtered.Condition == 'IFNg+Salmonella'][['gene_id', 'pvalue']].copy()
ifnsal.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_ifnsal = ifnsal.merge(ensembl_info, how='left', on=['gene_id'])
symbol_ifnsal[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/alasoo/genes_IfnSal.txt', index=False, header=False, sep='\t')

# IFNg
ifn = filtered.loc[filtered.Condition == 'IFNg'][['gene_id', 'pvalue']].copy()
ifn.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_ifn = ifn.merge(ensembl_info, how='left', on=['gene_id'])
symbol_ifn[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/alasoo/genes_IFN.txt', index=False, header=False, sep='\t')

# Salmonella
sal = filtered.loc[filtered.Condition == 'Salmonella'][['gene_id', 'pvalue']].copy()
sal.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_sal = sal.merge(ensembl_info, how='left', on=['gene_id'])
symbol_sal[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/alasoo/genes_Salmonella.txt', index=False, header=False, sep='\t')

# Naive
naive = filtered.loc[filtered.Condition == 'Naive'][['gene_id', 'pvalue']].copy()
naive.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_naive = naive.merge(ensembl_info, how='left', on=['gene_id'])
symbol_naive[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/alasoo/genes_naive.txt', index=False, header=False, sep='\t')
```
