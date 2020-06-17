# Explore overlap between Neanderthal-introgressed SNPs and recomputed Nedelec eQTLs from EMBL-EBI eQTL Catalogue

This document contains code to explore if any Neanderthal-introgressed SNPs are present in the list of recomputed eQTLs for Nedelec *et al.* (2016) from the EMBL-EBI eQTL catalogue. Overlapping recomputed eQTLs were filtered for p-value < 10<sup>-8</sup> and saved in CSV format. Gene symbols of genes associated with these eQTLs were obtained and saved to text files.

First, recomputed eQTLs for Nedelec *et al.* (2016) were downloaded from the EMBL-EBI eQTL catalogue [ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Nedelec_2016/ge/] and saved under: `/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/`.  

The recomputed Nedelec eQTLs and Neanderthal-introgressed SNPs were compared to look for overlap. Overlapping eQTLs were filtered for p-value < 10<sup>-8</sup> and saved in CSV format under: `/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/overlap_all_filtered.csv`. The number of eQTLs expressed per condition was computed and shown in Table 1.

```python
# Import modules
import pandas as pd

# Load Neanderthal SNPs with hg38 coordinates
neanderthal = pd.read_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/comparison_df_hg38.csv')
# Overlap between recomputed Nedelec eQTLs and Neanderthal SNPs

# Listeria
lis = pd.read_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/Nedelec_2016_ge_macrophage_Listeria.all.tsv.gz', sep='\t', compression='gzip')
lis.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_lis = neanderthal.merge(lis, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_lis.to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/overlap_listeria.csv', index=False)
overlap_lis['rsid'].nunique()  # 148527

# Salmonella
sal = pd.read_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/Nedelec_2016_ge_macrophage_Salmonella.all.tsv.gz', sep='\t', compression='gzip')
sal.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_sal = neanderthal.merge(sal, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_sal.to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/overlap_salmonella.csv', index=False)
overlap_sal['rsid'].nunique()  # 149414

# Naive
naive = pd.read_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/Nedelec_2016_ge_macrophage_naive.all.tsv.gz', sep='\t', compression='gzip')
naive.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_naive = neanderthal.merge(naive, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_naive.to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/overlap_naive.csv', index=False)
overlap_naive['rsid'].nunique()  # 148212

# Concatenate dfs
lis = pd.read_csv('overlap_listeria.csv')
lis['Condition'] = 'Listeria'
sal = pd.read_csv('overlap_salmonella.csv')
sal['Condition'] = 'Salmonella'
naive = pd.read_csv('overlap_naive.csv')
naive['Condition'] = 'Naive'
combined = pd.concat([lis, sal, naive])
combined.to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/overlap_all.csv', index=False)
combined['rsid'].nunique()  # 150199

# Keep only eQTLs with p < 10e-8
filtered = combined.loc[combined['pvalue'] < 10e-8].copy()
filtered.to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/overlap_all_filtered.csv', index=False)
filtered['rsid'].nunique()  # 1924

# Get number of eQTLs per condition
condition = filtered[['rsid', 'Condition']]
condition = condition.groupby(['rsid']).agg(lambda col: ','.join(col))
for i, row in condition.iterrows():
    condition.at[i, 'Condition'] = ','.join(set(condition.at[i, 'Condition'].split(',')))
condition.Condition.str.split(',', expand=True).stack().value_counts()
```

1,924 recomputed eQTLs (p < 10<sup>-8</sup>) were also identified as Neanderthal-introgressed SNPs.

**Table 1:** Number of condition-specific eQTLs per condition
| Condition | Number |
| --------- | ------ |
| Listeria-infected | 1,052 |
| Salmonella-infected | 979 |
| Naive | 954 | 

---

The gene symbols of genes associated with the eQTLs were than obtained based on the Ensembl gene IDs provided in the eQTL catalogue data and save as text files under: `/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/` for pathway analysis.

```python
# Import modules
import pandas as pd

# Load recomputed eQTLs
filtered = pd.read_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/overlap_all_filtered.csv')
filtered.sort_values('pvalue', inplace=True)

# Save list of gene symbols to text file
ensembl_info = pd.read_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/hg38.ensGene.txt', usecols=['gene_id', 'symbol'])

# All conditions
genes = filtered[['gene_id', 'pvalue']].copy()
genes.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol = genes.merge(ensembl_info, how='left', on=['gene_id'])
symbol[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/genes.txt', index=False, header=False, sep='\t')

# Listeria
LIS = filtered.loc[filtered.Condition == 'Listeria'][['gene_id', 'pvalue']].copy()
LIS.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_lis = LIS.merge(ensembl_info, how='left', on=['gene_id'])
symbol_lis[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/genes_listeria.txt', index=False, header=False, sep='\t')

# Salmonella
SAL = filtered.loc[filtered.Condition == 'Salmonella'][['gene_id', 'pvalue']].copy()
SAL.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_sal = SAL.merge(ensembl_info, how='left', on=['gene_id'])
symbol_sal[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/genes_salmonella.txt', index=False, header=False, sep='\t')

# Naive
Naive = filtered.loc[filtered.Condition == 'Naive'][['gene_id', 'pvalue']].copy()
Naive.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_naive = Naive.merge(ensembl_info, how='left', on=['gene_id'])
symbol_naive[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/genes_naive.txt', index=False, header=False, sep='\t')
```
