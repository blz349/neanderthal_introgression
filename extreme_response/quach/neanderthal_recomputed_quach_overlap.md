# Explore overlap between Neanderthal-introgressed SNPs and recomputed Quach eQTLs from EMBL-EBI eQTL catalogue

The Quach *et al.* (2016) study contains monocyte eQTLs for the 5 conditions: Pam3CSK4 (TLR1/2 agonist), R848 (TLR 7/8 agonist), LPS, Influenza A virus (IAV) and Naive. This document contains code to explore if any Neanderthal-introgressed SNPs are present in the list of recomputed eQTLs for Quach *et al.* (2016) from the EMBL-EBI eQTL catalogue. Overlapping recomputed eQTLs were filtered for p-value < 10<sup>-8</sup> and saved in CSV format. Gene symbols of genes associated with these eQTLs were obtained and saved to text files.

First, recomputed eQTLs for Quach *et al.* (2016) were downloaded from the EMBL-EBI eQTL catalogue [ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Quach_2016/ge/] and saved under: `/well/jknight/shiyao/data/quach/`.  

The recomputed Fairfax eQTLs and Neanderthal-introgressed SNPs were compared to look for overlap. Overlapping eQTLs were filtered for p-value < 10<sup>-8</sup> and saved in CSV format under: `/well/jknight/shiyao/data/quach/overlap_all_filtered.csv`. The number of eQTLs expressed per condition was computed and shown in Table 1.

```python
# Import modules
import pandas as pd

# Load Neanderthal SNPs with hg38 coordinates
neanderthal = pd.read_csv('/well/jknight/shiyao/data/quach/comparison_df_hg38.csv')

# Overlap between recomputed Quach eQTLs and Neanderthal SNPs
# IAV
IAV = pd.read_csv('/well/jknight/shiyao/data/quach/Quach_2016_ge_monocyte_IAV.all.tsv.gz', sep='\t', compression='gzip')
IAV.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_IAV = neanderthal.merge(IAV, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_IAV.to_csv('/well/jknight/shiyao/data/quach/overlap_IAV.csv', index=False)

# LPS
LPS = pd.read_csv('/well/jknight/shiyao/data/quach/Quach_2016_ge_monocyte_LPS.all.tsv.gz', sep='\t', compression='gzip')
LPS.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_LPS = neanderthal.merge(LPS, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_LPS.to_csv('/well/jknight/shiyao/data/quach/overlap_LPS.csv', index=False)

# Pam3CSK4
Pam = pd.read_csv('/well/jknight/shiyao/data/quach/Quach_2016_ge_monocyte_Pam3CSK4.all.tsv.gz', sep='\t', compression='gzip')
Pam.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_Pam = neanderthal.merge(Pam, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_Pam.to_csv('/well/jknight/shiyao/data/quach/overlap_Pam3CSK4.csv', index=False)

# R848
R848 = pd.read_csv('/well/jknight/shiyao/data/quach/Quach_2016_ge_monocyte_R848.all.tsv.gz', sep='\t', compression='gzip')
R848.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_R848 = neanderthal.merge(R848, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_R848.to_csv('/well/jknight/shiyao/data/quach/overlap_R848.csv', index=False)

# Naive
Naive = pd.read_csv('/well/jknight/shiyao/data/quach/Quach_2016_ge_monocyte_naive.all.tsv.gz', sep='\t', compression='gzip')
Naive.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_Naive = neanderthal.merge(Naive, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_Naive.to_csv('/well/jknight/shiyao/data/quach/overlap_naive.csv', index=False)

# Concatenate dfs
IAV = pd.read_csv('overlap_IAV.csv')
IAV['Condition'] = 'IAV'
LPS = pd.read_csv('overlap_LPS.csv')
LPS['Condition'] = 'LPS'
Pam = pd.read_csv('overlap_Pam3CSK4.csv')
Pam['Condition'] = 'Pam3CSK4'
R848 = pd.read_csv('overlap_R848.csv')
R848['Condition'] = 'R848'
Naive = pd.read_csv('overlap_naive.csv')
Naive['Condition'] = 'Naive'
combined = pd.concat([IAV, LPS, Pam, R848, Naive])
combined.to_csv('/well/jknight/shiyao/data/quach/overlap_all.csv', index=False)

# Keep only eQTLs with p < 10e-8
filtered = combined.loc[combined['pvalue'] < 10e-8].copy()
filtered.to_csv('/well/jknight/shiyao/data/quach/overlap_all_filtered.csv', index=False)
filtered['rsid'].nunique()

# Get number of eQTLs per condition
condition = filtered[['rsid', 'Condition']]
condition = condition.groupby(['rsid']).agg(lambda col: ','.join(col))
for i, row in condition.iterrows():
    condition.at[i, 'Condition'] = ','.join(set(condition.at[i, 'Condition'].split(',')))
condition.Condition.str.split(',', expand=True).stack().value_counts()
```

3,192 recomputed eQTLs (p < 10<sup>-8</sup>) were also identified as Neanderthal-introgressed SNPs.

**Table 1:** Number of condition-specific eQTLs per condition
| Condition | Number |
| --------- | ------ |
| Pam3CSK4 | 1,728 |
| R848 | 1,693 |
| LPS | 1,613 | 
| Naive | 1,543 | 
| IAV | 1,318 |

---

The gene symbols of genes associated with the eQTLs were than obtained based on the Ensembl gene IDs provided in the eQTL catalogue data and save as text files under: `/well/jknight/shiyao/data/quach/`.

```python
# Import modules
import pandas as pd

# Load recomputed eQTLs
filtered = pd.read_csv('/well/jknight/shiyao/data/quach/overlap_all_filtered.csv')
filtered.sort_values('pvalue', inplace=True)

# Save list of gene symbols to text file
filtered = pd.read_csv('/well/jknight/shiyao/data/quach/overlap_all_filtered.csv')
filtered.sort_values('pvalue', inplace=True)
ensembl_info = pd.read_csv('/well/jknight/shiyao/data/quach/hg38.ensGene.txt', usecols=['gene_id', 'symbol'])

# All conditions
genes = filtered[['gene_id', 'pvalue']].copy()
genes.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol = genes.merge(ensembl_info, how='left', on=['gene_id'])
symbol[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/quach/genes.txt', index=False, header=False, sep='\t')

# IAV
iav = filtered.loc[filtered.Condition == 'IAV'][['gene_id', 'pvalue']].copy()
iav.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_iav = iav.merge(ensembl_info, how='left', on=['gene_id'])
symbol_iav[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/quach/genes_IAV.txt', index=False, header=False, sep='\t')

# LPS
lps = filtered.loc[filtered.Condition == 'LPS'][['gene_id', 'pvalue']].copy()
lps.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_lps = lps.merge(ensembl_info, how='left', on=['gene_id'])
symbol_lps[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/quach/genes_LPS.txt', index=False, header=False, sep='\t')

# Pam3CSK4
pam = filtered.loc[filtered.Condition == 'Pam3CSK4'][['gene_id', 'pvalue']].copy()
pam.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_pam = pam.merge(ensembl_info, how='left', on=['gene_id'])
symbol_pam[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/quach/genes_Pam3CSK4.txt', index=False, header=False, sep='\t')

# R848
r848 = filtered.loc[filtered.Condition == 'R848'][['gene_id', 'pvalue']].copy()
r848.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_r848 = r848.merge(ensembl_info, how='left', on=['gene_id'])
symbol_r848[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/quach/genes_R848.txt', index=False, header=False, sep='\t')

# Naive
naive = filtered.loc[filtered.Condition == 'Naive'][['gene_id', 'pvalue']].copy()
naive.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_naive = naive.merge(ensembl_info, how='left', on=['gene_id'])
symbol_naive[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/quach/genes_naive.txt', index=False, header=False, sep='\t')
```
