# Explore overlap between Neanderthal-introgressed SNPs and recomputed Fairfax eQTLs from EMBL-EBI eQTL catalogue

This document contains code to explore if any Neanderthal-introgressed SNPs are present in the list of recomputed eQTLs for Fairfax *et al.* (2014) from the EMBL-EBI eQTL catalogue. Overlapping recomputed eQTLs were filtered for p-value < 10<sup>-8</sup> and saved in CSV format. Gene symbols of genes associated with these eQTLs were obtained and saved to text files.

First, recomputed eQTLs for Fairfax *et al.* (2014) were downloaded from the EMBL-EBI eQTL catalogue [ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Fairfax_2014/microarray/] and saved under: `/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/`.  

As the EMBL-EBI eQTL catalogue uses GRCh38 coordinates, coordinates of Neanderthal-introgressed SNPs were first converted from GRCh37 to GRCh38 and saved in CSV format under: `/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/comparison_df_hg38.csv`. This requires the `hg19ToHg38.over.chain.gz` file downloaded from [USBC](https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/). 

```python
# Import modules
import pandas as pd
from pyliftover import LiftOver

# Load Neanderthal-introgressed SNPs
neanderthal = pd.read_csv('/well/jknight/shiyao/data/comparison_df.csv')
neanderthal.drop('Unnamed: 0', axis=1, inplace=True)

# Liftover from hg19 to hg38
lo = LiftOver('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/hg19ToHg38.over.chain.gz')
chr_hg38 = []
pos_hg38 = []
for chr, pos in zip(neanderthal['Chromosome'], neanderthal['Position']):
    temp = lo.convert_coordinate('chr' + str(chr), int(pos))
    if len(temp) < 1:
        chr_hg38.append(float('nan'))
        pos_hg38.append(float('nan'))
    else:
        new_chr = temp[0][0]
        if '_' in new_chr:
            new_chr = new_chr.split('_', 1)[0].lstrip('chr')
        else:
            new_chr = new_chr.lstrip('chr')
        new_pos = temp[0][1]
        chr_hg38.append(int(new_chr))
        pos_hg38.append(int(new_pos))
neanderthal['chr_hg38'] = chr_hg38
neanderthal['pos_hg38'] = pos_hg38

# Save to csv file
neanderthal.to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/comparison_df_hg38.csv', index=False)
```

---

Next, the recomputed Fairfax eQTLs and Neanderthal-introgressed SNPs were compared to look for overlap. Overlapping eQTLs were filtered for p-value < 10<sup>-8</sup> and saved in CSV format under: `/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_all_filtered.csv`. The number of eQTLs expressed per condition was computed and shown in Table 1.

```python
# Import modules
import pandas as pd

# Load Neanderthal SNPs with hg38 coordinates
neanderthal = pd.read_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/comparison_df_hg38.csv')

# Overlap between recomputed Fairfax eQTLs and Neanderthal SNPs
# IFN
IFN = pd.read_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/Fairfax_2014_IFN24.all.tsv.gz', sep='\t', compression='gzip')
IFN.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_IFN = neanderthal.merge(IFN, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_IFN.to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_IFN.csv', index=False)
overlap_IFN['rsid'].nunique()

# LPS2
LPS2 = pd.read_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/Fairfax_2014_LPS2.all.tsv.gz', sep='\t', compression='gzip')
LPS2.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_LPS2 = neanderthal.merge(LPS2, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_LPS2.to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_LPS2.csv', index=False)
overlap_LPS2['rsid'].nunique()

# LPS24
LPS24 = pd.read_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/Fairfax_2014_LPS24.all.tsv.gz', sep='\t', compression='gzip')
LPS24.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_LPS24 = neanderthal.merge(LPS24, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_LPS24.to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_LPS24.csv', index=False)
overlap_LPS24['rsid'].nunique() 

# Naive
Naive = pd.read_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/Fairfax_2014_naive.all.tsv.gz', sep='\t', compression='gzip')
Naive.rename(columns={'chromosome': 'chr_hg38', 'position': 'pos_hg38'}, inplace=True)
overlap_Naive = neanderthal.merge(Naive, how='inner', on=['chr_hg38', 'pos_hg38'])
overlap_Naive.to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_naive.csv', index=False)
overlap_Naive['rsid'].nunique() 

# Concatenate dfs
IFN = pd.read_csv('overlap_IFN.csv')
IFN['Condition'] = 'IFN'
LPS2 = pd.read_csv('overlap_LPS2.csv')
LPS2['Condition'] = 'LPS2'
LPS24 = pd.read_csv('overlap_LPS24.csv')
LPS24['Condition'] = 'LPS24'
Naive = pd.read_csv('overlap_naive.csv')
Naive['Condition'] = 'Naive'
combined = pd.concat([IFN, LPS2, LPS24, Naive])
combined.to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_all.csv', index=False)
combined['rsid'].nunique()  # 184528

# Keep only eQTLs with p < 10e-8
filtered = combined.loc[combined['pvalue'] < 10e-8].copy()
filtered.to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_all_filtered.csv', index=False)
filtered['rsid'].nunique()  # 12047

# Get number of eQTLs per condition
condition = filtered[['rsid', 'Condition']]
condition = temp.groupby(['rsid']).agg(lambda col: ','.join(col))
for i, row in condition.iterrows():
    condition.at[i, 'Condition'] = ','.join(set(condition.at[i, 'Condition'].split(',')))
condition.Condition.str.split(',', expand=True).stack().value_counts()
```

12,047 recomputed eQTLs (p < 10<sup>-8</sup>) were also identified as Neanderthal-introgressed SNPs.

**Table 1:** Number of condition-specific eQTLs per condition
| Condition | Number |
| --------- | ------ |
| Naive | 8,282 |
| IFN | 7,018 |
| LPS24 | 5,620 | 
| LPS2 | 4,212 | 

---

To minimise overlapping signals, for eQTL-gene pairs expressed under multiple conditions, only the condition for which the eQTL has the lowest p-value was kept. The gene symbols of genes associated with the eQTLs were than obtained based on the Ensembl gene IDs provided in the eQTL catalogue data and save as text files under: `/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/`.

```python
# Import modules
import pandas as pd

# Clean up eQTLs expressed in multiple conditions
filtered = pd.read_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_all_filtered.csv')
filtered.sort_values('pvalue', inplace=True)
cleaned = filtered.drop_duplicates(subset=['rsid', 'gene_id'], keep='first').copy()

# Save list of gene symbols to text file
ensembl_info = pd.read_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/hg38.ensGene.txt', usecols=['gene_id', 'symbol'])

# IFN
ifn_c = cleaned.loc[cleaned.Condition == 'IFN'][['rsid', 'gene_id', 'pvalue']].copy()
ifn_c.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_ifn_c = ifn_c.merge(ensembl_info, how='left', on=['gene_id'])
symbol_ifn_c[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/genes_cleaned_ifn.txt', index=False, header=False, sep='\t')

# LPS2
lps2_c = cleaned.loc[cleaned.Condition == 'LPS2'][['rsid', 'gene_id', 'pvalue']].copy()
lps2_c.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_lps2_c = lps2_c.merge(ensembl_info, how='left', on=['gene_id'])
symbol_lps2_c[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/genes_cleaned_lps2.txt', index=False, header=False, sep='\t')

# LPS24
lps24_c = cleaned.loc[cleaned.Condition == 'LPS24'][['rsid', 'gene_id', 'pvalue']].copy()
lps24_c.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_lps24_c = lps24_c.merge(ensembl_info, how='left', on=['gene_id'])
symbol_lps24_c[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/genes_cleaned_lps24.txt', index=False, header=False, sep='\t')

# Naive
naive_c = cleaned.loc[cleaned.Condition == 'Naive'][['rsid', 'gene_id', 'pvalue']].copy()
naive_c.drop_duplicates(subset='gene_id', keep='first', inplace=True)
symbol_naive_c = naive_c.merge(ensembl_info, how='left', on=['gene_id'])
symbol_naive_c[['symbol', 'pvalue']].to_csv('/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/genes_cleaned_naive.txt', index=False, header=False, sep='\t')
```
