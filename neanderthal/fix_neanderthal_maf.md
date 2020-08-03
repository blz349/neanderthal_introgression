# Fix 1000 Genomes minor allele annotation for Neanderthal-introgressed SNPs

The minor allele of SNPs in the 1000 Genomes Project data is sometimes incorrectly annotate, resulting in flipped minor allele frequencies (MAFs). This document contains code for reannotating the minor allele of the Neanderthal-introgressed SNPs using the allele counts provided in the 1000 Genomes data, and correcting the MAFs accordingly.

---

The alleles, allele counts and MAFs of each SNP for each of the 26 1000 Genomes populations were saved in text files under: `/well/jknight/shiyao/data/1000genome/all_populations/`

First, the total allele counts were calculated from the individual population allele counts and saved as a csv file:
```python
# Import modules
import pandas as pd
from glob import glob

# Get list of populations
populations = []
files = glob('/well/jknight/shiyao/data/1000genome/all_populations/*')
for file_name in files:
    temp = file_name.split('/')
    populations.append(temp[-1])

# Get population allele counts
for i in range(6, 28, 7):			# Split into 4 dataframes due to limited memory
	pop_dfs = []
	for pop in populations[i-7:i]:
		if pop != 'CHS':
			dfs = []
			for i in range(1, 23):
				df = pd.read_csv('/well/jknight/shiyao/data/1000genome/all_populations/' + pop + '/' + pop + '_chr' + str(i) + '_freq_genome.txt', sep='\t', header=None, usecols=[1, 4, 5])
				dfs.append(df)
			concat = pd.concat(dfs)
			concat.rename(columns={1: 'ID', 4: pop + '_maj', 5: pop + '_min'}, inplace=True)
			pop_dfs.append(concat)
	pop_df = pd.concat([df.set_index('ID') for df in pop_dfs], axis=1, join='inner').reset_index()
	pop_df.to_csv('/well/jknight/shiyao/data/1000genome/pop_count' + str((i+1)/7) + '.csv', index=False)

# Sum allele counts for each df
dfs = []
for i in range(1, 5):
    df = pd.read_csv('/well/jknight/shiyao/data/1000genome/allele_count/pop_count' + str(i) + '.csv')
    maj_col = [x for x in list(df) if 'maj' in x]
    df['df' + str(i) + '_maj'] = df[maj_col].sum(axis=1)
    min_col = [x for x in list(df) if 'min' in x]
    df['df' + str(i) + '_min'] = df[min_col].sum(axis=1)
    pop_sum = df[['ID', 'df' + str(i) + '_maj', 'df' + str(i) + '_min']].copy()
    dfs.append(pop_sum)
pop_df = pd.concat([df.set_index('ID') for df in dfs], axis=1, join='inner').reset_index()

# Get total allele counts
maj_col = [x for x in list(pop_df) if 'maj' in x]
pop_df['Major count'] = pop_df[maj_col].sum(axis=1)
min_col = [x for x in list(pop_df) if 'min' in x]
pop_df['Minor count'] = pop_df[min_col].sum(axis=1)
pop_df.to_csv('/well/jknight/shiyao/data/1000genome/allele_count.csv', index=False)
```

The total allele counts across all populations were then used to fix the minor allele annotation and MAF:
```python
import pandas as pd

# Load total allele counts and neanderthal SNPs
count = pd.read_csv('/well/jknight/shiyao/data/1000genome/allele_count.csv', usecols=['ID', 'Major count', 'Minor count'])
allpop = pd.read_csv('/well/jknight/shiyao/data/allpop_df.csv')
allpop.drop('Unnamed: 0', axis=1, inplace=True)
fixed = allpop.merge(count, how='left', on='ID')

# Fix minor allele annotation
fixed.loc[fixed['Major count'] < fixed['Minor count'], ['Major', 'Minor']] = fixed.loc[fixed['Major count'] < fixed['Minor count'], ['Minor', 'Major']].values

# Fix population MAFs
cols = ['ACB', 'PUR', 'IBS', 'EUR', 'PJL', 'CHB', 'ESN', 'JPT', 'GBR', 'CEU', 'MXL', 'BEB', 'GWD', 'KHV', 'ASW', 'TSI',
        'LWK', 'YRI', 'GIH', 'FIN', 'PEL', 'ITU', 'MSL', 'STU', 'CLM', 'CDX']
for col in cols:
    fixed.loc[fixed['Major count'] < fixed['Minor count'], [col]] = 1 - fixed.loc[fixed['Major count'] < fixed['Minor count'], [col]].values

# Save Neanderthal SNPs with fixed population MAFs
fixed.to_csv('/well/jknight/shiyao/data/allpop_fixed.csv', index=False)
```

---

The summary stats before and after fixing MAFs are shown in Table 1 and Table 2 respectively.

**Table 1: Before fixing MAF**
|  | AFR | AFR admixed | AMR | EUR | SAS | EAS |
| -- | -- | -- | -- | -- | -- | -- |
| mean  | 0.020964 | 0.026115 | 0.051344 | 0.051390 | 0.051950 | 0.060114 |
| std  | 0.135506 | 0.132277 | 0.128224 | 0.126899 | 0.127333 | 0.138183 |
| min | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 25% | 0.00 | 0.00 | 0.002 | 0.002 | 0.0054 |0.00125 |
| 50% | 0.00 | 0.004 | 0.015 | 0.0144 | 0.0204 | 0.0135 |
| 75% | 0.001 | 0.0105 | 0.0485 | 0.051 | 0.0514 | 0.06275 |
| max | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 

**Table 2: After fixing MAF**
|  | AFR | AFR admixed | AMR | EUR | SAS | EAS |
| -- | -- | -- | -- | -- | -- | -- |
| mean  | 0.002280 | 0.008209 | 0.036729 | 0.037388 | 0.036630 | 0.044428 |
| std  | 0.011045 | 0.014684 | 0.055835 | 0.057866 | 0.047496 | 0.071073 |
| min | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 25% | 0.00 | 0.00 | 0.002 | 0.002 | 0.0054 |0.00125 |
| 50% | 0.00 | 0.004 | 0.01475 | 0.0142 | 0.0202 | 0.0125 | 
| 75% | 0.001 | 0.0105 | 0.04725 | 0.0498 | 0.0498 | 0.05875 | 
| max | 0.9972 | 0.9 | 0.887 | 0.8726 | 0.8204 | 0.91575 |
