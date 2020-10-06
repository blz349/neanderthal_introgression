# Get list of Neanderthal SNPs present in Chen introgressed sequences
This document contains code that identifies if any Neanderthal-derived SNPs identified by Dannemann *et al.* (2017) and Simonti *et al.* (2016) were present in the Neanderthal-introgressed sequences identified by Chen *et al.* (2020). 

The Neanderthal-derived SNPs from Dannemann and Simonti are saved under: `/well/jknight/shiyao/data/comparison_df.csv`  
The Chen introgressed sequences are saved under: `/well/jknight/shiyao/data/chen/Neanderthal_sequence_in_1000_genome.txt`

```python
# Import modules
import pandas as pd
import seaborn as sns

# Load datasets
chen = pd.read_csv("/well/jknight/shiyao/data/chen/Neanderthal_sequence_in_1000_genome.txt", sep=' ')

neanderthal = pd.read_csv("/well/jknight/shiyao/data/comparison_df.csv")
neanderthal.sort_values(['Chromosome', 'Position'], inplace=True)
neanderthal.reset_index(inplace=True)
neanderthal = neanderthal.drop(['index', 'Unnamed: 0'], axis=1)

# Check for Dannemann/Simonti Neanderthal SNPs in Chen introgressed sequences
neanderthal['Chen'] = 'No'
for i in range(1, 23):
    ori_nean = neanderthal.loc[neanderthal.Chromosome == i]
    chen_nean = chen.loc[chen.chr == i]
    for num_range in chen_nean.itertuples():
        matched = ori_nean['Position'].between(num_range.start, num_range.end)
        neanderthal.loc[matched[matched].index, 'Chen'] = 'Yes'
neanderthal.to_csv('/well/jknight/shiyao/data/chen/neanderthal_check.csv', index=False)
neanderthal.Chen.value_counts()
```

A total of 306,851 (87%) Neanderthal-derived SNPs were present in Neanderthal-introgressed sequences identified by Chen *et al.* (2020). These SNPs were saved under: `/well/jknight/shiyao/data/chen/neanderthal_check.csv`
