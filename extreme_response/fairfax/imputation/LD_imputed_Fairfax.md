# LD values for imputed SNPs with Fairfax eQTLs

This documents describes extracting linkage disequilibrium (LD) values for imputed SNPs that crossover with Neanderthal-introgressed SNPs versus Fairfax eQTLs. The r<sup>2</sup> statistic was used as a measure of LD and was obtained for European populations using the 1000 genomes dataset.

First, the 1000 Genomes Phase 3 VCF files for each chromosome were downloaded from [ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/] and saved under: `/well/jknight/shiyao/data/1000genome/vcfs/`.

The following script which is saved under: `/well/jknight/shiyao/data/1000genome/vcfs/concat_vcf.sh` was submitted to the cluster to concatenate the VCF files for chromosomes 1 to 22 into a single VCF file.

```bash
#!/bin/bash

#$ -cwd -V
#$ -N concat_vcf
#$ -P jknight.prjc
#$ -q short.qc
#$ -o /well/jknight/shiyao/data/1000genome/vcfs/concat_vcf_stdout.log
#$ -e /well/jknight/shiyao/data/1000genome/vcfs/concat_vcf_stderr.log
#$ -pe shmem 2

echo "************************"
echo "Job ID: "$JOB_ID
echo "Task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "************************"

/apps/well/bcftools/1.3/bin/bcftools concat ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -Oz -o  ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

echo "************************"
echo "Finished at: "`date`
echo "************************"
```

This VCF file was then converted into PLINK format (bed/bim) via the command below.

```bash
/apps/well/plink/1.90b3/plink --vcf ALL.autosomes.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out ../all_autosomes
```

---

European sample IDs was obtained from the sample panel and used to reduce the data to European individuals. The PLINK file was saved under: `/well/jknight/shiyao/data/1000genome/merged`

```bash
grep -e 'EUR' integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1, $1}' > ../eur_ids.txt

/apps/well/plink/1.90b3/plink --bfile merged/all_autosomes --keep eur_ids.txt --make-bed --out merged/all_autosomes_eur
```

Next, the variant IDs of imputed that crossover with Neanderthal-introgressed SNPs were obtained and saved under: `/well/jknight/shiyao/data/fairfax/imputed_LD/imputed_intersected_ids.txt`
```python
# Import modules
import pandas as pd

# Imputed SNPs (after removing high LD regions)
imputed_noLD = pd.read_csv('/well/jknight/shiyao/data/fairfax/imputation/info_0.9/fairfax_eqtl_genotyping_NoLongRangeLD.bim', sep='\t', header=None, usecols=[0, 3])
imputed_noLD.rename(columns={0: "Chromosome", 3: "Position"}, inplace=True)

# Neanderthal-introgressed SNPs
neanderthal = pd.read_csv('/well/jknight/shiyao/data/allpop_df.csv', usecols=['Chromosome', 'Position', 'ID'])

# Intersecting SNPs
intersect = neanderthal.merge(imputed_noLD, how="inner", on=["Chromosome", "Position"])
intersect.dropna(inplace=True)
intersect['ID'].to_csv('/well/jknight/shiyao/data/fairfax/imputed_LD/imputed_intersected_ids.txt',
                       index=False, header=False)
```

The data was then further reduced to only include SNPs that are present in both imputed and neanderthal datasets and saved under: `/well/jknight/shiyao/data/fairfax/LD_imputed/`
 
```bash
/apps/well/plink/1.90b3/plink --bfile /well/jknight/shiyao/data/1000genome/merged/all_autosomes_eur --extract imputed_intersected_ids.txt --make-bed --out all_autosomes_eur_intersected
```

---

The list of Fairfax eQTL variant IDs was obtained using the following command:

```bash
awk 'NR > 1 {print $1}' tab2_a_cis_eSNPs.txt | sort | uniq > LD_imputed/fairfax_eQTLs.txt
```

The LD scores for intersecting SNPs with Fairfax eQTLs was obtained using PLINK and saved under: `/well/jknight/shiyao/data/fairfax/LD_imputed/`. r2 values < 0.9 were excluded.

```bash
/apps/well/plink/1.90b3/plink --bfile all_autosomes_eur_intersected --r2 --ld-window 100000 --ld-window-kb 10000000 --ld-window-r2 0.9 --ld-snp-list fairfax_eQTLs.txt --out LD_intersected_fairfax_0.9
```

Finally, the .ld file was converted into a text file and saved under: `/well/jknight/shiyao/data/fairfax/LD_imputed/`
```bash
awk '{print $1, $2, $3, $4, $5, $6, $7}' LD_intersected_fairfax_0.9.ld > LD_intersected_fairfax_0.9.txt
```

---

Python was used to calculate the number of SNPs (out of 76,453 intersecting SNPs) that were in LD with a Fairfax eQTL.

```python
# Import libraries
import pandas as pd

# Load LD data
LD = pd.read_csv('/well/jknight/shiyao/data/fairfax/LD_imputed/LD_intersected_fairfax_0.9.txt', sep=' ')

# Keep highest R2 value for SNPs linked to >1 eQTLs
LD[LD.duplicated("SNP_B", keep=False)]
LD.sort_values("R2")
LD.drop_duplicates(subset="SNP_B", keep="last", inplace=True)

# Remove Fairfax eQTLs
LD = LD.query("BP_A != BP_B")

# Save df to csv file
LD.sort_values(by=["CHR_B", "BP_B"])
LD.to_csv('LD_intersected_fairfax_0.9.csv', index=False)
```

2,195 SNPs were in LD with a Fairfax eQTL. The r<sup>2</sup> values for these SNPs are saved under: `/well/jknight/shiyao/data/fairfax/LD_imputed/LD_intersected_fairfax_0.9.csv`

---

Information for the 92 Fairfax eQTLs in LD with imputed SNPs that crossover with Neanderthal-introgressed SNPs was obtained from the Fairfax dataset. The number of condition-specific eQTLs per condition was calculated (Table 1) and the genes associated with each eQTL was saved into a text file.

```python
# Import libraries
import pandas as pd

# Load LD data
LD = pd.read_csv('/well/jknight/shiyao/data/fairfax/LD_imputed/LD_intersected_fairfax_0.9.csv')

# Get eQTL info
fairfax = pd.read_csv("/well/jknight/shiyao/data/fairfax/tab2_a_cis_eSNPs.txt", sep="\t", usecols=["Gene", "SNP.Chrm", "SNP.pos", "Min.dataset"])
fairfax = fairfax.rename(columns={"SNP.Chrm": "CHR_A", "SNP.pos": "BP_A"})
data = LD.merge(fairfax, how='left', on=['CHR_A', 'BP_A'])

# Number of eQTLs per condition
condition = data[['SNP_A', 'Min.dataset']].copy()
condition.drop_duplicates(keep='first', inplace=True)
condition = condition.groupby(['SNP_A']).agg(lambda col: ','.join(col))
for i, row in condition.iterrows():
    condition.at[i, 'Min.dataset'] = ','.join(set(condition.at[i, 'Min.dataset'].split(',')))
condition['Min.dataset'].str.split(',', expand=True).stack().value_counts()

# Get associated genes
genes = data[['Gene']].copy()
genes.drop_duplicates(keep='first', inplace=True)
genes.Gene.to_csv('/well/jknight/shiyao/data/fairfax/LD_imputed/genes.txt', index=False, header=False)

ifn = data.loc[data['Min.dataset'] == 'IFN'][['Gene']].copy()
ifn.drop_duplicates(keep='first', inplace=True)
ifn.Gene.to_csv('/well/jknight/shiyao/data/fairfax/LD_imputed/genes_ifn.txt', index=False, header=False)

lps2 = data.loc[data['Min.dataset'] == 'LPS2'][['Gene']].copy()
lps2.drop_duplicates(keep='first', inplace=True)
lps2.Gene.to_csv('/well/jknight/shiyao/data/fairfax/LD_imputed/genes_lps2.txt', index=False, header=False)

lps24 = data.loc[data['Min.dataset'] == 'LPS24'][['Gene']].copy()
lps24.drop_duplicates(keep='first', inplace=True)
lps24.Gene.to_csv('/well/jknight/shiyao/data/fairfax/LD_imputed/genes_lps24.txt', index=False, header=False)

naive = data.loc[data['Min.dataset'] == 'Naive'][['Gene']].copy()
naive.drop_duplicates(keep='first', inplace=True)
naive.Gene.to_csv('/well/jknight/shiyao/data/fairfax/LD_imputed/genes_naive.txt', index=False, header=False)
```

**Table 1:** Number of condition-specific eQTLs per condition
| Condition | Number |
| --------- | ------ |
| LPS24 | 31 |
| IFN | 28 |
| LPS2 | 22 |
| Naive | 19 |
