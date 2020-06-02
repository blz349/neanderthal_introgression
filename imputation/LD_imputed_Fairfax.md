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

European sample IDs was obtained from the sample panel and used to reduce the data to European individuals. The PLINK file was saved under: `/well/jknight/shiyao/data/1000genome/`

```bash
grep -e 'EUR' integrated_call_samples_v3.20130502.ALL.panel | awk '{print $1, $1}' > ../eur_ids.txt

/apps/well/plink/1.90b3/plink --bfile all_autosomes --keep eur_ids.txt --make-bed --out all_autosomes_eur
```

Next, the data was further reduced to only include SNPs that are present in both imputed and neanderthal datasets and saved under: `/well/jknight/shiyao/data/extreme_response/LD_imputed/`
 
```bash
/apps/well/plink/1.90b3/plink --bfile /well/jknight/shiyao/data/1000genome/all_autosomes_eur --extract imputed_intersected_ids.txt --make-bed --out all_autosomes_eur_intersected
```

---

The list of Fairfax eQTL variant IDs was obtained using the following command:

```bash
awk 'NR > 1 {print $1}' tab2_a_cis_eSNPs.txt | sort | uniq > LD_imputed/fairfax_eQTLs.txt
```

The LD scores for intersecting SNPs with Fairfax eQTLs was obtained using PLINK and saved under: `/well/jknight/shiyao/data/extreme_response/LD_imputed/`. r2 values < 0.9 were excluded.

```bash
/apps/well/plink/1.90b3/plink --bfile all_autosomes_eur_intersected --r2 --ld-window 100000 --ld-window-kb 10000000 --ld-window-r2 0.9 --ld-snp-list fairfax_eQTLs.txt --out LD_intersected_fairfax_0.9
```

Finally, the .ld file was converted into a text file and saved under: `/well/jknight/shiyao/data/extreme_response/LD_imputed/`
```bash
awk '{print $1, $2, $3, $4, $5, $6, $7}' LD_intersected_fairfax_0.9.ld > LD_intersected_fairfax_0.9.txt
```

---

Python was used to calculate the number of SNPs (out of 76,453 intersecting SNPs) that were in LD with a Fairfax eQTL.

```python
# Import libraries
import pandas as pd

# Load LD data
LD = pd.read_csv('/well/jknight/shiyao/data/extreme_response/LD_imputed/LD_intersected_fairfax_0.9.txt', sep=' ')

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

2,195 SNPs were in LD with a Fairfax eQTL. The r<sup>2</sup> values for these SNPs are saved under: `/well/jknight/shiyao/data/extreme_response/LD_imputed/LD_intersected_fairfax_0.9.csv`

