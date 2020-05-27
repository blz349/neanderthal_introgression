# Neanderthal code

### Requirements:
Python 3

Python modules:
* pandas 1.0.3
* seaborn 0.10.0
* Matplotlib 3.1.1
* natsort 7.0.1
* scikit-allel 1.2.1
* SciPy 1.4.1

---

### Description of code
    compare_neanderthal_datasets.ipynb
This code investigates the overlap in Neanderthal-introgressed SNPs from 2 sources: Dannemann *et al.* (2017) and Simonti *et al.* (2016). Plots include: venn diagram showing overlap of Dannemann and Simonti datasets, and chromosomal distribution of Neanderthal-introgressed SNPs by source. A new csv file containing the chromosome, position and source of the Neanderthal-introgressed SNPs from both datasets was created and saved under: /well/jknight/shiyao/data/comparison_df.csv  
&nbsp;

    neanderthal_superpopulation_mafs.ipynb
This code extracts the MAFs for Neanderthal-introgressed SNPs in the 1000 genomes super populations dataset, which comprises a summarised list of 5026 SNPs for 5 super populations: African (AFR), Admixed American (AMR),  East Asian (EAS), European (EUR), South Asian (SAS). Violin plots were created to visualise the range of MAFs for the Neanderthal-introgressed SNPs in the super populations. Pairplots were created to compare the MAFs in African vs European population and correlation tests were performed. A new csv file containing the following information for the Neanderthal-introgressed SNPs: chromosome, position, source, ID, major allele, minor allele and MAFs in each of the 5 superpopulations was created, and saved under: /well/jknight/shiyao/data/superpop_df.csv  
&nbsp;

    get_neanderthal_allpopulation_mafs.md
This documents describes extracting MAFs for Neanderthal-introgressed SNPs from the 1000 genomes all populations dataset which comprises 26 populations. MAFs were first extracted for each population and saved in .csv format under: /well/jknight/shiyao/data/pop_dfs. The .csv files for all 26 populations were then combined into a single .csv file containing the following information for the Neanderthal-introgressed SNPs: chromosome, position, source, ID, major allele, minor allele and MAFs in each of the 26 populations. This is saved under /well/jknight/shiyao/data/allpop_df.csv  
&nbsp;

    neanderthal_allpopulation_mafs.ipynb
This code plots the range of MAFs for the Neanderthal-introgressed SNPs in the 26 1000 genomes populations, and in particular in the 5 African populations. The number of Neanderthal-introgressed SNPs with MAF > 0.01 in African populations was calculated and visualised using venn diagrams.

---

### Datasets
Neanderthal-introgressed SNPs from:
1. Dannemann M, Prufer K & Kelso J. Functional implications of Neandertal introgression in modern humans. Genome Biol 2017 18:61.
2. Simonti CN *et al.* The phenotypic legacy of admixture between modern humans and Neandertals. Science 2016 351:737-41.  
&nbsp;

1000 genomes super populations data and all populations data from:
* [1000 genomes](https://www.internationalgenome.org/data/)
