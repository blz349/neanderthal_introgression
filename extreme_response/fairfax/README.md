# Fairfax code

### Requirements:
Python 3 and modules:
* pandas 1.0.3
* seaborn 0.10.0
* Matplotlib 3.1.1
* SciPy 1.4.1
* plinkio 0.9.7
* networkx 2.4
* statsmodels 0.11.1

R and packages:
* XGR
* RCircos
* Tidyverse
* genio

---

### Description of code

    neanderthal_fairfax_overlap.ipynb
This code looks at the overlap between Dannemann/Simonti Neanderthal-introgressed SNPs (present in Chen introgressed sequences) and monocyte extreme response peak eQTLs identified by Fairfax *et al*. Fisher's exact test and bootstrap test were performed to test for enrichment.  
&nbsp;

    plot_fairfax_eQTLs.ipynb
This code plots the expression data for Fairfax eQTLs for which the same variant has also been identified as a Neanderthal-introgressed SNP. eQTL plots were divided by the 4 treatment groups: IFN, LPS 2h, LPS 24h and Naive.  
&nbsp;  

    compute_rs2066807_STAT2_eQTL.ipynb
This code plots the *STAT2* expression for rs2066807 and calculates the p-value and FDR of the association for each of the 4 contexts: IFN, LPS 2h, LPS 24h and Naive.  
&nbsp;  

    imputation
This folder contains the code for imputing Fairfax eQTLs using the Sanger Imputation Service and calculating the overlap with Neanderthal-introgressed SNPs and LD scores for imputed SNPs.  
&nbsp;

    neanderthal_recomputed_fairfax_overlap.md
This document contains code for determining the overlap between Dannemann/Simonti Neanderthal-introgressed SNPs (from Dannemann *et al.* and Simonti *et al.*) and recomputed Fairfax *et al.* eQTLs from the EMBL-EBI eQTL catalogue. Overlapping eQTLs were filtered for p < 10<sup>-8</sup> and saved as CSV format (coordinates in GRCh38) under: `/well/jknight/shiyao/data/fairfax/EMBL_recomputed/overlap_all_filtered.csv/`. Genes associated with the eQTLs were extracted and saved as text files under: `/well/jknight/shiyao/data/fairfax/EMBL_recomputed/genes_cleaned_*.txt`.  
&nbsp;

    plot_recomputed_fairfax_pathways.ipynb
This code plots the pathways in which genes associated with Neanderthal-introgressed eQTLs were enriched in, using pathway enrichment results obtained from the R package XGR, for each of the 4 conditions: IFN, LPS 2h, LPS 24h and Naive.  
&nbsp;

    plot_recomputed_fairfax_networks.ipynb
This code plots the network of genes associated with eQTLs, which have been identified as Neanderthal-introgressed, using subnet results obtained from the R package XGR, for each of the 4 conditions: IFN, LPS 2h, LPS 24h and Naive.  
&nbsp;

    introgressed_HLA_types.ipynb
This code calculates the number of individuals with each HLA type in the Fairfax *et al.* dataset and plots the HLA gene expression for putative Neanderthal-derived HLA types: HLA-A\*02:01, HLA-C\*07:02, HLA-B\*07:02 and HLA-B\*51:01. For HLA-C\*07:02, t-tests were performed to check if number of HLA-C\*07:02 alleles signficantly altered expression of *HLA-C* or *KIR2DL3*, and Fisher's exact test was used to test for enrichment of KIR2DL2/3 in individuals with HLA-C\*07:02.

---

### Datasets
Neanderthal-introgressed SNPs from:
1. Dannemann M, Prufer K & Kelso J. Functional implications of Neandertal introgression in modern humans. *Genome Biol* 2017 **18**:61.
2. Simonti CN *et al.* The phenotypic legacy of admixture between modern humans and Neandertals. *Science* 2016 **351**:737-41.  
&nbsp;

Neanderthal-introgressed sequences by Chen *et al.* from:
* Chen L *et al.* Identifying and interpreting apparent Neanderthal ancestry in African individuals. *Cell* 2020 **180**:677-687. 
&nbsp;

MAFs from:
* [1000 genomes](https://www.internationalgenome.org/data/)  
&nbsp;

Genotyped SNPs and monocyte, gene-expression data from:
* Fairfax BP *et al.* Innate immune activity conditions the effect of regulatory variants upon monocyte gene expression. *Science* 2014 **343**:1246949.  
&nbsp;

Recomputed Fairfax *et al.* eQTLs from:
* [EMBL-EBI eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/)
