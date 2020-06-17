# Fairfax code

### Requirements:
Python 3

Python modules:
* pandas 1.0.3
* seaborn 0.10.0
* Matplotlib 3.1.1
* SciPy 1.4.1
* plinkio 0.9.7

---

### Description of code

    neanderthal_fairfax_overlap.ipynb
This code looks at the overlap in SNPs between Neanderthal-introgressed SNPs (from Dannemann *et al.* and Simonti *et al.*) and monocyte extreme response peak eQTLs identified by Fairfax *et al*. Fisher's exact test and bootstrap test were performed to test for enrichment.  
&nbsp;

    plot_fairfax_eQTLs.ipynb
This code plots the expression data for Fairfax eQTLs for which the same variant has also been identified as a Neanderthal-introgressed SNP. eQTL plots were divided by the 4 treatment groups: IFN, LPS 2h, LPS 24h and Naive.  
&nbsp;  

    imputation
This folder contains the code for imputing Fairfax eQTLs using the Sanger Imputation Service and calculating the overlap with Neanderthal-introgressed SNPs and LD scores for imputed SNPs.  
&nbsp;

    neanderthal_recomputed_fairfax_overlap.md
This document contains code for determining the overlap between Neanderthal-introgressed SNPs (from Dannemann *et al.* and Simonti *et al.*) recomputed Fairfax *et al.* eQTLs from the EMBL-EBI eQTL catalogue. Overlapping eQTLs were filtered for p < 10<sup>-8</sup> and saved as CSV format (coordinates in GRCh38) under: `/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/overlap_all_filtered.csv/`. Genes associated with the eQTLs were extracted and saved as text files under: `/well/jknight/shiyao/data/extreme_response/EMBL_recomputed/genes_cleaned_*.txt`.  

---

### Datasets
Neanderthal-introgressed SNPs from:
1. Dannemann M, Prufer K & Kelso J. Functional implications of Neandertal introgression in modern humans. *Genome Biol* 2017 **18**:61.
2. Simonti CN *et al.* The phenotypic legacy of admixture between modern humans and Neandertals. *Science* 2016 **351**:737-41.  
&nbsp;

MAFs from:
* [1000 genomes](https://www.internationalgenome.org/data/)  
&nbsp;

Genotyped SNPs and monocyte, gene-expression data from:
* Fairfax BP *et al.* Innate immune activity conditions the effect of regulatory variants upon monocyte gene expression. *Science* 2014 **343**:1246949.  
&nbsp;

Recomputed Fairfax *et al.* eQTLs from:
* [EMBL-EBI eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/)
