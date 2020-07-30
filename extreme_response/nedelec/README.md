# Nedelec code

### Requirements:
Python 3 and modules:
* pandas 1.0.3
* seaborn 0.10.0
* Matplotlib 3.1.1
* SciPy 1.4.1
* plinkio 0.9.7
* networkx 2.4

R and packages:
* XGR
* RCircos

---

### Description of code

    neanderthal_nedelec_overlap.ipynb
This code looks at the overlap in SNPs between Neanderthal-introgressed SNPs (from Dannemann *et al.* and Simonti *et al.*) and cis eQTLs, reQTLs and asQTLs identified by Nedelec *et al* using macrophages that were *Listeria*-infected, *Salmonella*-infected or non-infected.  
&nbsp;

    plot_nedelec_eQTLs.ipynb
This code plots the expression data from African and European individuals for Nedelec condition-specific eQTLs for which the same variant has also been identified as a Neanderthal-introgressed SNP. Gene expression levels were compared between African and European individuals via a two-sided T-test. Plots were divided by the 3 experimental conditions: *Listeria*-infected, *Salmonella*-infected and Non-infected.   
&nbsp;

    neanderthal_recomputed_nedelec_overlap.md
This document contains code for determining the overlap between Neanderthal-introgressed SNPs (from Dannemann *et al.* and Simonti *et al.*) and recomputed Nedelec *et al.* eQTLs from the EMBL-EBI eQTL Catalogue. Overlapping eQTLs were filtered for p < 10<sup>-8</sup> and saved as CSV format (coordinates in GRCh38) under: `/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/overlap_all_filtered.csv/`. Genes associated with the eQTLs were extracted and saved as text files under: `/well/jknight/shiyao/data/nedelec_expression/EMBL_recomputed/genes_*.txt`.  
&nbsp;

    plot_recomputed_nedelec_networks.ipynb
This code plots the network of genes associated with eQTLs, which have been identified as Neanderthal-intorgressed, using subnet results obtained from the R package XGR, for each of the 3 conditions: *Listeria*, *Salmonella* and Naive.  
&nbsp;

---

### Datasets
Neanderthal-introgressed SNPs from:
1. Dannemann M, Prufer K & Kelso J. Functional implications of Neandertal introgression in modern humans. *Genome Biol* 2017 **18**:61.
2. Simonti CN *et al.* The phenotypic legacy of admixture between modern humans and Neandertals. *Science* 2016 **351**:737-41.  
&nbsp;

MAFs from:
* [1000 genomes](https://www.internationalgenome.org/data/)  
&nbsp;

Macrophage gene-expression data from:
* Nedelec Y *et al.* Genetic Ancestry and Natural Selection Drive Population Differences in Immune Responses to Pathogens. *Cell* 2016 **167**:657-669e21.  
&nbsp;

Recomputed Nedelec *et al.* eQTLs from:
* [EMBL-EBI eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/)  
