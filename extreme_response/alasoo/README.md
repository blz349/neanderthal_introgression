# Alasoo code

### Requirements:
Python 3 and modules:
* pandas 1.0.3
* Matplotlib 3.1.1
* networkx 2.4

R and packages:
* XGR
* RCircos

---

### Description of code

    neanderthal_recomputed_alasoo_overlap.md
This document contains code for determining the overlap between Neanderthal-introgressed SNPs (from Dannemann *et al.* and Simonti *et al.*) and recomputed Alasoo *et al.* eQTLs from the EMBL-EBI eQTL Catalogue. Overlapping eQTLs were filtered for p < 10<sup>-8</sup> and saved as CSV format (coordinates in GRCh38) under: `/well/jknight/shiyao/data/alasoo/overlap_all_filtered.csv/`. Genes associated with the eQTLs were extracted and saved as text files under: `/well/jknight/shiyao/data/alasoo/genes_*.txt`.  
&nbsp;

    plot_recomputed_alasoo_networks.ipynb
This code plots the network of genes associated with eQTLs, which have been identified as Neanderthal-intorgressed, using subnet results obtained from the R package XGR, for each of the 4 conditions: IFN, *Salmonella*, IFN+*Salmonella* and Naive.  
&nbsp;

### Datasets
Neanderthal-introgressed SNPs from:
1. Dannemann M, Prufer K & Kelso J. Functional implications of Neandertal introgression in modern humans. *Genome Biol* 2017 **18**:61.
2. Simonti CN *et al.* The phenotypic legacy of admixture between modern humans and Neandertals. *Science* 2016 **351**:737-41.  
&nbsp;

MAFs from:
* [1000 genomes](https://www.internationalgenome.org/data/)  
&nbsp;

Recomputed Alasoo *et al.* eQTLs from:
* [EMBL-EBI eQTL Catalogue](https://www.ebi.ac.uk/eqtl/Data_access/)  
