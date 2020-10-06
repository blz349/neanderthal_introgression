# Import modules
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
import numpy as np

## Get Chen Neanderthal SNPs in Fairfax recomputed eQTLs
chen = pd.read_excel('chen/Additional File 1.xlsx', 'Sheet1')
chen = chen.loc[(chen.Chen == 'Yes') & (chen['Fairfax recomputed eQTL'] == 'Yes')].copy()  # 9887
chen.drop(columns=['Chen', 'Fairfax eQTL', 'Fairfax recomputed eQTL', 'GWAS Catalog'], inplace=True)

# Keep most significant eQTL per gene
refairfax_nean = pd.read_csv('GWAS/overlap_filtered_fairfax.csv')
refairfax_nean.sort_values('pvalue', inplace=True)
refairfax_nean.drop_duplicates(subset='gene_id', keep='first', inplace=True)  # 418 genes, 390 rsids

# Get list of genes
ensembl_info = pd.read_csv('fairfax/hg38.ensGene.txt', usecols=['gene_id', 'symbol'])
symbol = refairfax_nean.merge(ensembl_info, how='left', on=['gene_id'])
symbol[['symbol', 'pvalue']].to_csv('fairfax/pathway/all_genes.txt', index=False, header=False, sep='\t')

ifn = symbol.loc[symbol.Condition == 'IFN'][['symbol', 'pvalue']].copy()
ifn[['symbol', 'pvalue']].to_csv('fairfax/pathway/ifn_genes.txt', index=False, header=False, sep='\t')

lps2 = symbol.loc[symbol.Condition == 'LPS2'][['symbol', 'pvalue']].copy()
lps2[['symbol', 'pvalue']].to_csv('fairfax/pathway/lps2_genes.txt', index=False, header=False, sep='\t')

lps24 = symbol.loc[symbol.Condition == 'LPS24'][['symbol', 'pvalue']].copy()
lps2[['symbol', 'pvalue']].to_csv('fairfax/pathway/lps24_genes.txt', index=False, header=False, sep='\t')

naive = symbol.loc[symbol.Condition == 'Naive'][['symbol', 'pvalue']].copy()
naive[['symbol', 'pvalue']].to_csv('fairfax/pathway/naive_genes.txt', index=False, header=False, sep='\t')

# Plot single pathways
enrich = pd.read_csv('fairfax/pathway/all_naive_enrichment.txt', sep='\t')
enrich['-log10(FDR)'] = -np.log10(enrich['adjp'])
enrich.sort_values("-log10(FDR)", ascending=False)
to_plot = enrich.iloc[0:15]
# sns.set_style("white")
# sns.set_style("ticks")
plt.figure(figsize=(2, 5))
ax = sns.barplot(x="-log10(FDR)", y="name", data=to_plot, dodge=False)  # palette='husl'
ax.set(xlabel='-log10(FDR)', ylabel='')
ax.set_yticklabels('\n'.join(textwrap.wrap(y.get_text(), 45)) for y in ax.get_yticklabels())
#ax.set_yticklabels((textwrap.fill(y.get_text(), 35) for y in ax.get_yticklabels()))
ax.set_title('Naive')
plt.savefig('fairfax/pathway/all_naive_pathways.pdf', bbox_inches='tight', pad_inches=0)
plt.show()
plt.clf()

# Plot all pathways
colors = ['tomato', 'mediumseagreen', 'cornflowerblue', 'goldenrod']
context = ['Naive', 'IFN', 'LPS2', 'LPS24']
i = 0
fig, axes = plt.subplots(nrows=2, ncols=2, constrained_layout=True)
axes = axes.flatten()
for ax in axes:
    grp = context[i]
    enrich = pd.read_csv('fairfax/pathway/all_' + str(grp).lower() + '_enrichment.txt', sep='\t')
    enrich['-log10(FDR)'] = -np.log10(enrich['adjp'])
    enrich.sort_values("-log10(FDR)", ascending=False)
    to_plot = enrich.iloc[0:15]
    sns.barplot(x="-log10(FDR)", y="name", data=to_plot, ax=ax, dodge=False, color=colors[i])
    ax.set(ylabel='')
    ax.set_title(grp)
    i += 1
plt.show()
plt.savefig('fairfax/pathway/all_pathways.pdf', bbox_inches='tight', pad_inches=0)

# Plot neanderthal pathways
colors = ['tomato', 'mediumseagreen']
context = ['Naive', 'IFN']
i = 0
fig, axes = plt.subplots(nrows=2, constrained_layout=True, gridspec_kw={'height_ratios': [2, 1]})
for ax in axes:
    grp = context[i]
    enrich = pd.read_csv('fairfax/pathway/nean_' + str(grp).lower() + '_enrichment.txt', sep='\t')
    enrich['-log10(FDR)'] = -np.log10(enrich['adjp'])
    enrich.sort_values("-log10(FDR)", ascending=False)
    to_plot = enrich.iloc[0:15]
    sns.barplot(x="-log10(FDR)", y="name", data=to_plot, ax=ax, dodge=False, color=colors[i])
    ax.set(ylabel='')
    ax.set_title(grp)
    i += 1
plt.show()
plt.savefig('fairfax/pathway/nean_pathways.pdf', bbox_inches='tight', pad_inches=0)

