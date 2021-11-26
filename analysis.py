#!/usr/bin/env python
# coding: utf-8
import numpy as np
import pandas as pd
import scanpy as sc
import sys

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

results_file = 'write/covtest.h5ad'  # the file that will store the analysis results
input_path = sys.argv[1]
adata = sc.read_h5ad(input_path)                              
# write a cache file for faster subsequent reading
# # Basic filtering:
adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=1)

# Let's assemble some information about mitochondrial genes, which are important for quality control.
# 
# Citing from "Simple Single Cell" workflows [(Lun, McCarthy & Marioni, 2017)](https://master.bioconductor.org/packages/release/workflows/html/simpleSingleCell.html#examining-gene-level-metrics):
# 
# > High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells. The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.
vir_index = [n.find("ENSG")!=0 for n in adata.var_names]
virus_genemanes = adata.var_names[vir_index]

# With `pp.calculate_qc_metrics`, we can compute many metrics very efficiently.

# A violin plot of some of the computed quality measures:
# 
# * the number of genes expressed in the count matrix
# * the total counts per cell
# * the percentage of counts in mitochondrial genes
adata.var['virus'] = adata.var_names.str.startswith('gene-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['virus'], percent_top=None, log1p=False, inplace=True)
adata.obs["viral_positive"] = adata.obs["total_counts_virus"] >0
adata.obs = adata.obs.astype({"viral_positive":int})

# Total-count normalize (library-size correct) the data matrix $\mathbf{X}$ to 10,000 reads per cell, so that counts become comparable among cells.
sc.pp.normalize_total(adata, target_sum=1e4)

# Logarithmize the data:
sc.pp.log1p(adata)

# Identify highly-variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0)
adata.raw = adata

# Actually do the filtering
adata = adata[:, adata.var.highly_variable]

# Scale each gene to unit variance. Clip values exceeding standard deviation 10. 
sc.pp.scale(adata, max_value=10)


# ## Principal component analysis

# Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
sc.tl.pca(adata, svd_solver='arpack')
adata.write(results_file)


# ## Computing the neighborhood graph

# Let us compute the neighborhood graph of cells using the PCA representation of the data matrix. You might simply use default values here. For the sake of reproducing Seurat's results, let's take the following values.
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# ## Embedding the neighborhood graph

# We suggest embedding the graph in two dimensions using UMAP ([McInnes et al., 2018](https://arxiv.org/abs/1802.03426)), see below. It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preserves trajectories. In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:
# 
# ```
# tl.paga(adata)
# pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
# tl.umap(adata, init_pos='paga')
# ```
sc.tl.umap(adata)
pd_meta = pd.read_csv("meta.csv",index_col=0)
adata.obs = adata.obs.merge(pd_meta,left_index=True, right_index=True,how="left")
adata.obs.loc[adata.obs["viral_positive"]==1].to_csv("viral_cells.csv")

# As we set the `.raw` attribute of `adata`, the previous plots showed the "raw" (normalized, logarithmized, but uncorrected) gene expression. You can also plot the scaled and corrected gene expression by explicitly stating that you don't want to use `.raw`.

# ## Clustering the neighborhood graph

# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by [Traag *et al.* (2018)](https://scanpy.readthedocs.io/en/latest/references.html#traag18). Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.
sc.tl.leiden(adata,resolution=0.08)


# Plot the clusters, which agree quite well with the result of Seurat.
sc.pl.umap(adata, color=['celltype','leiden','viral_positive'],save="umapplot.pdf")


# ## Finding marker genes

# The result of a [Wilcoxon rank-sum (Mann-Whitney-U)](https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test) test is very similar. We recommend using the latter in publications, see e.g., [Sonison & Robinson (2018)](https://doi.org/10.1038/nmeth.4612). You might also consider much more powerful differential testing packages like MAST, limma, DESeq2 and, for python, the recent diffxpy.
try:
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,save="deg.pdf")
finally:
    print("Cannot find degs")
adata.write(results_file)

