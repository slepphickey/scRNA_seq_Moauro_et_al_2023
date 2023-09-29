# Single-Cell RNA-seq analysis for Moauro et al 2023
This repository: ![](https://github.com/slepphickey/scRNA_seq_Moauro_et_al_2023/blob/main/zenodo_doi.svg)

These experiments are described in the paper:

[_OCT4 is expressed in extraembryonic endoderm stem (XEN) cell progenitors during somatic cell reprogramming_](link)  
Alexandra Moauro, Stephanie L. Hickey, Michael A. Halbisen, Anthony Parenti, and Amy Ralston*

*Correspondence: aralston@msu.edu

All files in the `data` folder, including the scRNA-seq raw count data, can be found on [Zenodo (doi.org/10.5281/zenodo.8388276)](https://doi.org/10.5281/zenodo.8388276)

## Processing Mouse Early Gastrulation scRNA-seq data from Mohammed et al, 2017

We compared data from [_Single-Cell Landscape of Transcriptional Heterogeneity and Cell Fate Decisions during Mouse Early Gastrulation, Mohammed et al. 2017_](http://www.sciencedirect.com/science/article/pii/S2211124717309610) with that from our 17 day OSKM reprogrammed MEFs. Preprocessing of the Mohammed et al data is described in the following documents:

**Rmarkdown**: `src/mohammed_preprocessing.Rmd`

**Main output**: `data/Mohammed_seurat.rdata`: Used in `Moauro_2023_Figure3.Rmd`

**Report**: `reports/mohammed_preprocessing.html`

## Quality control and filtering of 17 day OSKM reprogrammed MEF scRNA-seq data

Quality control and filtering of 17 day OSKM reprogrammed MEFs was performed as described in the following documents:

**Rmarkdown**: `src/scRNA_QC_filter.Rmd`

**Parameters**:

  * `working_dir`: Working directory  
  * `project_name`: Name of the project. A new directory with this name will be made in the working directory. The output files can be found here.  
  * `lab`: Name of the group that produced the scRNA-seq data, for the report header. 
  * `raw_data_dir`: Path to the directoy containing the cellranger output files:  
      + `matrix.mtx.gz` 
      + `barcodes.tsv.gz` 
      + `features.tsv.gz`  
  * `pMlo`: lower cutoff for percent reads from mitochondria genes  
  * `pMhi`: upper cutoff for percent reads from mitochondria genes    
  * `nFeatureHi`: upper cutoff for number of genes expressed per cell 
  * `nFeatureLo`: lower cutoff for number of genes expressed per cell 
  * `nCountLo`: lower cutoff for UMI counts per cell  
  * `nCountHi`: upper cutoff for UMI counts per cell

**Main Output**: `data/Day17_Reprogramming_scRNA_seq_QC_filtered_SeuratObject.Rdata`: a Seurat object containing the filtered UMI count data.  

**Report**: `reports/scRNA_QC_filter.html`  

## Normalization and clustering of 17 day OSKM reprogrammed MEF scRNA-seq data

17 day OSKM reprogrammed MEF scRNA-seq data were normalized using `scTransform` and clustered at multiple resolutions as described in the following documents. Either all cells were included in the analysis, or _Pou5f1_ or _Sox2_ positive cells were extracted before normalization and clustering.

**Rmarkdown**: `src/scRNA_scTransform_norm_and_cluster.Rmd`

**Parameters**:

  * `working_dir`: Working directory  
  * `project_name`: Name of the project. A new directory with this name will be made in the working directory if one does not already exist. The output files can be found here.  
  * `lab`: Name of the group that produced the scRNA-seq data, for the report header. 
  * `seuratOb`: QC filtered seurat object from scRNA_QC_filter.Rmd  
  * `goi`: Cells expressing this gene (UMI > 1) should be extracted. “none” if all cells should be included in the analysis.  

**Main output**:

  * `data/Day17_Reprogramming_scRNA_seq_QC_filtered_SeuratObject.Rdata`: a Seurat object containing the filtered UMI count data, normalized count data, and cluster assignments from all cells. Clusters at resolution 0.70 were renamed in this figure 3:
      + cluster 0 in this file = cluster 1 in Figure 3A-B. 
      + cluster 1 in this file = cluster 2 in Figure 3A-B. 
      + cluster 2 in this file = cluster 3 in Figure 3A-B. 
      + And so on 
  * `data/Day17_Reprogramming_scRNA_seq_Pou5f1_expressing_QC_filtered_SeuratObject.Rdata`: a Seurat object containing the filtered UMI count data, normalized count data, and cluster assignments from Pou5f1 expressing cells. Clusters at resolution 0.70 were renamed in this figure 3:  
      + cluster 0 in this file = cluster A in Figure 3E-F. 
      + cluster 1 in this file = cluster B in Figure 3E-F. 
      + cluster 2 in this file = cluster C in Figure 3E-F. 
      + And so on 
  * `data/Day17_Reprogramming_scRNA_seq_Sox2_expressing_QC_filtered_SeuratObject.Rdata`: a Seurat object containing the filtered UMI count data, normalized count data, and cluster assignments from Sox2 expressing cells. Clusters at resolution 0.70 were renamed in this figure 3:  
      + cluster 0 in this file = cluster R in Figure 3G-H. 
      + cluster 1 in this file = cluster S in Figure 3G-H. 
      + cluster 2 in this file = cluster T in Figure 3G-H. 
      + And so on   

**Reports**:

  * `reports/scRNA_scTransform_norm_and_cluster_all_cells_expressing_cells.html`
  * `reports/scRNA_scTransform_norm_and_cluster_Pou5f1_expressing_cells.html`
  * `reports/scRNA_scTransform_norm_and_cluster_Sox2_expressing_cells.html`

## Find cluster enirched genes at all cluster resolutions in 17 day OSKM reprogrammed MEF scRNA-seq data

**Rscript**: `scr/seuratFindAllMarkers_allRes.R`

**Arguments**:

  1. path to Seurat object 
  2. cluster assignment column prefix ie “SCT_snn_res” 
  3. adjusted pval filter  
  4. path to output directory  

**Main Output**:

  * `results/all_cells/SCT_snn_res.0.7_marker_genes.csv`: Using `data/Day17_Reprogramming_scRNA_seq_QC_filtered_SeuratObject.Rdata` as input. 
  * `results/Pou5f1_expressing_cells/SCT_snn_res.0.7_marker_genes.csv`: Using `data/Day17_Reprogramming_scRNA_seq_Pou5f1_expressing_QC_filtered_SeuratObject.Rdata` as input. 
  * `results/Sox2_expressing_cells/SCT_snn_res.0.7_marker_genes.csv`: Using `data/Day17_Reprogramming_scRNA_seq_Sox2_expressing_QC_filtered_SeuratObject.Rdata` as input. 

## Wamaitha et. al. mESC vs eXEN DE genes

Obtaining genes differentially expressed between mESCs and eXEN cell beadarray expression data generated as described in:

Wamaitha *et. al.*, 2015, *Gata6 potently initiates reprograming of pluripotent and differentiated cells to extraembryonic endoderm stem cells*, Genes & Dev. 2015. 29: 1239-1255, [doi:10.1101/gad.257071.114](https://genesdev.cshlp.org/content/29/12/1239.long#sec-9)

**Rmarkdown**: `src/Wamaitha_mESC_vs_eXEN.Rmd`

**Output**: `results/Wamaitha_mESC_vs_eXEN_DE.csv`

**Report**: `reports/Wamaitha_mESC_vs_eXEN.html`

## Figure 3 plots

Code for generating the plots in Figure 3 using the output files described above

**Rmarkdown**: `src/Moauro_2023_Figure3.Rmd`  

**Main Output**

  * `results/figures/fig3ab.png`: A) UMAP visualization of clusters with all cells. B) Over-representation analysis of cluster enriched genes with Mohammed et al cluster enriched genes  
  * `results/figures/fig3d.png`: D) Over-representation analysis of Mohammed et al cluster enriched genes with genes up-regulated in ES vs XEN cells and XEN vs ES cells from Wamaitha et al.  
  * `results/figures/fig3ef.png`: E) UMAP visualization of clusters with _Pou5f1_ expressing cells. F) Over-representation analysis of cluster enriched genes with Mohammed et al cluster enriched genes  
  * `results/figures/fig3ef.png`: G) UMAP visualization of clusters with _Sox2_ expressing cells. H) Over-representation analysis of cluster enriched genes with Mohammed et al cluster enriched genes  
  * `results/figures/all_cells_mohammed_shared_genes.csv`: Mohammed et al cluster enriched genes shared with 17 day OSKM reprogrammed MEF cluster enriched genes from all cells.  
  * `results/figures/Oct4_mohammed_shared_genes.csv`: Mohammed et al cluster enriched genes shared with 17 day OSKM reprogrammed MEF cluster enriched genes from _Pou5f1_ expressing cells.  
  * `results/figures/Sox2_mohammed_shared_genes.csv`: Mohammed et al cluster enriched genes shared with 17 day OSKM reprogrammed MEF cluster enriched genes from _Sox2_ expressing cells.  
  * `results/figures/SCT_snn_res.0.7_marker_genes_all_cells_figure3_clusters.csv`: same as `results/all_cells/SCT_snn_res.0.7_marker_genes.csv` with cluster names that match figure 3. 
  * `results/figures/SCT_snn_res.0.7_marker_genes_Pou5f1_expressing_cells_figure3_clusters.csv`: same as `results/Pou5f1_expressing_cells/SCT_snn_res.0.7_marker_genes.csv` with cluster names that match figure 3. 
  * `results/figures/SCT_snn_res.0.7_marker_genes_Sox2_expressing_cells_figure3_clusters.csv`: same as `results/Sox2_expressing_cells/SCT_snn_res.0.7_marker_genes.csv` with cluster names that match figure 3. 
