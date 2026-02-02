# Welcome to the SPACE!
# Spatial Profiling Across Cancer Ecosystems (SPACE)
This repository contains codes for the project: "Single-cell spatial profiling of the tumor core microenvironment reveals spatially-resolved cellular states, molecular programs and therapeutic targets".

For more details, see our publication: ""

Codes provided here are to perform:   
- Step1-Preprocessing_QC_and_cell_annotation.R: Data pre-processing, quality control, and cell annotation of Xenium 5K in situ spatial transcriptomics data.  
- Step2-Spatial_niche_identification.R: Identification of 5 distinct spatial niches.  
- Step3-Neighboring_cell_identification.R: Identify neighboring cell and calculate neighboring cell composition in each cancer type.  
- Step4-Cell_Colony_Analysis.R: Perform cell colony size analysis across cancer types using R package dbscan.  
- Step5-ICB_ORR_association.R: Perform correlation analysis between cell proportion and objective response rate in cancer immune checkpoint therapy.  
- Step6-Spatial_DEGs_identification.R: Identification of spatially enriched genes (SEGs) across cancer types, cell types, and spatial niches.  
- Step7-TF_analysis.R: Identification of spatially enriched TFs (SETFs) and differentially activated TFs (DATFs) in STAR-T cells compared to other-T cells.  
- Step8-Metabolism.R: Perform AUCell to evaluate metabolic activity based on KEGG dataset and Reactome dataset.  
- Step9-Cell_cell_communication.R: Identification of spatially enriched Ligand-Receptor pairs (SELRPs) across cell types and spatial niches, and differentially enriched Ligand-Receptor pairs in STAR-T cells compared to other-T cells.  
- Step10-NMF_analysis.ipynb: Identification of spatial clustering using neighbor vector-based clustering method.  
- Step11-Spatial_gradient_gene_expression_analysis.R: Identification of transcriptomic hallmarks of tumor cells linking to T-cell proximity.

# Downloading the data
- The raw data of Xenium 5K in situ spatial transcriptomics has been deposited in OMIX (https://ngdc.cncb.ac.cn/omix/releaseList) under .  
- The processed data (RDS) of Xenium 5K in situ spatial transcriptomics has been deposited in OMIX (https://ngdc.cncb.ac.cn/omix/releaseList) under OMIX014642.  
- The raw data of bulk RNA-seq profiling of INSR-overexpressing and non-targeting control HUVECs has been deposited in GSA-Human (https://ngdc.cncb.ac.cn/gsa-human/) under .
- The processed gene counts matrix of bulk RNA-seq profiling of INSR-overexpressing and non-targeting control HUVECs has been deposited in OMIX (https://ngdc.cncb.ac.cn/omix/releaseList) under OMIX014682.

- # Operation systems
- Linux version 5.4.0-42-generic x86_64
- Windows 10 64-bit
- # Hardware requirements
Due to the large number of cells in the dataset, the pipeline requires a Linux platform with sufficient RAM to support the computational workload. 
For optimal performance, we used a Linux server with the following specifications: an Intel® Xeon® Gold 5318Y CPU @ 2.10 GHz and 1024 GB of RAM.
# Installation guide
R packages required for the pipeline can be installed from CRAN (https://cran.r-project.org/) using the install.packages() function, or from Bioconductor (https://bioconductor.org/) using the BiocManager::install() function.
We used conda to manage python packages, python packages can be installed with conda install --- and pip install ---.
pySCENIC used for transcription factor analysis can be installed from https://github.com/aertslab/pySCENIC.
# Packages
# Pyhton and python packages:
- python (version 3.11, version 2.8.20 for pySCENIC)
- pySCENIC (version 0.12.1)
- matplotlib (version 3.8.0)
- numpy (version 1.24.4)
- pandas (version 1.5.3)
- nimfa (version 1.4.0)
- scikit-learn (version 1.3.1)
- scanpy (version 1.9.8)
# R base and packages:
- R-base (version 4.4.0)
- Seurat (version 5.3.0)
- ggsci (version 3.2.0)
- RColorBrewer (version 1.1.3)
- ggplot2 (version 3.5.2)
- dplyr (version 1.1.4)
- viridis (version 0.6.5)
- BiocParallel (version 1.40.2)
- patchwork (version 1.3.0)
- ComplexHeatmap (version 2.22.0)
- clustertree (version 0.5.1)
- tidyverse (version 2.0.0)
- paletteer (version 1.6.0)
- future (version 1.49.0)
- presto (version 1.0.0)
- ggpubr (version 0.6.0)
- SeuratDisk (version 0.0.0.9021)
- clusterProfiler (version 4.14.4)
- ReactomePA (version 1.50.0)
- ggrepel (version 0.9.6)
- AUCell (version 1.28.0)
- GSEABase (version 1.68.0)
- GSVA (version 1.52.3)
- KEGGREST (version 1.46.0)
- ggSCvis (version 0.0.3)
- CellChat (version 2.1.2)
- edgeR (version 4.4.2)
# Other softwares
- Xenium Explorer (version 4.4.0)
