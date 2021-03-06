# Hotspot Analysis Code

This is the analysis repository which accompanies the manuscript for [Hotspot](https://github.com/yoseflab/hotspot)

## Directory Contents

- `/data`: Downloading and formatting of public data
- `/pipelineScripts`: Shared scripts for Snakemake pipelines
- `/SlideSeq`: Analysis code for the Slide-Seq example (mouse cerebellum)
- `/Lineage`: Analysis code for the Lineage example (mouse embryogenesis)
- `/Transcriptomics`: Analysis code for the CD4 T cell and Simulation examples
- `/Notebooks`: Example usage notebooks

## Notes

Pipelines here make heavy use of [Snakemake](https://snakemake.readthedocs.io/en/stable/).

The `Snakefile` in directories describes the scripts and steps needed to run the pipeline.

Outputs can be recreated in each directory by running `snakemake all`

## Figures

Code to re-create figures can be found in various `Figures` directories and depends on prior execution of Snakemake pipelines.

- Figure 1 - Algorithm Diagram

- Figure 2
    - Panel A: `/Transcriptomics/Figures/EvaluateFeatureSelection/plotRelevance.py`
    - Panel B: `/Transcriptomics/Figures/EvaluateFeatureSelection/compareLatentSpaces.py`
    - Panel C: `/Transcriptomics/Figures/EvaluateFeatureSelection/compare_local_expression.py`
    - Panel D: `/Transcriptomics/Figures/CD4_Correlation/plot.py`
    - Panel E: `/Transcriptomics/Figures/CD4_Correlation/ModuleConsistency_CD4.py`

- Figure 3
    - Panel A: `/SlideSeq/Figures/MainFigure/moduleHeatmap.py`
    - Panel B: `/SlideSeq/Figures/MainFigure/moduleCellTypes.py`
    - Panel C: `/SlideSeq/Figures/MainFigure/moduleScores.py`

- Figure 4
    - Panel A: `/Lineage/Figure/plotCorrelations.py`
    - Panel B: `/Lineage/Figure/plotUMAPs.py`
    - Panel C: `/Lineage/Figure/plotKernels.py`
    - Panel D: `/Lineage/Figure/plotKernelsTx.py`
    - Panel E: `/Lineage/Figure/plotAngioblasts.py`

- Figure S1
    - Panel A: `/Transcriptomics/Figures/Simulation/plotTSNEs.py`
    - Panel B: `/Transcriptomics/Figures/Simulation/plotAUC_PR.py`
    - Panel C: `/Transcriptomics/Figures/Simulation/plotModuleAssignment.py`

- Figure S2
    - Panel A: `/SlideSeq/Figures/Supp1/plotMeanVar.py`, `/SlideSeq/Figures/Supp_Autocorr/compare_local_expression.py`
    - Panel B: `/SlideSeq/Figures/Supp1/plotPR.py`
    - Panel C: `/SlideSeq/Figures/Supp1/plotTimings.py`
    - Panel D: `/SlideSeq/Figures/Supp1/plotIDR.py`

- Figure S3
    - Panel A: `/SlideSeq/Figures/Supp2/comparePairwiseZScores.py`
    - Panel B: `/SlideSeq/Figures/Supp2/compareModules.py`
    - Panel C: `/SlideSeq/Figures/Supp2/compareModuleAssignments.py`

- Figure S4
    - Panel A: `/SlideSeq/Figures/Supp3/compareModulesSpatialDE.py`
    - Panel B: `/SlideSeq/Figures/Supp3/compareTiming.py`
    - Panel C: `/SlideSeq/Figures/Supp3/compareModulesSpatialDE.py`

- Figure S5
    - Panel A: `/SlideSeq/Figures/Supp4/plotPValues.py`
    - Panel B: `/Transcriptomics/Figures/Supp_Pvals/plotPValues.py`
    - Panel C: `/Lineage/Figure/plotPValues.py`

- Figure S6
    - All Panels: `/Transcriptomics/Figures/EvaluateFeatureSelection/hvg_vs_hs.py`

- Figure S7
    - Panel A: `/Transcriptomics/Figures/Simulation/downsampling_correlation.py`
    - Panel B: `/Transcriptomics/Figures/CD4_Correlation/plot_downsampled.py`
    - Panel C: `/Transcriptomics/Figures/CD4_Correlation/plot_downsampled.py`
    - Panel D: `/Transcriptomics/Figures/Simulation/downsampling_correlation.py`

- Figure S8
    - All Panels: `/SlideSeq/Figures/Supp5_HVG_vs_HS/hvg_vs_hs.py`

- Figure S9
    - Panel A: `/SlideSeq/Figures/Supp6_NegBinom_vs_Bernoulli/plotPR.py`
    - Panel B: `/SlideSeq/Figures/Supp6_NegBinom_vs_Bernoulli/compareModuleAssignments.py`

- Figure S10
    - Column 1: `/Transcriptomics/Figures/Simulation/plotAUC_PR_k_sensitivity.py`
    - Column 2: `/Transcriptomics/Figures/EvaluateFeatureSelection/plotRelevance_k_sensitivity.py`
    - Column 3: `/SlideSeq/Figures/Supp1/plotPR_k_sensitivity.py`

- Figure S11
    - Panel A: `/SlideSeq/Figures/Supp7_Bernoulli/plot.py`
    - Panel B: `/Transcriptomics/Figures/CD4_Correlation/ModuleConsistency_Monocytes.py`

- Figure S12
    - All Panels: `/Transcriptomics/Figures/EvaluateFeatureSelection/compare_local_expression.py`

- Figure S13
    - All Panels: `/SlideSeq/Figures/Supp_Autocorr/compare_local_expression.py`


## Software Versions

**Python 3.6.8**

`biopython==1.72`  
`ete3==3.1.1`  
`feather-format==0.4.0`  
`h5py==2.9.0`  
`loompy==2.0.17`  
`matplotlib==3.1.0`  
`naivede==1.2.0`  
`numba==0.45.0`  
`numpy==1.16.5`  
`pandas==0.25.1`  
`scanpy==1.4`  
`scipy==1.2.1`  
`scvi==0.2.4`  
`seaborn==0.9.0`  
`scikit-learn==0.21.2`  
`spatialde==1.1.1`  
`snakemake==5.5.3`  
`statsmodels==0.10.0`  
`tqdm==4.32.2`  
`umap-learn==0.3.6`  
`hotspot==0.9.0` ([github.com/yoseflab/hotspot](https://github.com/yoseflab/hotspot))  
`bio_utils` (included in /packages directory)  
`gene_enrich` (included in /packages directory)  

**R 3.5.1**

`edgeR==3.22.3`  
`feather==0.3.5`  
`ggplot2==3.1.0`  
`idr==1.2`  
`jsonlite==1.6`  
`loomR==0.2.0`  
`M3Drop==3.10.4`  
`matrixStats==0.54.0`  
`NMF==0.21.0`  
`openxlsx==4.1.0`  
`pbmcapply==1.3.1`  
`Rtsne==0.15`  
`Seurat==2.3.4`  
`SymSim==0.0.0.9000` ([github.com/yoseflab/symsim](https://github.com/yoseflab/symsim))  
`VISION==2.0.0` ([github.com/yoseflab/VISION](https://github.com/yoseflab/VISION))  
`DropSeq.util==2.0` (Available via [DropViz](dropviz.org))
