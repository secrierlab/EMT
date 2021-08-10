# Epithelial-to-mesenchymal trajectory quantification in cancer
### Guidantonio Malagoli Tagliazucchi, Maria Secrier 
UCL Genetics Institute, Department of Genetics,  Evolution and Environment,  University College London,  UK 

<p align="center">
  <img width="700" height="400" src="https://github.com/secrierlab/EMT/blob/main/figrepo.png">
</p>

# Table of contents

## System Requirements
Operating system(s): Unix (linux, mac)
Programming Language: R
All the analysis have been run on a server with 64 CPU (Intel(R) Xeon(R) Gold 5218 CPU @ 2.30GHz) and 376G of RAM.

## Installation
- Install R >= 3.6.1, available on https://cran.r-project.org/
- Install RStudio (Optional), available on https://www.rstudio.com/products/rstudio/download/
- An install file with the R-libraries to install and dependencies is provided (INSTALL.R). The user can copy and paste the following line on a R shell.
- See section "Data" to get more information about the sources of data and how retrieve them

Here the procedure:

```r
# Move in the directory in which there is a script
setwd("/Users/username/Desktop/EMT-main")

# run the script
source("install.R")
```
On a normal laptop/desktop computer (2.4 GHz Quad-Core Intel Core i5, 8 GB 2133 MHz LPDDR3), the installation could take 20-30 minutes.

## Data
The majority of the data used in this work are in the folde /data. In addition, we provide a script to download the data from several repositories.
The references of the datasets used in this work are reported in the section "References" of this page, and in the main manuscript. 

Steps:

1) run the bash scripts inside the /data folder as follow:
```console
foo@bar:~$ sh decompress_data.sh
foo@bar:~$ sh download_data.sh
```
## How to run a script

1) Check that all the depencies are installed (INSTALL.R script)
2) Open RStudio or shell
3) Move in the path with the script of interest as follow. 

For example if the user wants to run the segmentation (run_hmm.R) on a toy dataset of gene-expression (TCGA) these are the steps to use:

```r
# Move in the directory in which there is a script
setwd("/Users/username/Desktop/EMT-main/02_HMM_macrostates_EMT")

# run the script
source("run_hmm.R")
```
The example reported above employs less than 5 min to run on a normal laptop/desktop computer (2.4 GHz Quad-Core Intel Core i5, 8 GB 2133 MHz LPDDR3).
The expected output file for this toy example are several charts (with the results of the segmentation) and tables. 

**Note:**, this is toy example, the biological assignments (e.g. pEMT) of the HMM states could be different respect the pan-cancer analysis. See material and methods for more details.
**Note:** change the string "EMT-main" to "EMT-EMTquant.v*.*", if you download the data from the "Releases" section (menu on the left of this page)


**Notes:**: 
- The scripts are able to upload the data (in the folder /data) that they must to use automatically.
- All the output of the scripts are saved in the folder (output_dir) 
 
## ReconstructionEMTbulk
This folder contains the code to reconstruct the EMT trajectory of the bulk RNA-seq data from TCGA
- **projection_EMT_trajectory_MCF10A_to_TCGA.R:** main script to reconstruct the EMT trajectory of TCGA samples using as reference single-cell MCF10 data and their defined EMT pseudospace (trajectory). This script performs also several exploratory analysis. 
- **run_analysis_pseudospace_TCGA_MET500.R:** procedure to quantify the EMT trajectory using TCGA and MET500 data
- **projection_EMT_trajectory_MCF10A_to_CCLE.R:** main script to reconstruct the EMT trajectory of CCLE samples using as reference single-cell MCF10 data and their defined EMT pseudospace (trajectory). This script performs also several exploratory analysis. This script combine also EMT states with the information of the metastatic potential.
- **knnForCCLE_MCF.R:** functions useful to run projection_EMT_trajectory_MCF10A_to_CCLE.R
- **petalChartGMT.R:** functions useful to run projection_EMT_trajectory_MCF10A_to_CCLE.R
 

## HMM_macrostates_EMT
This folder contains the code to perform the segmentation of the EMT trajectory and to identify the macrostates of EMT
- **run_hmm.R:** This script contain the code to run the segmentation and perform the identification of the EMT states
- **HMM_noise.R:** This script implements a procedure to add noise (through the jitter function) in the original datasets and quantify the robustness of the HMM states 
  
## extrinsic_hallmarks
This folder contains the code to quantify the extrinsic hallmarks of EMT
- **run_fibroblast_deconvolution.R:** This is the script used to quantify the levels of fibroblast infiltration, contains also the code to run the multinomial logistic regression with the fibroblasts signals and the HMM states.
- **run_tme_deconvolution.R:** This is the script used to quantify the cell components of Tumour Micro Environment (TME), contains also the code to run the multinomial logistic regression with the TME signals and the HMM states.
- **run_stemness_broad_analysis.R:** script to compute the stemness scores using different gene-sets.
- **characterize_extrinsic_hallmarks.R:** In this script there is the procedure to quantify the hypoxia scores (according to Buffa et al.). This script create also a data.frame that combine the quantification of several others hallmarks (e.g. centromeric amplification (CA20), aneuploidy, etc)

## intrinsic_hallmarks
This folder contains the code to quantify the intrinsic hallmarks of EMT
- **mutational_signatures.R:** Script to perform the identification of the mutational signatures in each EMT states with mixed-effect models (Buffa et al.)

## biomarkers_EMT
This folder contains the pipeline used to identify the biomarkers of EMT
- **run_pipeline_ml.sh:** This script contains the pipeline to perform the identification of the biomarkers. Please open the .sh and infofile.txt to get more details about the procedure.
- **multiOmicsCharts.sh:** This script plots the results of the biomarkers discovered

## tissue_specific_trajectories
This folder contains the pipeline used to identify tissue specifics trajectories in LUAD and BRCA, and the code to define the markers associated with EMT states
- **ts_analysis_find_driver_events_CNV.R:** find drivers copy number alteration in LUAD and BRCA
- **ts_analysis_find_driver_events_MUT.R:** find drivers mutations in LUAD and BRCA
- **plot_cook_states_dndscv_Selected.RL:** plot only the relevant mutational events in each EMT state

## clinical_EMT
This folder contains the code to perform the characterization of the clinical features related with EMT. The folder contain also the code to analyze pharmacogenomic datasets (GDSC, depmap, POG570)

- **endpoints_EMT.R:** Code to run the overall survival analysis and cox proportional hazard models.
- **run_demographic_table.R:** Create a demographic table and perform other statitics.
- **integrate_EMT_with_drugs.R:** Integrate the HMM states, EMT scores with data from GDSC
- **aov_depmap_ml_features.R:** Integrate the HMM states, EMT scores with data from DepMap
- **endpoints_response_to_theraphyApril.R:** Script to create the Kaplan Meier of PFI for significant therapies
- **CoxPh_model_mutations.R:** Script to find the significant mutations impacting the overall survival of patients
- **knnForPOG_MCF.R.R:** kNN procedure with POG dataset
- **knnForGDSC_MCF.R:** kNN procedure with GDSC dataset
- **compare_before_after_treatment_TCGA.R:** compare the levels of EMT of samples naive treated (TCGA) and treated (POG)

## Validation
This folder contains the code to perform the validation using external datasets.
- **validate_biomarkers_with_CCLE_and_METMAP_MUT_v2.R:** Compare the metastatic potential of samples with and without mutations in the biomarkers
- **validate_biomarkers_with_CCLE_and_METMAP_CNV_v2.R:** Compare the metastatic potential of samples with and without copy number alterations in the biomarkers
- **CRISPR_validation_MUT.R:** Get the Ceres Scores of the biomarkers associated with mutations and plot heatmaps
- **CRISPR_validation_CNV.R:** Get the Ceres Scores of the biomarkers associated with copy number alterations and plot heatmaps
- **validation_Schaller_MUT.R:** This script check if putative biomarkers of EMT with mutations are targets of TFs important for EMT
- **validation_Schaller_CNV.R:** This script check if putative biomarkers of EMT associated with CNV are targets of TFs important for EMT

# How to cite
At present, a version of this manuscript is available on bioarXiv: https://www.biorxiv.org/content/10.1101/2021.07.23.453584v1.

# Reference datasets

The results published here are in part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga. The authors would like to acknowledge the American Association for Cancer Research and its financial and material support in the development of the AACR Project GENIE registry, as well as members of the consortium for their commitment to data sharing. Interpretations are the responsibility of study authors.

-	AACR Project GENIE Consortium. AACR Project GENIE: Powering Precision Medicine through an International Consortium. Cancer Discov 7, 818–831 (2017).
- Bhandari, V. et al. Molecular landmarks of tumor hypoxia across cancer types. Nature Genetics 51, 308–318 (2019).
- Bhandari, V., Li, C. H., Bristow, R. G., Boutros, P. C., & PCAWG Consortium. Divergent mutational processes distinguish hypoxic and normoxic tumours. Nat Commun 11, 737 (2020).
- Iorio, F. et al. A Landscape of Pharmacogenomic Interactions in Cancer. Cell 166, 740–754 (2016).
- Jin, X. et al. A metastasis map of human cancer cell lines. Nature 588, 331–336 (2020).
- Miranda, A. et al. Cancer stemness, intratumoral heterogeneity, and immune response across cancers. Proc Natl Acad Sci U S A 116, 9020–9029 (2019).
- Robinson, D. R. et al. Integrative clinical genomics of metastatic cancer. Nature 548, 297–303 (2017).
-	Rheinbay, E. The genomic landscape of advanced cancer. Nat Cancer 1, 372–373 (2020).
- Taylor, A. M. et al. Genomic and Functional Approaches to Understanding Cancer Aneuploidy. Cancer Cell 33, 676-689.e3 (2018).
- de Almeida, B. P., Vieira, A. F., Paredes, J., Bettencourt-Dias, M. & Barbosa-Morais, N. L. Pan-cancer association of a centrosome amplification gene expression signature with genomic alterations and clinical outcome. PLoS Comput Biol 15, e1006832 (2019).
-	Zehir, A. et al. Mutational landscape of metastatic cancer revealed from prospective clinical sequencing of 10,000 patients. Nat Med 23, 703–713 (2017).


# Copyright
This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.

