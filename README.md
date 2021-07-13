# Epithelial-to-mesenchymal quantification in cancer
### Guidantonio Malagoli Tagliazucchi, Maria Secrier 
UCL Genetics Institute, Department of Genetics,  Evolution and Environment,  University College London,  UK 

<p align="center">
  <img width="700" height="400" src="https://github.com/secrierlab/EMT/blob/main/figrepo.png">
</p>

# Table of contents

## ReconstructionEMTbulk
This folder contains the code to reconstruct the EMT trajectory of the bulk RNA-seq data from TCGA
- **namescript1.R:**
- **namescript2.R:**

## HMM_macrostates_EMT
This folder contains the code to perform the segmentation of the EMT trajectory and to identify the macrostates of EMT
- **run_hmm.R:** This script contain the code to run the segmentation and perform the identification of the EMT states
 
## extrinsic_hallmarks
This folder contains the code to quantify the extrinsic hallmarks of EMT
- **run_fibroblast_deconvolution.R:** This is the script used to quantify the levels of fibroblast infiltration, contains also the code to run the multinomial logistic regression with the fibroblasts signals and the HMM states.
- **run_tme_deconvolution.R:** This is the script used to quantify the cell components of Tumour Micro Environment (TME), contains also the code to run the multinomial logistic regression with the TME signals and the HMM states.
- **run_tme_deconvolution.R:** This is the script used to quantify the cell components of Tumour Micro Environment (TME), contains also the code to run the multinomial logistic regression with the TME signals and the HMM states.
- **characterize_extrinsic_hallmarks.R:** In this script there is the procedure to quantify the hypoxia scores (according to Buffa et al.). This script create also a data.frame that combine the quantification of several others hallmarks (e.g. centromeric amplification (CA20), aneuploidy, etc)

## intrinsic_hallmarks
This folder contains the code to quantify the intrinsic hallmarks of EMT
- **mutational_signatures.R:** Script to perform the identification of the mutational signatures in each EMT states with mixed-effect models (Buffa et al.)

## biomarkers_EMT
This folder contains the pipeline used to identify the biomarkers of EMT
- **run_pipepline_ml.sh:** This script contains the pipeline to perform the identification of the biomarkers. Please open the .sh and infofile.txt to get more details about the procedure.


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
- **validation_Schaller_MUT.R:** This script check if putative biomarkers of EMT associated with CNV are targets of TFs important for EMT

# How to cite

# Copyright
This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.

