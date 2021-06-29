# Genomic events shaping epithelial-to-mesenchymal trajectrories in cancer [under construction]
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
- **run_pipepline_ml.sh:** The bash script call several Rscript to perform the identification of the biomarkers. Please open the .sh and infofile.txt to get more details about the procedure.


## tissue_specific_trajectories
This folder contains the pipeline used to identify tissue specifics trajectories in LUAD and BRCA, and the code to define the markers associated with EMT states

## clinical_EMT
This folder contains the code to perform the characterization of the clinical features related with EMT. The folder contain also the code to analyze pharmacogenomic datasets (GDSC, depmap, POG570)

# Copyright

This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.

