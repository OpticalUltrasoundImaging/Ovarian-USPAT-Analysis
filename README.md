# Ovarian-USPAT-Analysis
Algorithms for extracting quantitative imaging features from co-registered US-PAT images.  
Machine learning models for classifying ovarian lesions based on quantitative imaging features.  
by Yixiao Lin (https://opticalultrasoundimaging.wustl.edu/)

<figure>
  <img src="https://github.com/OpticalUltrasoundImaging/Ovarian-USPAT-Analysis/blob/main/System-schematic.png" alt="Imaging system">
  <figcaption>Fig 1. Co-registered US-PAT imaging system schematic</figcaption>
</figure>  
<br/>  
<figure>
  <img src="https://github.com/OpticalUltrasoundImaging/Ovarian-USPAT-Analysis/blob/main/Analysis-Pipeline.png" alt="Analysis pipeline">
  <figcaption>Fig 2. Procedure for lesion characterization and classification from a co-registered US-PAT scan.</figcaption>
</figure>

# Paper abstract
Ovarian-adnexal lesions are conventionally assessed with ultrasound (US) under the guidance of the Ovarian-Adnexal Reporting and Data System (O-RADS). However, the low specificity of O-RADS results in many unnecessary surgeries. Here, we use co-registered US and photoacoustic tomography (PAT) to improve the diagnostic accuracy of O-RADS.  Physics-based parametric algorithms for US and PAT were developed to estimate the acoustic and photoacoustic properties of 93 ovarian lesions. Additionally, statistics-based radiomic algorithms were applied to quantify differences in the lesion texture on US-PAT images. A machine learning model (US-PAT KNN model) was developed based on an optimized subset of eight US and PAT imaging features to classify a lesion as either cancer, one of four subtypes of benign lesions, or a normal ovary. The model achieved an area under the receiver operating characteristic curve (AUC) of 0.969 and a balanced six-class classification accuracy of 86.0%. For the first time , we demonstrate that the combination of parametric and radiomic US and PAT image analysis significantly improves the diagnostic accuracy in assessing ovarian-adnexal lesions.

# Usage
processOnePatientMPUS.m: <br />            &nbsp;&nbsp;&nbsp;&nbsp; function takes patient identifier as input, and computes the <ins>US</ins> imaging features from the raw data stored in the patient data folder. <br />
processOnePatientMPPA.m:  <br />           &nbsp;&nbsp;&nbsp;&nbsp; function takes patient identifier as input, and computes the <ins>PAT</ins> imaging features from the raw data stored in the patient data folder. <br />
PAT_clinical_data_cropROI.m: script <br /> &nbsp;&nbsp;&nbsp;&nbsp; allows manual cropping of the lesion on a co-registered US-PAT B scan and stores the imaging features computed inside the ROI in a struct. <br />
PAT_clinical_data_classification.m: <br /> &nbsp;&nbsp;&nbsp;&nbsp; script contains code for constructing models for classifying lesion types based on the computed imaging features. <br />
Feature-Extraction:  <br />                &nbsp;&nbsp;&nbsp;&nbsp; folder containing helper functions needed to compute imaging features.  

<figure>
  <img src="https://github.com/OpticalUltrasoundImaging/Ovarian-USPAT-Analysis/blob/main/Model-predictions-examples.png" alt="Example model predictions">
  <figcaption>Fig 3. Representative test images unused during model training, and their calibrated malignancy risk predicted by the US-PAT KNN model. The top left and bottom right corners of the lesion are marked with blue ‘+’ symbols on each scan. In (c), the stromal calcification is indicated by the pink arrow. In the plots for model prediction scores, the six numbered categories are C = cancer, BC = benign cystic, BS = benign solid, T = teratoma, E = endometriosis, N = normal.</figcaption>
</figure>

# Contact
For any questions, please contact Yixiao Lin at lin.yixiao@wustl.edu.
