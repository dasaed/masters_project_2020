# Exploring Interpretability Models for Machine Learning in Prioritizing Bona Fide Bacterial Small RNAs Supplementary Material
## This project was created in partial fulfillment of the requirements for the degree of Master of Science in Computer Science at Memorial University of Newfoundland.
### Author: Carlos Dasaed Salcedo Carreno
### Supervisor: Dr. Lourdes Pena-Castillo, PhD. 

This project is based of the work found [here](https://github.com/BioinformaticsLabAtMUN/sRNARanking)(ie. this is where we obtained the RF\_classifier4sRNA.rds file), and the corresponding article can be found in [Eppenhof EJ, Peña-Castillo L. (2019) Prioritizing bona fide bacterial small RNAs with machine learning classifiers. PeerJ 7:e6304](https://peerj.com/articles/6304/)

## Abstract - as is in my report
While there is a lot of research focused on creating black box machine learning (ML)models with high performance metrics, understanding what these models are learningor how they are arriving at their predictions has received less attention until recently.However, understanding what these black box models are actually learning can be justas important as the model’s performance.  Here, we applied Interpretability Models(IM)  to  the  random  forest  (RF)  model  developed  by  Eppenhof  et  al.   [1]  in  thepursuit of new knowledge about bona fide bacterial small RNAs, or sRNAs, and theusefulness of IMs.  Sequencing-based studies have generated an overwhelming amountof  putative  sRNAs  available  in  the  literature.   Eppenhof  et  al.’s  model  is  designedto  assist  researchers  in  narrowing  down  putative  sRNAs  by  prioritizing  sequencesbased on the predicted likelihood of them being a bona fide sRNA. We used PartialDependence Plots (PDP), Local Interpretable Model-agnostic Explanations (LIME),and SHapley Additive Explanations (SHAP) to obtain robust results and assess theagreement between the IMs.  In particular, the PDPs showed that the lower the freeenergy  of  the  predicted  secondary  structure  of  a  genomic  sequence,  the  higher  thelikelihood of it being bona fide sRNA. The distance to the opening reading frames(ORF) was also identified as one of the most important features in the identificationof  bona  fide  sRNAs.   While  LIME  and  SHAP  served  to  confirm  findings  from  thePDPs, they may have also hinted at the existence of a possible bias in the RF modelwhen the distance to the ORFs is within a couple of nucleotides.  Thanks to the IMs,not only was Eppenhof et al.’s model explained, but possible improvements may havealso been identified, which may lead to new models with better performance metrics


## About the code

The code is written in R, and I recommend opening up the files in R Studio to run the code, rather than the terminal, as the main file is written as an R notebook (IM\_sRAN.Rmd) file rather than a normal R script.
For this project, 3 different IMs were used: PDPs, LIME, and SHAP. 
However, while I did try to keep all the code in a single notebooks, there were conflicts between the different IMs' libraries, and I was forced to split them. 
  * The *RFE_sRNA.R* contains the code using the Random Forest Explainer library and the PDP R libray. (ie, this is the code that generated my PDPs)
  * The *shapper\_sRNA.R* library contains the code generated the SHAP explanations table. This code in particular generates a csv.  
  * The main file is written as an R Notebook, and it contains the output from the other files. 

All the R files contain comments, which may serve one of following purposes: 
  1. Provide explanations about the code or the reasoning behind certain functions
  1. To prevent certain lines from being run. For example, every time LIME runs, it generates different explanations, and so to prevent the lime\_explanations.pdf from being overwritten at every run, the corresponding lines were commented out.

I recommend running the code from R studio in the following order:
  1. RFE\_sRNA.R
  2. shapper\_sRNA.R
  3. IM\_sRNA.Rmd

The lines that install the libraries have been commented out, so if you're missing libraries, I recommend modifying the corresponding file and uncommenting the lines so that the libraries are installed. 



