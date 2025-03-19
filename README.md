***AutScore.r***

- *AutScore.r* is a statistical model generated prediction score to prioritize the ASD candidate rare variants (SNVs and INDELS) based on the Whole Exome Sequencing (WES) data from simplex/multiplex ASD affected families from the Israeli population.
  
- After preparing the list of rare LP/P/LGD variants with annotations from different tools (e.g., InterVar, Psi-Variant, etc.) and databases (SFARI, DOMINO, DisGeNET, ClinVar, etc.) as described elsewhere, AutScore.r can be implemented by a dedicated R script: *AutScore.r_implementation.R*. Here we transform the modules of AutScore (I, P, D, etc.) first and then fit a logit model and generate the prediction score, defined as *AutScore.r*.

- In this reporsitory we have provide a sample dataset to annotate *AutScore.r* using the script: *AutScore.r_implementation.R*. 
  
- Soon we will publish a 2nd version of *AutScore.r* using which one can QC and annotate all the required annotations to .vcf file and compute such prediction scores based on AI models to prioritize ASD candidate variants. 

Please refer to our published paper on AutScore to learn more about specific strategies and modules used. For any further query, please contact: apurba.shil316@gmail.com or idanmen@bgu.ac.il (https://www.idanme.com/).
