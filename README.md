# Diet Mixture Model: Estimation Procedure for Maximum Likelihood Estimates of Prey Contributions to Predator Diets

Pamela Moriarty, Tim Essington and Eric Ward

###Purpose:
This is an online supplement for:
 Moriarty, P.E., T.E. Essington, E.J. Ward. 2016. A Novel Method to Estimate Prey Contributions to Predator Diets. Canadian Journal of Fishery and Aquatic Sciences. *In review.* 

This document demonstrates an example of the mixture model described in the above article. In this supplement, we provide a simulated dataset and an actual dataset used in the article above. We use the provided R code to run our mixture model in order to estimate the model parameters for the datasets. Code was built and tested under R version 3.1.0.

###File Descriptions:
- ‘README.md’: Explanation of the contents of this supplement and directions on running the model.
- ‘PredatorExample.Rdata’: Contains 1 simulated example dataset and 1 real dataset (lingcod) from Moriarty et al. (in review). 
- ‘MixtureModelExample.R’: Contains R code to run the simulated and real examples.
- Main R functions: see descriptions in table at end of this document

###Packages Needed:
plotrix  
stats4  
optimx  
numDeriv 
MASS  
fitdistrplus  

###Model Execution:
*Required Data Format*: The data must be in two columns, where one is a vector of the proportion of each stomach that is the prey type of interest (prey type i) and the other column is the corresponding total stomach contents mass for that stomach. There must be at least one stomach for each of the three cases: proportion of the prey type (p<sub>i</sub>) is 0, p<sub>i</sub> = 1, 0 < p<sub>i</sub> < 1. 

*Run Mixture Model Examples*: 
To see examples of the mixture model, open the file called ‘MixtureModelExample.R’. Once the necessary packages are installed and loaded, this file will load the example datasets. As demonstrated in the example, to run the mixture model, use the function called ‘run.model.R’. The arguments are the proportion of each stomach that is the prey type (called ‘prop’) and the total stomach contents mass for each stomach (called ‘total.mass’). Inputting these two columns of data into ‘run.model.R’ returns maximum likelihood estimates for all parameters in the mixture model. 
Additionally, to compare the estimates and error from the mixture model to other methods, use the function ‘model.comparison.R’. The arguments are the same as for ‘run.model.R’. 

###Acknowledgements:
Thanks to Anne Beaudreau for use of the lingcod data. 

###Support, Liability and Copyright: 
Support in using this software is available from Pamela Moriarty, pmoriart@uw.edu. The code and files contained in this package may be copied and distributed for non-commercial purposes. If it is used in a publication, please cite Moriarty et al. (2016) (full citation at top of this document).  The authors, University of Washington, and National Oceanic and Atmospheric Administration accept no liability for the use of this software.  

###Contact Information
Please contact Pamela Moriarty with questions and/or comments   
Email: pmoriart@uw.edu  
Mailing Address: School of Aquatic and Fishery Sciences  
                PO Box 355020  
                Seattle, WA 98195  
