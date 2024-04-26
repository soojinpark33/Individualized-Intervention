# Individualized_Intervention

Soojin Park<sup>1</sup>, Namhwa Lee<sup>2</sup>, and Rafael Quintana<sup>3</sup>

<sup>1</sup> School of Education, University of California, Riverside  
<sup>2</sup> Department of Statistics, University of California, Riverside
<sup>3</sup> School of Education and Human Sciences, The University of Kansas


## Overview

A key objective of decomposition analysis is to identify a factor (the `mediator') contributing to disparities in an outcome between social groups. In decomposition analysis, a scholarly interest often centers on estimating how much the disparity (e.g., health disparities between Black women and White men) would be reduced/remain if we set the mediator (e.g., education) distribution of one social group equal to another. However, causally identifying disparity reduction and remaining depends on the no omitted mediator-outcome confounding assumption, which is not empirically testable. Therefore, we propose a set of sensitivity analyses to assess the robustness of disparity reduction to possible unobserved confounding. We derived general bias formulas for disparity reduction, which can be used beyond a particular statistical model and do not require any functional assumptions. Moreover, the same bias formulas apply with unobserved confounding measured before and after the group status. Based on the formulas, we provide sensitivity analysis techniques based on regression coefficients and $R^2$ values by extending the existing approaches. The $R^2$-based sensitivity analysis offers a straightforward interpretation of sensitivity parameters and a standard way to report the robustness of research findings. Although we introduce sensitivity analysis techniques in the context of decomposition analysis, they can be utilized in any mediation setting based on interventional indirect effects when the exposure is randomized (or conditionally ignorable given covariates).

For more details of our proposed methods, see [our paper](XX). 
Here, we provide `R` codes to reproduce our simulation study and replicate our data analysis. 

## Case Study

* Data
  
The HSLS:09 data used for the case study can be downloaded from the MIDUS portal by clicking [here](https://nces.ed.gov/surveys/hsls09/hsls09_data.asp). 

* `individualized_HSLS.R` 
 
   This `R` file replicates Tables 1 and 2 of our study.

## Simulation Study

* `Pop_simulation.R`  

   This `R` file contains the simulation codes for evaluating the performance of the proposed estimators for each individualized effect. This code replicates our results in Appendix D of our paper.


These supplementary materials are provided solely for the purpose of reproducibility and must be used in compliance with academic ethical guidelines. If you reference these materials in your own work, please ensure proper citation of the original sources.
