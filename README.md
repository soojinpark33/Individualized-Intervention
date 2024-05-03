# Individualized_Intervention

Soojin Park<sup>1</sup>, Namhwa Lee<sup>2</sup>, and Rafael Quintana<sup>3</sup>

<sup>1</sup> School of Education, University of California, Riverside  
<sup>2</sup> Department of Statistics, University of California, Riverside
<sup>3</sup> School of Education and Human Sciences, The University of Kansas


## Overview

Causal decomposition analysis aims to identify risk factors (referred to as ‘mediators’) that contribute to social disparities in an outcome. Despite promising developments in causal decomposition analysis, current methods are limited to addressing a time-fixed mediator and outcome only, which has restricted our understanding of the causal mechanisms underlying social disparities. In particular, existing approaches largely overlook individual characteristics when designing (hypothetical) interventions to reduce disparities. To address this issue, we extend current longitudinal mediation approaches to the context of disparities research. Specifically, we develop a novel decomposition analysis method that addresses individual characteristics by (1) using optimal dynamic treatment regimes (DTRs) and (2) conditioning on a selective set of individual characteristics. Incorporating optimal DTRs into the design of interventions can be used to strike a balance between equity (reducing disparities) and excellence (improving individuals' outcomes). We illustrate the proposed method using the High School Longitudinal Study data.

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
