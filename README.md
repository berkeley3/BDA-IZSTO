# Bayesian Data Analysis for Medical Data

This is the material for the two days course on **Bayesian Data Analysis for Medical Data** at the Istituto Zooprofilattico Sperimentale del Piemonte, Liguria e Valle d'Aosta, 13 - 14 June 2017
http://www.izsto.it/index.php/formazione/corsi-e-convegni/2025-analisi-bayesiana-dei-dati-biomedici


## Objectives

Bayesian methods are getting increasingly popular since they can allow for fitting complex models and providing richer inference from empirical observations without reference to p-values. Bayesian data analysis applies flexibly and seamlessly to complex hierarchical models and realistic data structures, including small samples, unbalanced designs, missing data, censored data while Bayesian analysis software, which is now widely available, can be used in an extremely broad variety of data models. This course will show how to carry out and interpret Bayesian statistical analysis, hands on, with free software R. The first part will be aimed to describe the basic of Bayesian inference by examining simple Bayesian models. Models that are more complicated will also be explored, including logistic regression as well as hierarchical and auto-regressive models. Bayesian computational methods, particularly Markov Chain Monte Carlo methods, will be progressively introduced as motivated by the models discussed. Emphasis will also be placed on model checking and model diagnostics.

## Course Audience

The course is intended for students, researchers from all disciplines and clinicians who want a ground-floor introduction to Bayesian data analysis for medical data.

## Program
- Introduction to Bayesian Inference
- Hierarchical Modeling
    - with examples in STAN
- Bayesian meta-analysis
- Introduction to MCMC
- Predictive Distribution Model Checking
    - missing values handling
- Bayesian disease mapping
- Prior calibration

## Course Prerequisites

No specific mathematical expertise is required. Some familiarity with statistical methods such as t-test and linear regression can be helpful, as well as some previous experience with programming in R, but they are not critical.

## Software
- R version 3.2.0 or later http://www.r-project.org/
- RStan https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
    - make sure to carefully follow all the instructions and setup correctly the toolchain 
        - for Windows, you need to install Rtools on Windows (pay attention to install the correct version according to your R version) https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows
        - for MAC check the prerequisite installation instructions https://github.com/stan-dev/rstan/wiki/RStan-Mac-OS-X-Prerequisite-Installation-Instructions
- Although it is not required, it is recommend installing RStudio https://www.rstudio.com

- R packages arm, rstan, rstanarm, glmer2stan, lme4, loo, shiny and shinystan are required
- Please note that glmer2stan is not available on CRAN. You can install it from source by using the following lines of R code 
   
   options(repos=c(getOption('repos'), glmer2stan="http://xcelab.net/R"))
   install.packages('glmer2stan', type='source')
    
- Other R packages that will be used: knitr, xtable, plotrix, ggplot2, gridExtra, LearnBayes, nleqslv, metafor, maptools, sp, spdep, RColorBrewer, classInt
 
### Contributors

* Paola Berchialla

<!----_Material is under development and subject to change._----> 
