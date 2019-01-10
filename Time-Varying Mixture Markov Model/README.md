# Time-Varying Mixture Markov Model
This folder includes scripts for estimating a time-varying mixture Markov model. Model formulation can be found in DissertationProposal_Jingyue.pdf.

## Data
Data used for this study is obtained by combining multiyear National Household Travel Survey. Script for data preparation can be found in Data Preparation folder

## Script
* [TimeVaryingMMM_Estimation.R](TimeVaryingMMM_Estimation.R): EM algorithm is implemented in R for estimated time-varying mixture Markov model
* [MarkovChainSimulation.R](MarkovChainSimulation.R): R script for simulating sequential data based on predefined time-varying transition probabilities
* [DescriptiveAnalysis.R](DescriptiveAnalysis.R): R script for conducting descriptive analysis of activity-travel behavior of the elderly population
  * Plot activity-travel pattern
  * Plot activity-travel patterns of working and non-working; younger (65 to 75 years old) and older (above 75 years) elderly
  * Plot count of transition pairs over time
