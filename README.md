# PhD-Dissertation
This repository includes Jingyue Zhang's Ph.D. dissertation entitled "Advanced Modeling and Forecasting Methodologies for Activity-travel Behavior Analysis". The dissertation includes three studies, abstract and scripts will be added for each of them
## Dissertation documents include:
* [DissertationProposal_Jingyue.pdf](DissertationProposal_Jingyue.pdf)
* [DivideandCombineStudy_Abstract.txt](DivideandCombineStudy_Abstract.txt)

## Scripts:
### Data Preparation
*  [NHTS2017_DataPreparation.ipynb](NHTS2017_DataPreparation.ipynb)
    * Download National Household Travel Survey Data at https://nhts.ornl.gov/
    * Perform data cleanse to remove the records with invalid trip purpose; trip start and end time, travel mode etc. 
    * Characterize individual's daily activity-travel behavior as categorical time series data and create new data files
* [CHTS2010_TourFormation.ipynb](CHTS2010_TourFormation.ipynb)
    * Data used for this script is 2010 California Household Travel Survey https://www.nrel.gov/transportation/secure-transportation-data/tsdc-california-travel-survey.html
### Divide and Combine Study
* [DivideData.R](DivideData.R): R script for randomly split a large dataset to multiple subsamples
* [MMM_Estimation.R]([MMM_Estimation.R): R script for estimating Mixture Markov models for each subsample
* [CombineResult.R](CombineResult.R): R script for combining model estimation results from all subsamples and then apply hierarchical clustering method for obtaining final clustering solution
* [DataforPlot.R](DataforPlot.R): Prepare data for plotting activity-travel pattern
* [Plot.R](Plot.R): R script for plotting daily activity-travel patterns
* [PostAnalysis.R](PostAnalysis.R): R script for exploring influence factors of activity-travel behavior
### Time Varying Mixture Markov Model
 * [TimeVaryingMMM_Estimation.R](TimeVaryingMMM_Estimation.R): EM algorithm is implemented in R for estimated time-varying mixture Markov model
* [MarkovChainSimulation.R](MarkovChainSimulation.R): R script for simulating sequential data based on predefined time-varying transition probabilities
* [DescriptiveAnalysis.R](DescriptiveAnalysis.R): R script for conducting descriptive analysis of activity-travel behavior of the elderly population
  * Plot activity-travel pattern
  * Plot activity-travel patterns of working and non-working; younger (65 to 75 years old) and older (above 75 years) elderly
  * Plot count of transition pairs over time

### Contact
* Email: jingyue.zhang@uconn.edu

### License
* This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
