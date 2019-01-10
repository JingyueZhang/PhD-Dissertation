## A Novel Divide and Combine based Approach to Estimating Mixture Markov Model for Analyzing Large Categorical Time Series Data
Please see the abstract and dissertation proposal of this study for details
### Data
* Data used for case study are obtained by combining five waves of National Household Travel Survey (1990, 1995, 2001, 2009, 2017)
* Script for data preparation can be found under Data Preparation folder

### Scripts
* [DivideData.R](DivideData.R): R script for randomly split a large dataset to multiple subsamples
* [MMM_Estimation.R]([MMM_Estimation.R): R script for estimating Mixture Markov models for each subsample
* [CombineResult.R](CombineResult.R): R script for combining model estimation results from all subsamples and then apply hierarchical clustering method for obtaining final clustering solution
* [DataforPlot.R](DataforPlot.R): Prepare data for plotting activity-travel pattern
* [Plot.R](Plot.R): R script for plotting daily activity-travel patterns
* [PostAnalysis.R](PostAnalysis.R): R script for exploring influence factors of activity-travel behavior
