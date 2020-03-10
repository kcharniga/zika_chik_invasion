# zika_chik_invasion

The models folder in this repository contains the data files and models necessary to reproduce the main analysis in our paper "Spatial and temporal invasion dynamics of the 2014-2017 Zika and chikungunya epidemics in Colombia." See details on the functions and inputs below.

The figures folder contains the data and code necessary to reproduce the figures in the main text. For a description of the datasets, see the supplementary information.

The animations folder contains the data and code necessary to reproduce the animated hexagon maps in the supplement. For a description of the datasets, see the supplementary information.

## What does it do?

The model functions include:
* **gravity_model**: implements a gravity model using Metropolis-Hastings Markov chain Monte Carlo 
* **get_params**: extracts the posterior distribution of parameters, log likelihood, and acceptances for each of three chains
* **calc_accept**: calculates acceptance rates of parameters
* **calc_dic**: calculates the Deviance Information Criteria (DIC) for each of three MCMC chains
* **calc_med_est**: calculates median parameter estimates 95% credible intervals 
* **plot_res**: plots histograms of the posterior distributions of parameters and MCMC traces

## Data inputs for models

The following data are provided to run the gravity models:
data1, list type containing:
* **d_cities**: numeric square matrix with N x N elements representing the great circle distance between two cities in km. Matrix is symmetric, matrix where rows and columns are cities. Values are the great circle distance between each city in km. Zeroes, which represent the distance between one city and itself, on the diagonal have been replaced by NA to avoid dividing by 0 in the calculation of the distance kernel. 
* **I**: integer matrix with weeks as rows and cities as columns where the matrix element in row i and column j equals 0 if the city corresponding to column j is susceptible at time i, and equals 1 if the city corresponding to column j is invaded at time i. This matrix was generated from estimates of invasion time. Week of invasion was estimated for each city based on when the first cases were reported.
* **infectious**: integer matrix with weeks as rows and cities as columns where the matrix element in row i and column j equals 0 if the city corresponding to column j is susceptible at time i, and equals 1 if the city corresponding to column j is infectious at time i. This matrix is based off of the I matrix. All of the values are shifted down by one row to represent a one week time lag from being invaded to becoming infectious. 
* **J**: integer matrix with weeks as rows and cities as columns. This is based off of the infectious matrix. J equals NA when infectious equals 1 and J equals 1 when infectious equals 0. A value of 1 for J indicates that a city contributes to the likelihood during that week. The first infected cities have NA for all time points and never contribute to the likelihood
* **N**: vector with population sizes for each city. Data comes from DANE 2016 estimates derived from the 2005 Colombia Census. Population sizes were divided by 10,000.
* **H**: matrix with weeks as rows and cities as columns. Values are the weekly case counts weighted by each infection’s generation time.

## Model specifications

The following model options must be specified:
* **iterations**: positive integer number of samples for MCMC
* **burnIn**: non-negative integer number of iterations that will be discarded from the final Metropolis-Hastings MCMC routine, referred to as the burn-in period
* **number_of_parameters**: the number of parameters to be estimated by MCMC
* **vals**: list containing starting values and standard deviations of the proposal function for each of 3 MCMC chains
