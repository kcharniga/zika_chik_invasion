# zika_chik_invasion

The models folder in this repository contains the data files and models necessary to reproduce the main analysis in our paper "Spatial and temporal invasion dynamics of the 2014-2017 Zika and chikungunya epidemics in Colombia." The figures folder contains the data and code necessary to reproduce the figures in the main text. The animations folder contains the data and code necessary to reproduce the animated hexagon maps (S1 and S2 Movies). For a description of the datasets, see the bottom of this README file.

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
* **d_cities**: numeric square matrix with N x N elements representing the geographic distance between two cities in km. Matrix is symmetric, matrix where rows and columns are cities. Values are the geographic distance between each city in km. Zeroes, which represent the distance between one city and itself, on the diagonal have been replaced by NA to avoid dividing by 0 in the calculation of the distance kernel. 
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

## Explanation of datasets

### Data for models (main text)
* **S1: chik_d_cities.txt** 
Matrix where rows and columns are cities. Values are the geographic distance between each city in km. Zeroes, which represent the distance between one city and itself, on the diagonal have been replaced by NA to avoid dividing by 0 in the calculation of the distance kernel.
* **S2: chik_I.txt**
Matrix where rows are weeks and columns are cities. The matrix element in row i and column j equals 0 if the city corresponding to column j is susceptible at time i, and equals 1 if the city corresponding to column j is invaded at time i. This matrix was generated from invasion time. Week of invasion was determined for each city based on when the first cases were reported.
* **S3: chik_infectious.txt**
Matrix with weeks as rows and cities as columns where the matrix element in row i and column j equals 0 if the city corresponding to column j is susceptible at time i, and equals 1 if the city corresponding to column j is infectious at time i. This matrix is based off of the I matrix: all of the values are shifted down by one row to represent a one-week time lag from being invaded to becoming infectious.
* **S4: chik_J.txt**
Matrix with weeks as rows and cities as columns. This is based off of the infectious matrix. J equals NA when infectious equals 1 and J equals 1 when infectious equals 0. A value of 1 for J indicates that a city contributes to the likelihood during that week. The first invaded cities have NA for all time points and never contribute to the likelihood.
* **S5: chik_N.txt**
Vector with population sizes for each city. Data comes from DANE 2016 projections derived from the 2005 Colombia Census. Population sizes were divided by 10,000.
* **S6: chik_H.txt**
Matrix with weeks as rows and cities as columns. Values are the weekly case counts weighted by each infection’s generation time.
* **S7: zika_d_cities.txt**
Matrix where rows and columns are cities. Values are the geographic distance between each city in km. Zeroes, which represent the distance between one city and itself, on the diagonal have been replaced by NA to avoid dividing by 0 in the calculation of the distance kernel.
* **S8: zika_I.txt**
Matrix where rows are weeks and columns are cities. The matrix element in row i and column j equals 0 if the city corresponding to column j is susceptible at time i, and equals 1 if the city corresponding to column j is invaded at time i. This matrix was generated from invasion time. Week of invasion was determined for each city based on when the first cases were reported.
* **S9: zika_infectious.txt**
Matrix with weeks as rows and cities as columns where the matrix element in row i and column j equals 0 if the city corresponding to column j is susceptible at time i, and equals 1 if the city corresponding to column j is infectious at time i. This matrix is based off of the I matrix: all of the values are shifted down by one row to represent a one-week time lag from being invaded to becoming infectious.
* **S10: zika_J.txt**
Matrix with weeks as rows and cities as columns. This is based off of the infectious matrix. J equals NA when infectious equals 1 and J equals 1 when infectious equals 0. A value of 1 for J indicates that a city contributes to the likelihood during that week. The first invaded cities have NA for all time points and never contribute to the likelihood.
* **S11: zika_N.txt**
Vector with population sizes for each city. Data comes from DANE 2016 projections derived from the 2005 Colombia Census. Population sizes were divided by 10,000.
* **S12: zika_H.txt**
Matrix with weeks as rows and cities as columns. Values are the weekly case counts weighted by each infection’s generation time.

### Data for figures (main text)
* **S1-S12** (see above for details)
* **S13: lat_lon_complete.txt**
Latitudes and longitudes corresponding to the geographic center of each city in Colombia.
7 variables,
admin2: city code,
admin1_code: department code,
admin1_name: department name,
admin2_name: city name,
admin2_longitude: longitude of city,
admin2_latitude: latitude of city,
pop: 2016 population projections from DANE
* **S14: chik_first_reported_cases_338.txt**
invasion weeks for 338 cities based on first reported cases of CHIKV in each city.
admin2: city code,
first_report: week before cases were first reported
* **S15: zika_first_reported_cases_288.txt**
invasion weeks for 288 cities based on first reported cases of ZIKV in each city.
admin2: city code,
first_report: week before cases were first reported
* **S16: distance_matrix_km.txt**
Geographic distance between all cities in km.
rows: 1,122 cities and 
columns: 1,122 cities
* **S17: chik_admin1_series.txt**
Time series of reported CHIKV cases at admin level 1 consisting of 6 variables,
disease: CHIKV,
admin1: department code,
week: epidemiological week,
year: year,
cases: number of reported cases,
admin1_name: department name
* **S18: zika_admin1_series.txt**
Time series of reported ZIKV cases at admin level 1 consisting of 6 variables,
disease: ZIKV,
admin1: department code,
week: epidemiological week,
year: year,
cases: number of reported cases,
admin1_name: department name
* **S19: admin1_names.txt**
Names of departments and codes consisting of 3 variables,
admin1_code: department code,
admin1_name: department name as it is in disease data,
admin1_map: department name with correct Spanish accents
* **S20: pop_wt_centroids.txt**
Centroids and population weighted centroids for all countries. 9 variables,
V1: country (admin 0),
V2: department or state (admin 1),
V3: country id,
V4: admin 1 id,
V5: longitude of admin 1 centroid,
V6: latitude of admin 1 centroid,
V7: longitude of admin 1 population weighted centroid,
V8: latitude of admin 1 population weighted centroid,
V9: population of admin 1
* **S21: chik_mcmc_chain1.txt**
Posterior distribution of parameters for MCMC chains from best-fitting CHIKV model
* **S22: zika_mcmc_chain1.txt**
Posterior distribution of parameters for MCMC chains from best-fitting ZIKV model
* **S23: chik_I_100.txt**
Similar to S2 except with more rows added to the bottom of the matrix in order to run the epidemic simulations for more weeks.
* **S24: chik_infectious_100.txt**
Similar to S3 except with more rows added to the bottom of the matrix in order to run the epidemic simulations for more weeks.
* **S25: chik_J_100.txt** 
Similar to S4 except with more rows added to the bottom of the matrix in order to run the epidemic simulations for more weeks.
* **S26: chik_H_100.txt**
Similar to S6 except with more rows added to the bottom of the matrix in order to run the epidemic simulations for more weeks.
* **S27: zika_I_50.txt**
Similar to S8 except with more rows added to the bottom of the matrix in order to run the epidemic simulations for more weeks.
* **S28: zika_infectious_50.txt**
Similar to S9 except with more rows added to the bottom of the matrix in order to run the epidemic simulations for more weeks.
* **S29: zika_J_50.txt** 
Similar to S10 except with more rows added to the bottom of the matrix in order to run the epidemic simulations for more weeks.
* **S30: zika_H_50.txt**
Similar to S12 except with more rows added to the bottom of the matrix in order to run the epidemic simulations for more weeks.

### Data for animations (supplement)
* **S31: chik_monthly_inc_for_animation.txt**
10 variables,
month_no: number of months (time points) for which we will run animation,
cases: number of reported cases in each department, 
admin1_name: department name according to DANE classification,
admin1_code: administrative code associated to each department,
dep_short: abbreviated department name,
pop_2016_DANE: population projections for 2016 from DANE for each department,
month_num: calendar month number,
month_name: calendar month name,
year: year,
inc_rate: number of cases divided by population of department times 100,000
* **S32: zika_monthly_inc_for_animation.txt**
10 variables,
month_no: number of months (time points) for which we will run animation,
cases: number of reported cases in each department, 
admin1_name: department name according to DANE classification,
admin1_code: administrative code associated to each department,
dep_short: abbreviated department name,
pop_2016_DANE: population projections for 2016 from DANE for each department,
month_num: calendar month number,
month_name: calendar month name,
year: year,
inc_rate: number of cases divided by population of department times 100,000
