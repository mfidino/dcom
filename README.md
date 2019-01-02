# dcom

This is the data and code for:

Fidino, M., Simonis, J. L., and Magle, S. B. (2018). A multi-state dynamic occupancy model to estimate local colonization-extinction rates and patterns of co-occurrence between two or more interacting species. *Methods in Ecology and Evolution*, 0:1--12.

This repository consists of 4 scripts used for analysis, model summary, and plotting. They are located within the working directory. These scripts include:

**fit_softmax_model.R:** This script  1) reads in the raw detection non-detection data for coyote, opossum, and raccoon and converts it to community states at each site that can be supplied to the dynamic co-occurrence model, 2) reads in the covariate data and generates the urbanization covariate URB, 3) fits the 6 different models outlined in the manuscript, and 4) compares the fit of each of these models.

**fit_models_utility_functions.R:** This script contains a number of useful functions that are sourced for fit_softmax_model.R. The functions are commented out within this script and describe what they do. It include three functions: two to generate initial values for JAGS models (those with and without species interactions) and one to generate the vectors to organize the species interactions within a matrix in JAGS.

**calculate_cpo.R:** This function calculates the conditional predictive ordinate for each data point within a model when given the model output as a matrix and the data list supplied to JAGS. It returns the summary statistic:〖-Σ〗_(k,t)  log⁡(〖CPO〗_kt ) for site k and time t, where lower values indicate better model fit. It is sourced within fit_softmax_model.R
plotting_script.R: This is the code used to generate the figures in the manuscript. Some additional calls to imageMagick were done in the console to generate the final figures.

**plotting_script_utility_functions.R:** This script included functions to help generate the figures in plotting_script.R. It includes 4 functions: one to generate a species colonization and extinction rate given the presence or absence of another species, two to label the x and y axes for the figures, and one to generate each sub-figure for figure 3 (the occupancy plot).
calculate_steady_state.R: This function calculates the expected occupancy rate along an environmental gradient for the urbanization covariate with the best fit model and the data list used to fit this model. It is sourced and called within plotting_script.R to generate figure 3.



**Data S1 also has 2 JAGS dynamic co-occurrence occupancy models that are used in this analysis. They should be placed within the jags_models sub-folder of the working directory. These include:**

**dcom_inxs.R:** This is the dynamic co-occurrence occupancy model that is used for models 1 through 4 and includes interactions between species. It is general enough such that any number of covariates could be fit to each process (e.g., include 2 covariates on all species interactions for colonization). The parameters within each linear predictor are named as they are within the manuscript. 

**dcom_no_inxs.R:** This is the dynamic co-occurrence occupancy model that is used for models 5 and 6 and does not include interactions between species. It is general enough such that any number of covariates could be fit to each process. The parameters within each linear predictor are named as they are within the manuscript. 

**There are two data files within the data sub-folder which are used in this analysis. They include:**

**fidino_sp_data.csv:** This is the detection / non-detection data for coyote, opossum, and raccoon in long format. It is 44772 records and 9 columns. The ordering of this data frame is an important part of how the data gets sorted and summarized within fit_softmax_model.R. In R, these data should be read such that character strings are not automatically turned into factors:

`data <- read.csv(“fidino_sp_data.csv”, header = TRUE, stringsAsFactors = FALSE)`

## Structure of `fidino_sp_data.csv`

| Column header  | Data type  | Description  |
|---|---|---|
| Season  | Categorical  |This contains information on which season and year the detection / non-detection data is associated it. Each element has four characters. The first two characters describe which season (SP = spring, SU = summer, FA = fall, WI = winter) while the final two characters describe the year (10 = 2010, 11 = 2011, etc.). The data is currently sorted so that the season and year continue temporally (i.e., SP10 data, SU10 data, FA10 data, WI11  data, etc.)   |
|  Week | Categorical  | Which week of sampling do the 7 detection / non-detection days fall within per season. Represented as ‘week 1’, ‘week 2’, ‘week 3’ and ‘week 4’.  |
| Date  | Date  | The day associated to the detection / non-detection data. Format is `%Y-%M-%D`  |

