####################################
#
# fit_softmax_model
#
# Code written by Mason Fidino
#
####################################

# This script prepares the raw detection / non-detection data for analysis
# and fits varying models to these data. Most of the code is data manipulation.
# Hopefuuly, this should provide those who are interested in fitting this
# model to their own data some tips on how to get it ready for analysis.

# folders required:
#  data: stores 'fidino_sp_data.csv' and 'fidino_covariate_data.csv'
#  jags_models: stores the two dcom models used in this analysis
#  model_output: where the mcmc output for each jags model is stored as an RDS
#  file

# load packages
library(dplyr)
library(LaplacesDemon)
library(runjags)
library(coda)
library(rjags)

###########################################################
# Note: snake_case for objects, CamelCase for column names.
###########################################################

# Source utility functions for fitting these models.
# the functions that supply initial values to the jags model
# require you to have the object 'data_list' in your environment,
# which is the list of data supplied to JAGS, as well as the object
# 'sp_inits' which is the maximum observed community state that 
# the community was observed within.
source("./fit_models_utility_functions.R")
### This file contains 3 functions:
#
### inits: Used to supply initial values to models with species interactions
#
### inits_no_inxs: Used to supply initial values to models without species \
###                interactions
#
### make_inxs: We store interactive parameters in a matrix within the model.
###            This function creates vectors that do the nested indexing
###            within the model. Give it an integer of the number of species
###            and it will make a vector for the row and column indexing.

# Source function to calculate the cpo score for each model
source("./calculate_cpo.R")
### This file contains 1 function:
#
### calculate_cpo: This requires the jags output of all the parameters of the
###                model as the first object, which needs to be a matrix. The
###                second object is the data_list supplied to that specific
###                model. It calculates the likelihood of each data point
###                per step in the MCMC chain and from there calculates
###                the cpo score.

# bring in the raw species data
raw_dat <- read.csv( "./data/fidino_sp_data.csv", 
                     header = TRUE, stringsAsFactors = FALSE )
# make seasonal column a factor
raw_dat$Season <- factor(raw_dat$Season, levels = unique(raw_dat$Season))

# We are going to reduce down to weekly detections instead of daily
# First step: Determine if week is entirely NA
# Returns new column, DaysSampled, which takes a value between 0 and 7 for 
# number of days a camera trap was active during theat week
days_per_week_sampled <- raw_dat %>% group_by(Season, SurveyID, SeasonWeek) %>% 
  summarise(DaysSampled = length(Coyote) - sum(is.na(Coyote)),
            StationID = unique(StationID))

# Condense the daily species detection data to weekly 
week_dat <- raw_dat %>% group_by(Season, SurveyID, SeasonWeek) %>% 
  summarise(Coyote = sum(Coyote, na.rm = TRUE),
            Opossum = sum(Opossum, na.rm = TRUE),
            Raccoon = sum(Raccoon, na.rm = TRUE),
            StationID = unique(StationID))

# If days_per_week_sampled$DaysSampled is zero, make species columns NA
# in week_dat because we do not know if the species was there or not
week_dat$Coyote[days_per_week_sampled$DaysSampled== 0] <- NA
week_dat$Opossum[days_per_week_sampled$DaysSampled== 0] <- NA
week_dat$Raccoon[days_per_week_sampled$DaysSampled== 0] <- NA

# Make the species columns binary, but create a different object of week_dat
# to do this. 
wd_bin <- week_dat
wd_bin$Coyote[wd_bin$Coyote>0] <- 1
wd_bin$Opossum[wd_bin$Opossum>0] <- 1
wd_bin$Raccoon[wd_bin$Raccoon>0] <- 1

# Combine the coyote, opossum, and raccoon states. This will help us determine
# the observed community state. For example, 100 would mean coyote is there,
# opossum are not, and raccoon are not given the ordering.
wd_bin$State <- paste0( wd_bin$Coyote,
                        wd_bin$Opossum,
                        wd_bin$Raccoon )

# convert the binary states to numeric states, but make the NA
# values remain NA. This data.frame is the conversion between the binary
# switches and the numeric categories we will use with the
# categorical distribution.
categories <- data.frame(State = c( "000", "100", "010", 
                                    "001", "110", "011",
                                    "101", "111", "NANANA" ),
                         Lvls = c( seq( 1:8 ), NA ) )

# combine lvls with dat with a left join so now
# there is a new column 'lvls' that reflect the categories in
# dat$state
wd_bin <- dplyr::left_join(wd_bin, categories, by = "State")

# get the stationIDs that are all NA. These are sites
# that have absolutely 0 data (i.e., sites introduced
# after the dataset used in this analysis). 
NA_all <- wd_bin %>% group_by( StationID ) %>% 
  summarise( missing = sum( is.na( Coyote ) ) ) %>% data.frame( ) %>% 
  filter( missing == 52 ) # 52 = 4 weeks * 13 sampling seasons

# remove sites that have more than 11 missing time steps (i.e., seasons)
NA_survey <- wd_bin %>% group_by( SurveyID, StationID ) %>% 
  summarise( Missing = sum( is.na( Coyote ) ) ) %>% 
  mutate( NGone = sum( Missing == 4 ) ) %>% # equals 1 if 4 weeks of no sampling
  group_by( StationID ) %>% 
  summarise( NSeasonGone =sum( NGone )) %>% 
  filter( NSeasonGone > 11 )

# There are just the sampling periods of sites that wherein sampling
# did not occur
NA_missing <- wd_bin %>% group_by( SurveyID, StationID ) %>% 
  summarise( Missing = sum( is.na( Coyote ) ) ) %>% 
  mutate( NGone = sum( Missing == 4 ) ) %>% 
  summarise( NSeasonGone =sum( NGone ) )

# remove these from the data, subset on rows
wd_bin <- wd_bin[-which( wd_bin$StationID %in% NA_all$StationID ), ]
wd_bin <- wd_bin[-which( wd_bin$StationID %in% NA_survey$StationID ), ]

# do the same thing for days_per_week_sampled (removing station IDs)
days_per_week_sampled <- 
  days_per_week_sampled[-which(days_per_week_sampled$StationID %in% 
                               NA_all$StationID), ]
days_per_week_sampled <- 
  days_per_week_sampled[-which(days_per_week_sampled$StationID %in% 
                                 NA_survey$StationID), ]
# Compile some info to convert lvls into a 3d array (weeks x nsite x nyear)
# Number of unique sites
usi <- length(unique(wd_bin$StationID))
# Number of unique seasons
uyr <- length(unique(wd_bin$Season))
# Number of weeks each site was observed
udays <- max(table(wd_bin$SurveyID))

# convert the levels to a 3 dimensional array
sp_ob_state <- array( wd_bin$Lvls, dim = c( udays, usi, uyr ) )
# have dimensions go nsite x nyear x n_observations
sp_ob_state <- aperm( sp_ob_state, c( 2,3,1 ) )

# make a days per week sampled array as well for the detection level
days_samp <- array( days_per_week_sampled$DaysSampled, dim = c( udays, usi, uyr ) )
# give it the same order as sp_ob_state
days_samp <- aperm( days_samp, c( 2, 3, 1 ) )

# we need to determine the maximum state that a site was in each survey
# to supply as an initial value in a JAGS model. Conversely, you could 
# just make the maxiumum state being the starting value. For example,
# if we have 3 observations that go spA, spA, spB then we know that
# we are at least in state AB (assuming closure).

# make SurveyID a factor to keep ordering correct
wd_bin$SurveyID <- factor( wd_bin$SurveyID, levels = unique( wd_bin$SurveyID ) )

# grouping by SurveyID, which denotes the (potential) 4 weeks of observation
# at a site on a given season. This is basically the total number of
# weeks each species was detected per site per season
dat_sum <- wd_bin %>% group_by( SurveyID ) %>% 
  summarise( Coyote = sum( Coyote, na.rm = TRUE ),
             Opossum = sum( Opossum, na.rm = TRUE ),
             Raccoon = sum( Raccoon, na.rm = TRUE ) )

# turn counts over 1 to 1
dat_sum$Coyote[ dat_sum$Coyote > 0 ] <- 1
dat_sum$Opossum[ dat_sum$Opossum > 0 ] <- 1
dat_sum$Raccoon[ dat_sum$Raccoon > 0 ] <- 1

# If we have no observations (lvls = NA, we make it state 1).
# state 1 is "no species present". This is where we use the 
# NA_missing object from above
dat_sum$State <- paste0( dat_sum$Coyote, dat_sum$Opossum, dat_sum$Raccoon )

# if no sampling done change State to NANANA
dat_sum$State[which(dat_sum$SurveyID %in% 
                        NA_missing$SurveyID[NA_missing$NSeasonGone==1])] <- "NANANA"

# convert the binary states to numeric states (i.e., categories),
# but make the NA values remain NA.
# combine lvls with dat via left join
dat_sum <- left_join( dat_sum, categories, by = "State" )
# make NA lvls 8 (the maximum possible state, which assumes all species present)
dat_sum$Lvls[is.na(dat_sum$Lvls)] <- NA

# These are the initial values to supply in the JAGS
# inits function in the init_functions.R script.
# it contains the maximum observed state assuming closure within
# a sampling period.
sp_inits <- matrix( dat_sum$Lvls, nrow = usi, ncol = uyr )

# this is the number of times each state was observed
# per season. nsite x nyear x total number of states. Each cell holds
# the number of times that state was observed. This is used
# to generate the binomial coefficient for model selection.
state_count <- array(0, dim = c(103,13,8))
for(i in 1:103){
  for(j in 1:13){
    a <- table(sp_ob_state[i,j,])
    if(length(a)==0) {
      next 
    } else {
      state_count[i,j, as.numeric(names(a))] <- as.numeric(a)
    }
  }
}

# The number of observations that occured at a site and season
n_obs <- apply( state_count, c( 1, 2 ), sum )
# The numerator of the binomial coefficient
numerator <- factorial( n_obs )
# an indicator variable that will convert an NA to 0
# Which will, in turn zero out the binomial coefficent
# when no samples were taken so that estimated states do not influence the 
# likelihood of the model.
ind_var <- n_obs 
ind_var[ind_var > 0] <- 1
# zero out numerator if NA
numerator[is.na(numerator)] <- 0
# give denominator the same dimensions as numerator
denominator <- n_obs

# the product of the factorial of state_count is the denominator
for(i in 1:103){
  for( j in 1:13){
    denominator[i, j] <- prod( factorial(state_count[i,j, ] ) )
  }
}

# zero it out if no samples were taken
binco <- (numerator / denominator) * ind_var

# Bring in the covariate data, the sites are ordered
# the same exact way as the species data above.
cdat <- read.csv("./data/fidino_covariate_data.csv")

# Do a pca with these data
cdat_pca <- prcomp(cdat[,-1], scale. = TRUE)

# make covariate matrices, multiplying by negative 1 
# so that positive values = higher housing denisty
x <- cbind(1, cdat_pca$x[,1] * -1)

# data_list to be supplied to JAGS. This includes:

#----------------------------#
#       SYMBOL NAMES         #
#----------------------------#
### psi   = initial occupancy
### gam   = colonization
### eps   = extinction
### rho   = detection
### pi    = inxs on colonization
### tau   = inxs on extinction
### delta = inxs on detection
#----------------------------#
#    OBJECT QUALIFIERS       # X = SYMBOL NAME (e.g., psi_cov)
#----------------------------#
### X_cov  = matrix of covaraties for a given process (including intercept)
### ncov_X = number of covariates for a given process (including intercept)
#----------------------------#
#      ADDITIONAL DATA       # 
#----------------------------#
### n_inxs = number of interactions in species pool
### nspec = number of species
### nyear = number of years sampled
### nsite = number of sites sampled
### nsurvey = number of days in secondary sampling period
### n_out = number of unique species categories
### binco = binomial coefficient calculated above
### state_count = number of times each state was observed by site and season
### rows_vec = sp_inxs indexing from make_inxs function
### cols_vec = sp_inxs indexing from make_inxs function
### y.sim = the observed data

# This is an intercept only matrix, used across a number of models
int_mat <- matrix(1, ncol = 1, nrow = 103)

# load glm module for JAGS
load.module("glm")
# This stores the cpo scores from each of the 6 models
cpo_scores <- data.frame(model = c("full", "inxs_URB_on_inxs",
                                   "inxs_URB_on_intercept",
                                   "inxs", "no_inxs_URB", "null"), cpo = 0)

##################################################
# Model 1: Full model, interactions with covariate
##################################################
data_list <- list(psi_cov = x, gam_cov = x, eps_cov = x, pi_cov = x,
                  tau_cov = x, delta_cov = int_mat,
                  rho_cov = x, 
                  ncov_psi = ncol(x), ncov_gam = ncol(x), 
                  ncov_eps = ncol(x), ncov_pi = ncol(x), 
                  ncov_tau = ncol(x), ncov_delta = 1, ncov_rho = 2,
                  n_inxs = (3^2)-3, nspec = 3, nyear = 13, nsite = nrow(x), 
                  nsurvey = 4, nout = 8, rows_vec = make_inxs(3)$rows_vec, 
                  cols_vec = make_inxs(3)$cols_vec, 
                  y = sp_ob_state, binco = binco, 
                  state_count = state_count,
                  pvec = rep(1/8, 8))

# track the z matrix so that we can calculate the likelihood
mout_full <- run.jags(model = "./jags_models/dcom_inxs.R",
  monitor = c("a", "b", "d", "f", "h", "g", "l", "b0", "d0", "z"), 
  data = data_list,
  n.chains = 7,
  inits = inits,
  adapt = 400,
  burnin = 100000,
  sample = ceiling(50000/7),
  thin = 10,
  summarise = FALSE,
  plots = FALSE,
  method = "parallel")

# save the output
saveRDS(mout_full, "./model_output/model_1_full.rds")
# Convert the jags output to a matrix
model_matrix <- as.matrix(as.mcmc.list(mout_full), chains = TRUE)
# calculate the cpo score
cpo_scores$cpo[1] <- calculate_cpo(model_matrix, data_list = data_list)
# remove these objects because they are large
rm(model_matrix)
rm(mout_full)
############################################
# model 2: inxs, but URB only on inxs
############################################
data_list <- list(psi_cov = int_mat, gam_cov = int_mat, eps_cov = int_mat,
                  pi_cov = x, tau_cov = x, delta_cov = int_mat,
                  rho_cov = x, 
                  ncov_psi = ncol(int_mat), ncov_gam = ncol(int_mat), 
                  ncov_eps = ncol(int_mat), ncov_pi = ncol(x), 
                  ncov_tau = ncol(x), ncov_delta = 1, ncov_rho = 2,
                  n_inxs = (3^2)-3, nspec = 3, nyear = 13, nsite = nrow(x), 
                  nsurvey = 4, nout = 8,binco = binco, 
                  state_count = state_count, rows_vec = make_inxs(3)$rows_vec, 
                  cols_vec = make_inxs(3)$cols_vec, 
                  y = sp_ob_state)

mout_inxs_inxcov <- run.jags(model = "./jags_models/dcom_inxs.R",
  monitor = c("a", "b", "d", "f", "h", "g", "l", "b0", "d0", "z"), 
  data = data_list,
  n.chains = 7,
  inits = inits,
  adapt = 400,
  burnin = 100000,
  sample = ceiling(50000/7),
  thin = 10,
  summarise = FALSE,
  plots = FALSE,
  method = "parallel")
# save output
saveRDS(mout_inxs_inxcov, "./model_output/model_2_inxs_inxcov.rds")
# convert to matrix
model_matrix <- as.matrix(as.mcmc.list(mout_inxs_inxcov), chains = TRUE)
# calculate cpo
cpo_scores$cpo[2] <- calculate_cpo(model_matrix, data_list = data_list)
# remove objects
rm(model_matrix)
rm(mout_inxs_inxcov)

#########################################
# Model 3: inxs, but no URB on inxs
#########################################
data_list <- list(psi_cov = x, gam_cov = x, eps_cov = x, pi_cov = int_mat,
                  tau_cov = int_mat, delta_cov = int_mat,
                  rho_cov = x, 
                  ncov_psi = ncol(x), ncov_gam = ncol(x), 
                  ncov_eps = ncol(x), ncov_pi = ncol(int_mat), 
                  ncov_tau = ncol(int_mat), ncov_delta = 1, ncov_rho = ncol(x),
                  n_inxs = (3^2)-3, nspec = 3, nyear = 13, nsite = nrow(x), 
                  nsurvey = 4, nout = 8,binco = binco, 
                  state_count = state_count, rows_vec = make_inxs(3)$rows_vec, 
                  cols_vec = make_inxs(3)$cols_vec, 
                  y = sp_ob_state)

mout_inxs_intcov <- run.jags(model = "./jags_models/dcom_inxs.R",
                  monitor = c("a", "b", "d", "f", "h", "g", "l", "b0", "d0", "z"), 
                  data = data_list,
                  n.chains = 7,
                  inits = inits,
                  adapt = 400,
                  burnin = 100000,
                  sample = ceiling(50000/7),
                  thin = 10,
                  summarise = FALSE,
                  plots = FALSE,
                  method = "parallel")
# save output
saveRDS(mout_inxs_intcov, "./model_output/model_3_inxs_intcov.rds")
# convert to matrix
model_matrix <- as.matrix(as.mcmc.list(mout_inxs_intcov), chains = TRUE)
# calculate cpo
cpo_scores$cpo[3] <- calculate_cpo(model_matrix, data_list = data_list)
# remove objects
rm(model_matrix)
rm(mout_inxs_intcov)

##############################################################
# Model 4: species interactions, but no urbanization covariate
##############################################################
data_list <- list(psi_cov = int_mat, gam_cov = int_mat, eps_cov = int_mat, 
                  pi_cov = int_mat, tau_cov = int_mat, 
                  delta_cov = int_mat,
                  rho_cov = int_mat,
                  ncov_psi = 1, ncov_gam = 1, 
                  ncov_eps = 1, ncov_pi = 1, 
                  ncov_tau = 1, ncov_delta = 1, ncov_rho = 1,
                  n_inxs = (3^2)-3, nspec = 3, nyear = 13, 
                  nsite = nrow(int_mat), 
                  nsurvey = 4, nout = 8,
                  rows_vec = make_inxs(3)$rows_vec, 
                  cols_vec = make_inxs(3)$cols_vec, 
                  y = sp_ob_state, binco = binco, 
                  state_count = state_count)

mout_int <- run.jags(model = "./jags_models/dcom_inxs.R",
                  monitor = c("a", "b", "d", "f", "h", "g", "l", "b0", "d0", "z"), 
                  data = data_list,
                  n.chains = 7,
                  inits = inits,
                  adapt = 400,
                  burnin = 100000,
                  sample = ceiling(50000/7),
                  thin = 10,
                  summarise = FALSE,
                  plots = FALSE,
                  method = "parallel")
# save the output
saveRDS(mout_int, "./model_output/model_4_intinxs.rds")
# convert to matrix
model_matrix <- as.matrix(as.mcmc.list(mout_int), chains = TRUE)
# calculate cpo
cpo_scores$cpo[4] <- calculate_cpo(model_matrix, data_list = data_list)
# remove objects
rm(model_matrix)
rm(mout_int)

#######################################################
# Model 5: no inxs, covariate on species specific rates
#######################################################
data_list <- list(psi_cov = x, gam_cov = x, eps_cov = x,
                  delta_cov = matrix(1, ncol = 1, nrow = 103),
                  ncov_psi = ncol(x), ncov_gam = ncol(x), 
                  ncov_eps = ncol(x), ncov_rho = ncol(x),
                  ncov_pi = 1, rho_cov = x, ncov_tau = 1, ncov_delta = 1,
                  pi_cov = int_mat, tau_cov = int_mat, 
                  n_inxs = (3^2)-3, nspec = 3, nyear = 13, nsite = nrow(x), 
                  nsurvey = 4, nout = 8,binco = binco, 
                  state_count = state_count, rows_vec = make_inxs(3)$rows_vec, 
                  cols_vec = make_inxs(3)$cols_vec, 
                  y = sp_ob_state)

mout_pca_noinxs <- run.jags(model = "./jags_models/dcom_no_inxs.R",
                monitor = c("a", "b", "d", "f", "h", "g", "l", "b0", "d0", "z"), 
                data = data_list,
                n.chains = 7,
                inits = inits_no_inxs,
                adapt = 400,
                burnin = 100000,
                sample = ceiling(50000/7),
                thin = 10,
                summarise = FALSE,
                plots = FALSE,
                method = "parallel")
# Save output
saveRDS(mout_pca_noinxs, "./model_output/model_5_pca_noinxs.rds")
# convert to matrix
model_matrix <- as.matrix(as.mcmc.list(mout_pca_noinxs), chains = TRUE)
# calculate cpo
cpo_scores$cpo[5] <- calculate_cpo(model_matrix, data_list = data_list)
# remove objects
rm(model_matrix)
rm(mout_pca_noinxs)


######################################################
# Model 6: Full intercept only model (no interactions)
######################################################
data_list <- list(psi_cov = int_mat, gam_cov = int_mat, eps_cov = int_mat, 
                  pi_cov = int_mat, tau_cov = int_mat, 
                  delta_cov = int_mat,
                  rho_cov = int_mat,
                  ncov_psi = 1, ncov_gam = 1, 
                  ncov_eps = 1, ncov_pi = 1, 
                  ncov_tau = 1, ncov_delta = 1, ncov_rho = 1,
                  n_inxs = (3^2)-3, nspec = 3, nyear = 13, 
                  nsite = nrow(int_mat), 
                  nsurvey = 4, nout = 8, 
                  rows_vec = make_inxs(3)$rows_vec, 
                  cols_vec = make_inxs(3)$cols_vec, 
                  y = sp_ob_state, binco = binco, 
                  state_count = state_count)

mout_int_noinxs <- run.jags(model = "./jags_models/dcom_no_inxs.R",
                monitor = c("a", "b", "d", "f", "h", "g", "l", "b0", "d0", "z"), 
                data = data_list,
                n.chains = 7,
                inits = inits_no_inxs,
                adapt = 400,
                burnin = 100000,
                sample = ceiling(50000/7),
                thin = 10,
                summarise = FALSE,
                plots = FALSE,
                method = "parallel")
# save model output
saveRDS(mout_int_noinxs, "./model_output/model_6_int_no_inxs.rds")
# convert to matrix
model_matrix <- as.matrix(as.mcmc.list(mout_int_noinxs), chains = TRUE)
# calculate cpo
cpo_scores$cpo[6] <- calculate_cpo(model_matrix, data_list = data_list)
# remove objects
rm(model_matrix)
rm(mout_int_noinxs)



# order the models best to worst and calculate delta cpo
final_cpo <- data.frame(model = cpo_scores$model, 
                        cpo = cpo_scores$cpo,
                        delta_cpo = cpo_scores$cpo - min(cpo_scores$cpo))
final_cpo <- final_cpo[order(final_cpo$cpo),]
write.csv(final_cpo, "fidino_cpo_scores.csv")