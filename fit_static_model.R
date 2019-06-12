library(AHMbook)
library(rjags)
library(runjags)

# bring in the data
data(HubbardBrook)
str(HubbardBrook)

# Get letter codes for the 13 species
speclist <- dimnames(HubbardBrook$counts)[[4]]

sel.species <- c(4, 7, 13) # Selecting BLBW, BTNW, and REVI

# Restrict all data to last 10 years: 2009-2018
year <- 2009:2018
str(counts <- HubbardBrook$counts[,,11:20,]) # Counts
str(dates <- HubbardBrook$dates[,,11:20])    # Survey dates
str(hours <- HubbardBrook$times[,,11:20])    # Survey hour (in hours)

# Date of full canopy leaf expansion as a phenological measure (see Lany et al. 2016)
budburst <- c(138.4, 133.8, 143.6, 133.9, 137.3, 146.6, 142.0, 148.3, 143.3, 141.0)
speclist <- speclist[sel.species] # Update species list

# get counts for these species and for 2016
counts <- counts[,,'2016',sel.species]
dates <- dates[,,'2016']
hours <- hours[,,'2016']

# standardize the detection covariates
DATES <- standardize(dates)
DATES[is.na(DATES)] <- 0
HOURS <- standardize(hours)
HOURS[is.na(HOURS)] <- 0

# Get sample sizes
nsite <- dim(counts)[1]
nrep <-  dim(counts)[2]
nspec <- dim(counts)[3]
ncat <- 2^nspec # number of possible community states


# collect and standardize the environmental covariates
elev <- HubbardBrook$sitecov$Elev
aspect <- HubbardBrook$sitecov$Aspect
elev <- (elev - 500) / 100 # center on 500 and scale by 100
north <- abs(HubbardBrook$sitecov$Aspect-180)/180 * pi

# model matrix for first order occupancy (psi)
#  there should include a column of 1's for the intercept
psi_cov <- matrix(NA, ncol = 5, nrow = nsite )
psi_cov[,1] <- 1 # intercept
psi_cov[,2] <- elev # elevation
psi_cov[,3] <- elev^2 # elevation squared
psi_cov[,4] <- north # north
psi_cov[,5] <- north^2 # north squared

# model matrix for second order occupancy (psi), letting all of the
#  covariates in psi_cov vary based off the first order ones
psi_inxs_cov <- psi_cov

# model matrix for first order detection (rho)
rho_cov <- array(NA, dim = c(nsite, nrep, 3))
rho_cov[,,1] <- 1 # intercept
rho_cov[,,2] <- DATES # date of count
rho_cov[,,3] <- HOURS # time duration of count

# model matrix for second order detection (rho)
#  Just assuming an intercept here as I do not think
#  the other covariates likely influence detection
#  via species inxs.
rho_inxs_cov <- rep(1, nsite)

# we need to change the species detections into categories
#  for each potential community state.

y <- counts
y[y>0] <- 1

# this condenses the count array to be site x survey
#  the values represent the species seen.
ycat <- apply(y, c(1,2), paste, collapse = "")

ycat[ycat =="000"] <- 1
ycat[ycat =="100"] <- 2
ycat[ycat =="010"] <- 3
ycat[ycat =="001"] <- 4
ycat[ycat =="110"] <- 5
ycat[ycat =="011"] <- 6
ycat[ycat =="101"] <- 7
ycat[ycat =="111"] <- 8
ycat[ycat =="NANANA"] <- NA # for surveys that did not happen

# convert each column to a numeric
ycat <- apply(ycat, 2, as.numeric)

# get the maximum possible state across all 3 potential
#  surveys at a site.

# returns a site x species matrix
#  values are number of surverys a species
#  was seen.
zinit <- apply(y, c(1,3), sum, na.rm = TRUE)

# make binary
zinit[zinit>1] <- 1

# convert to a category
zcat <- apply(zinit, 1, paste, collapse = '')

zcat[zcat =="000"] <- 1
zcat[zcat =="100"] <- 2
zcat[zcat =="010"] <- 3
zcat[zcat =="001"] <- 4
zcat[zcat =="110"] <- 5
zcat[zcat =="011"] <- 6
zcat[zcat =="101"] <- 7
zcat[zcat =="111"] <- 8

# make numeric again
zcat <- as.numeric(zcat)

# initial values
inits <- function() list(z = zcat)

# Parameters monitored
params <- c('betaA', 'betaB', 'betaC', 'betaAB', 'betaBC', 'betaAC',
            'alphaA', 'alphaB', 'alphaC', 'alphaAB', 'alphaBA',
            'alphaBC', 'alphaCB', 'alphaAC', 'alphaCA')

# MCMC settings
na <- 1000  ;  nc <- 3  ;  ni <- 30000  ;  nb <- 20000  ;  nt <- 1

bdata <- list(y = ycat, psi_cov = psi_cov,
              psi_inxs_cov = psi_inxs_cov, rho_cov = rho_cov,
              rho_inxs_cov = rho_inxs_cov, nsite = nsite,
              nrep = nrep, nfirst_order_psi = ncol(psi_cov),
              nsecond_order_psi = ncol(psi_inxs_cov),
              nfirst_order_rho = dim(rho_cov)[3],
              nsecond_order_rho = 1, ncat = ncat)

load.module('glm')


start_time <- Sys.time()
static_inxs <- run.jags('static_model.R',
                        monitor = params,
                        data = bdata,
                        n.chains = nc,
                        inits = inits,
                        burnin = nb,
                        sample = floor(ni/nc),
                        adapt = na,
                        thin = nt,
                        module = 'glm')
end_time <- Sys.time()

