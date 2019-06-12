library(AHMbook)
library(rjags)
library(runjags)

# bring in the data
data(HubbardBrook)

# Get letter codes for the 13 species
speclist <- dimnames(HubbardBrook$counts)[[4]]

sel.species <- c(4, 7, 13) # Selecting BLBW, BTNW, and REVI

# Restrict all data to last 10 years: 2009-2018
year <- 2009:2018
str(counts <- HubbardBrook$counts[,,11:20,]) # Counts
str(dates <- HubbardBrook$dates[,,11:20])    # Survey dates
str(hours <- HubbardBrook$times[,,11:20])    # Survey hour (in hours)

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



cat(file = 'static_categorical.txt', "
model{
  ######
  # Detection model
  ######
  for(i in 1:nsite) {
      for(j in 1:nrep) {
        y[i, j] ~ dcat( rdm[i, j, ( 1:ncat ) , z[i] ] )
      }
  }
  #######
  # Latent state model
  #######
  for(i in 1:nsite) {
    z[i] ~ dcat( lsp[i, ( 1:ncat )] )
  }
  for( i in 1:nsite ) {
    ######
    # Latent state probabilities (lsp)
    ######
    #  Probabilities for each state
    lsp[i, 1] <- 1 #----------------------------------------|U
    lsp[i, 2] <- exp( psiA[i] ) #---------------------------|A
    lsp[i, 3] <- exp( psiB[i] ) #---------------------------|B
    lsp[i, 4] <- exp( psiC[i] ) #---------------------------|C
    lsp[i, 5] <- exp( psiAB[i] ) #--------------------------|AB
    lsp[i, 6] <- exp( psiBC[i] ) #--------------------------|BC
    lsp[i, 7] <- exp( psiAC[i] ) #--------------------------|AC
    lsp[i, 8] <- exp( psiABC[i]  ) #------------------------|ABC
    for(j in 1:nrep){
    ######
    # detection matrix (OS = observed state, TS = true state)
    # rdm = rho detection matrix. Each row sums to 1.
    # OS along rows, TS along columns
    ######
    # TS = U
    rdm[i, j, 1, 1] <- 1 #----------------------------------|OS = U
    rdm[i, j, 2, 1] <- 0 #----------------------------------|OS = A
    rdm[i, j, 3, 1] <- 0 #----------------------------------|OS = B
    rdm[i, j, 4, 1] <- 0 #----------------------------------|OS = C
    rdm[i, j, 5, 1] <- 0 #----------------------------------|OS = AB
    rdm[i, j, 6, 1] <- 0 #----------------------------------|OS = BC
    rdm[i, j, 7, 1] <- 0 #----------------------------------|OS = AC
    rdm[i, j, 8, 1] <- 0 #----------------------------------|OS = ABC
    # TS = A
    rdm[i, j, 1, 2] <- 1 #----------------------------------|OS = U
    rdm[i, j, 2, 2] <- exp( rhoA[i, j] ) #------------------|OS = A
    rdm[i, j, 3, 2] <- 0 #----------------------------------|OS = B
    rdm[i, j, 4, 2] <- 0 #----------------------------------|OS = C
    rdm[i, j, 5, 2] <- 0 #----------------------------------|OS = AB
    rdm[i, j, 6, 2] <- 0 #----------------------------------|OS = BC
    rdm[i, j, 7, 2] <- 0 #----------------------------------|OS = AC
    rdm[i, j, 8, 2] <- 0 #----------------------------------|OS = ABC
    # TS = B
    rdm[i, j, 1, 3] <- 1 #----------------------------------|OS = U
    rdm[i, j, 2, 3] <- 0 #----------------------------------|OS = A
    rdm[i, j, 3, 3] <- exp( rhoB[i, j] ) #------------------|OS = B
    rdm[i, j, 4, 3] <- 0 #----------------------------------|OS = C
    rdm[i, j, 5, 3] <- 0 #----------------------------------|OS = AB
    rdm[i, j, 6, 3] <- 0 #----------------------------------|OS = BC
    rdm[i, j, 7, 3] <- 0 #----------------------------------|OS = AC
    rdm[i, j, 8, 3] <- 0 #----------------------------------|OS = ABC
    #TS = C
    rdm[i, j, 1, 4] <- 1 #----------------------------------|OS = U
    rdm[i, j, 2, 4] <- 0 #----------------------------------|OS = A
    rdm[i, j, 3, 4] <- 0 #----------------------------------|OS = B
    rdm[i, j, 4, 4] <- exp( rhoC[i, j] ) #------------------|OS = C
    rdm[i, j, 5, 4] <- 0 #----------------------------------|OS = AB
    rdm[i, j, 6, 4] <- 0 #----------------------------------|OS = BC
    rdm[i, j, 7, 4] <- 0 #----------------------------------|OS = AC
    rdm[i, j, 8, 4] <- 0 #----------------------------------|OS = ABC
    # TS = AB
    rdm[i, j, 1, 5] <- 1 #----------------------------------|OS = U
    rdm[i, j, 2, 5] <- exp( rhoAB[i, j] ) #-----------------|OS = A
    rdm[i, j, 3, 5] <- exp( rhoBA[i, j] ) #-----------------|OS = B
    rdm[i, j, 4, 5] <- 0 #----------------------------------|OS = C
    rdm[i, j, 5, 5] <- exp( rhoAB[i, j] + rhoBA[i, j]) #----|OS = AB
    rdm[i, j, 6, 5] <- 0 #----------------------------------|OS = BC
    rdm[i, j, 7, 5] <- 0 #----------------------------------|OS = AC
    rdm[i, j, 8, 5] <- 0 #----------------------------------|OS = ABC
    # TS = BC
    rdm[i, j, 1, 6] <- 1 #----------------------------------|OS = U
    rdm[i, j, 2, 6] <- 0 #----------------------------------|OS = A
    rdm[i, j, 3, 6] <- exp( rhoBC[i, j] ) #-----------------|OS = B
    rdm[i, j, 4, 6] <- exp( rhoCB[i, j] ) #-----------------|OS = C
    rdm[i, j, 5, 6] <- 0 #----------------------------------|OS = AB
    rdm[i, j, 6, 6] <- exp( rhoBC[i, j] + rhoCB[i, j] ) #---|OS = BC
    rdm[i, j, 7, 6] <- 0 #----------------------------------|OS = AC
    rdm[i, j, 8, 6] <- 0 #----------------------------------|OS = ABC
    # TS = AC
    rdm[i, j, 1, 7] <- 1 #----------------------------------|OS = U
    rdm[i, j, 2, 7] <- exp( rhoAC[i, j] ) #-----------------|OS = A
    rdm[i, j, 3, 7] <- 0 #----------------------------------|OS = B
    rdm[i, j, 4, 7] <- exp( rhoCA[i, j] ) #-----------------|OS = C
    rdm[i, j, 5, 7] <- 0 #----------------------------------|OS = AB
    rdm[i, j, 6, 7] <- 0 #----------------------------------|OS = BC
    rdm[i, j, 7, 7] <- exp( rhoAC[i, j] + rhoCA[i, j] ) #---|OS = AC
    rdm[i, j, 8, 7] <- 0 #----------------------------------|OS = ABC
    # TS = ABC
    rdm[i, j, 1, 8] <- 1 #----------------------------------|OS = U
    rdm[i, j, 2, 8] <- exp( rhoABC[i, j] ) #----------------|OS = A
    rdm[i, j, 3, 8] <- exp( rhoBAC[i, j] ) #----------------|OS = B
    rdm[i, j, 4, 8] <- exp( rhoCAB[i, j] ) #----------------|OS = C
    rdm[i, j, 5, 8] <- exp( rhoABC[i, j] + rhoBAC[i, j] ) #-|OS = AB
    rdm[i, j, 6, 8] <- exp( rhoBAC[i, j] + rhoCAB[i, j] ) #-|OS = BC
    rdm[i, j, 7, 8] <- exp( rhoABC[i, j] + rhoCAB[i, j] ) #-|OS = AC
    rdm[i, j, 8, 8] <- exp( rhoABC[i, j] + rhoBAC[i, j]  +
                            rhoCAB[i, j] ) #----------------|OS = ABC
    }
    ######
    # Occupancy linear predictors...
    ######
    # ...for states A, B, and C
    psiA[i]  <- inprod( betaA, psi_cov[i, ] )
    psiB[i]  <- inprod( betaB, psi_cov[i, ] )
    psiC[i]  <- inprod( betaC, psi_cov[i, ] )
    # ...for states AB, BC, and AC (in that order)
    psiAB[i] <-  psiA[i] + psiB[i] + inprod( betaAB,
                                             psi_inxs_cov[i, ] )
    psiBC[i] <-  psiB[i] + psiC[i] + inprod( betaBC,
                                             psi_inxs_cov[i, ] )
    psiAC[i] <-  psiA[i] + psiC[i] + inprod( betaAC,
                                             psi_inxs_cov[i, ] )
    # ...for state ABC
    psiABC[i] <-  psiA[i] + psiB[i] + psiC[i] +
      inprod( betaAB, psi_inxs_cov[i, ] ) +
      inprod( betaBC, psi_inxs_cov[i, ] ) +
      inprod( betaAC, psi_inxs_cov[i, ] )
    #######
    # Detection linear predictors
    #######
    for(j in 1:nrep){
    # These are the baseline detection linear predictors
    #  that do not incorporate interactions.
    rhoA[i, j] <- inprod( alphaA, rho_cov[i, j, ] )
    rhoB[i, j] <- inprod( alphaB, rho_cov[i, j, ] )
    rhoC[i, j] <- inprod( alphaC, rho_cov[i, j, ] )
    # These are the asymmetric interactions between all 3 species
    rhoAB[i, j] <- rhoA[i, j] + inprod( alphaAB, rho_inxs_cov[i] )
    rhoAC[i, j] <- rhoA[i, j] + inprod( alphaAC, rho_inxs_cov[i] )
    rhoBA[i, j] <- rhoB[i, j] + inprod( alphaBA, rho_inxs_cov[i] )
    rhoBC[i, j] <- rhoB[i, j] + inprod( alphaBC, rho_inxs_cov[i] )
    rhoCA[i, j] <- rhoC[i, j] + inprod( alphaCA, rho_inxs_cov[i] )
    rhoCB[i, j] <- rhoC[i, j] + inprod( alphaCB, rho_inxs_cov[i] )
    # These are the aymmetric interactions when all 3 species
    #  are present
    rhoABC[i, j] <- rhoA[i, j] + inprod(alphaAB, rho_inxs_cov[i]) +
      inprod(alphaAC, rho_inxs_cov[i])
    rhoBAC[i, j] <- rhoB[i, j] + inprod(alphaBA, rho_inxs_cov[i]) +
      inprod(alphaBC, rho_inxs_cov[i])
    rhoCAB[i, j] <- rhoC[i, j] + inprod(alphaCA, rho_inxs_cov[i]) +
      inprod(alphaCB, rho_inxs_cov[i])
    }
  } # closes the j loop for sites ALL the way at the top of the model
  #####
  # Priors
  ######
  # For first order psi priors
  for(fo_psi in 1:nfirst_order_psi){
  betaA[fo_psi] ~ dnorm(0, 0.1)
  betaB[fo_psi] ~ dnorm(0, 0.1)
  betaC[fo_psi] ~ dnorm(0, 0.1)
  }
  # for second order psi priors
  for(so_psi in 1:nsecond_order_psi){
    betaAB[so_psi] ~ dnorm(0, 0.1)
    betaAC[so_psi] ~ dnorm(0, 0.1)
    betaBC[so_psi] ~ dnorm(0, 0.1)
  }
  # for first order rho priors
  for(fo_rho in 1:nfirst_order_rho){
    alphaA[fo_rho] ~ dnorm(0, 0.1)
    alphaB[fo_rho] ~ dnorm(0, 0.1)
    alphaC[fo_rho] ~ dnorm(0, 0.1)
  }
  # for second order rho priors
  for(so_rho in 1:nsecond_order_rho){
    alphaAB[so_rho] ~ dnorm(0, 0.1)
    alphaBA[so_rho] ~ dnorm(0, 0.1)
    alphaAC[so_rho] ~ dnorm(0, 0.1)
    alphaCA[so_rho] ~ dnorm(0, 0.1)
    alphaBC[so_rho] ~ dnorm(0, 0.1)
    alphaCB[so_rho] ~ dnorm(0, 0.1)
  }
}
  "
  )

load.module('glm')

static_inxs <- run.jags('static_categorical.txt',
                        monitor = params,
                        data = bdata,
                        n.chains = nc,
                        inits = inits,
                        burnin = nb,
                        sample = floor(ni/nc),
                        adapt = na,
                        thin = nt,
                        module = 'glm')
