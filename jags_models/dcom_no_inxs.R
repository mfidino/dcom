model{
######
# Detection model
######
for(j in 1:nsite) {
  for(ti in 1:nyear) {
    for(week in 1:nsurvey) {
      y[j, ti, week] ~ dcat( rdm[j, ( 1:nout ) , z[j, ti] ] )
    }
  }
}
#######
# Latent state model
#######
for(j in 1:nsite) {
  # for first season
  z[j, 1] ~ dcat( fsm[j, ( 1:nout )] )
  for(t in 2:nyear){
    z[j, t] ~ dcat( tpm[j, ( 1:nout ) , z[ j, t-1]] )
  }
}
for( j in 1:nsite ) {
######
# Fill in all of the transition probabilities
######
######
# Latent state for first season, fsm = transition vector for season 1
######
# First season probabilities for each state
fsm[j, 1] <- 1 #------------------------------------------------------------|U
fsm[j, 2] <- exp( psinit[1, j] ) #------------------------------------------|A
fsm[j, 3] <- exp( psinit[2, j] ) #------------------------------------------|B
fsm[j, 4] <- exp( psinit[3, j] ) #------------------------------------------|C
fsm[j, 5] <- exp( psinit[1, j] + psinit[2, j] ) #---------------------------|AB
fsm[j, 6] <- exp( psinit[2, j] + psinit[3, j] ) #---------------------------|BC
fsm[j, 7] <- exp( psinit[1, j] + psinit[3, j] ) #---------------------------|AC
fsm[j, 8] <- exp( psinit[1, j] + psinit[2, j] + psinit[3, j] ) #------------|ABC
######
# Latent state for dynamic part of model
# tpm = transition probability matrix. All columns sum to 1.
# dim(tpm)[1] = site j 
# dim(tpm)[2] = state at time t
# dim(tpm)[3] = state at time t-1
######
# U to ...
tpm[j, 1, 1] <- 1 #---------------------------------------------------------|U
tpm[j, 2, 1] <- exp( gam[1, j] ) #------------------------------------------|A
tpm[j, 3, 1] <- exp(             gam[2, j] ) #------------------------------|B
tpm[j, 4, 1] <- exp(                         gam[3, j] ) #------------------|C
tpm[j, 5, 1] <- exp( gam[1, j] + gam[2, j] ) #------------------------------|AB
tpm[j, 6, 1] <- exp( gam[1, j] +             gam[3, j] ) #------------------|AC
tpm[j, 7, 1] <- exp(             gam[2, j] + gam[3, j] ) #------------------|BC
tpm[j, 8, 1] <- exp( gam[1, j] + gam[2, j] + gam[3, j] ) # -----------------|ABC
# A to ...
tpm[j, 1, 2] <- exp( eps[1, j] ) #------------------------------------------|U
tpm[j, 2, 2] <- 1 #---------------------------------------------------------|A
tpm[j, 3, 2] <- exp( eps[1, j] + gam[2, j] ) #------------------------------|B
tpm[j, 4, 2] <- exp( eps[1, j] +             gam[3, j] ) #------------------|C
tpm[j, 5, 2] <- exp(             gam[2, j] ) #------------------------------|AB
tpm[j, 6, 2] <- exp(                         gam[3, j] ) #------------------|AC
tpm[j, 7, 2] <- exp( eps[1, j] + gam[2, j] + gam[3, j] ) #------------------|BC
tpm[j, 8, 2] <- exp(             gam[2, j] + gam[3, j] ) #------------------|ABC
# B to ...
tpm[j, 1, 3] <- exp(             eps[2, j] ) #------------------------------|U
tpm[j, 2, 3] <- exp( gam[1, j] + eps[2, j] ) #------------------------------|A
tpm[j, 3, 3] <- 1 #---------------------------------------------------------|B
tpm[j, 4, 3] <- exp(             eps[2, j] + gam[3, j] ) #------------------|C
tpm[j, 5, 3] <- exp( gam[1, j] ) #------------------------------------------|AB
tpm[j, 6, 3] <- exp( gam[1, j] + eps[2, j] + gam[3, j] ) #------------------|AC
tpm[j, 7, 3] <- exp(                         gam[3, j] ) #------------------|BC
tpm[j, 8, 3] <- exp( gam[1, j] +             gam[3, j] ) #------------------|ABC
# C to ...
tpm[j, 1, 4] <- exp(                         eps[3, j] ) #------------------|U
tpm[j, 2, 4] <- exp( gam[1, j] +             eps[3, j] ) #------------------|A
tpm[j, 3, 4] <- exp(             gam[2, j] + eps[3, j] ) #------------------|B
tpm[j, 4, 4] <- 1 #---------------------------------------------------------|C
tpm[j, 5, 4] <- exp( gam[1, j] + gam[2, j] + eps[3, j] ) #------------------|AB
tpm[j, 6, 4] <- exp( gam[1, j]) #-------------------------------------------|AC
tpm[j, 7, 4] <- exp(             gam[2, j] ) #------------------------------|BC
tpm[j, 8, 4] <- exp( gam[1, j] + gam[2, j] ) #------------------------------|ABC
# AB to ..
tpm[j, 1, 5] <- exp( eps[1, j] + eps[2, j] ) #------------------------------|U
tpm[j, 2, 5] <- exp(             eps[2, j] ) #------------------------------|A
tpm[j, 3, 5] <- exp( eps[1, j] ) #------------------------------------------|B
tpm[j, 4, 5] <- exp( eps[1, j] + eps[2, j] + gam[3, j] ) #------------------|C
tpm[j, 5, 5] <- 1 #---------------------------------------------------------|AB
tpm[j, 6, 5] <- exp(             eps[2, j] + gam[3, j] ) #------------------|AC
tpm[j, 7, 5] <- exp( eps[1, j] +             gam[3, j] ) #------------------|BC
tpm[j, 8, 5] <- exp(                         gam[3, j] ) #------------------|ABC
# AC to ...
tpm[j, 1, 6] <- exp( eps[1, j] +             eps[3, j] ) #------------------|U
tpm[j, 2, 6] <- exp(                         eps[3, j] ) #------------------|A
tpm[j, 3, 6] <- exp( eps[1, j] + gam[2, j] + eps[3, j] ) #------------------|B
tpm[j, 4, 6] <- exp( eps[1, j] ) #------------------------------------------|C
tpm[j, 5, 6] <- exp(             gam[2, j] + eps[3, j] ) #------------------|AB
tpm[j, 6, 6] <- 1 #---------------------------------------------------------|AC
tpm[j, 7, 6] <- exp( eps[1, j] + gam[2, j] ) #------------------------------|BC
tpm[j, 8, 6] <- exp(             gam[2, j] ) #------------------------------|ABC
# BC to ...
tpm[j, 1, 7] <- exp(             eps[2, j] + eps[3, j] ) #------------------|U
tpm[j, 2, 7] <- exp( gam[1, j] + eps[2, j] + eps[3, j] ) #------------------|A
tpm[j, 3, 7] <- exp(                         eps[3, j] ) #------------------|B
tpm[j, 4, 7] <- exp(             eps[2, j] ) #------------------------------|C
tpm[j, 5, 7] <- exp( gam[1, j] +             eps[3, j] ) #------------------|AB
tpm[j, 6, 7] <- exp( gam[1, j] + eps[2, j] ) #------------------------------|AC
tpm[j, 7, 7] <- 1 #---------------------------------------------------------|BC
tpm[j, 8, 7] <- exp( gam[1, j] ) #------------------------------------------|ABC
# ABC to ...
tpm[j, 1, 8] <- exp( eps[1, j] + eps[2, j] + eps[3, j] ) #------------------|U
tpm[j, 2, 8] <- exp(             eps[2, j] + eps[3, j] ) #------------------|A
tpm[j, 3, 8] <- exp( eps[1, j] +             eps[3, j] ) #------------------|B
tpm[j, 4, 8] <- exp( eps[1, j] + eps[2, j] ) #------------------------------|C
tpm[j, 5, 8] <- exp(                         eps[3, j] ) #------------------|AB
tpm[j, 6, 8] <- exp(             eps[2, j] ) #------------------------------|AC 
tpm[j, 7, 8] <- exp( eps[1, j] ) #------------------------------------------|BC
tpm[j, 8, 8] <- 1 #---------------------------------------------------------|ABC
######
# detection matrix (OS = observed state, TS = true state)
# rdm = rho detection matrix. Each row sums to 1.
# OS along rows, TS along columns
######
# TS = U
rdm[j, 1, 1] <- 1 #----------------------------------------------------|OS = U
rdm[j, 2, 1] <- 0 #----------------------------------------------------|OS = A
rdm[j, 3, 1] <- 0 #----------------------------------------------------|OS = B
rdm[j, 4, 1] <- 0 #----------------------------------------------------|OS = C
rdm[j, 5, 1] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 1] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 1] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 1] <- 0 #----------------------------------------------------|OS = ABC
# TS = A
rdm[j, 1, 2] <- 1 #----------------------------------------------------|OS = U
rdm[j, 2, 2] <- exp( rho[1, j] ) #-------------------------------------|OS = A
rdm[j, 3, 2] <- 0 #----------------------------------------------------|OS = B
rdm[j, 4, 2] <- 0 #----------------------------------------------------|OS = C
rdm[j, 5, 2] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 2] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 2] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 2] <- 0 #----------------------------------------------------|OS = ABC
# TS = B
rdm[j, 1, 3] <- 1 #----------------------------------------------------|OS = U
rdm[j, 2, 3] <- 0 #----------------------------------------------------|OS = A
rdm[j, 3, 3] <- exp(             rho[2, j] ) #-------------------------|OS = B
rdm[j, 4, 3] <- 0 #----------------------------------------------------|OS = C
rdm[j, 5, 3] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 3] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 3] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 3] <- 0 #----------------------------------------------------|OS = ABC
#TS = C
rdm[j, 1, 4] <- 1 #----------------------------------------------------|OS = U
rdm[j, 2, 4] <- 0 #----------------------------------------------------|OS = A
rdm[j, 3, 4] <- 0 #----------------------------------------------------|OS = B
rdm[j, 4, 4] <- exp(                         rho[3, j] ) #-------------|OS = C
rdm[j, 5, 4] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 4] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 4] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 4] <- 0 #----------------------------------------------------|OS = ABC
# TS = AB
rdm[j, 1, 5] <- 1 #----------------------------------------------------|OS = U
rdm[j, 2, 5] <- exp( rho[1, j] ) #-------------------------------------|OS = A
rdm[j, 3, 5] <- exp(             rho[2, j] ) #-------------------------|OS = B
rdm[j, 4, 5] <- 0 #----------------------------------------------------|OS = C
rdm[j, 5, 5] <- exp( rho[1, j] + rho[2, j] ) #-------------------------|OS = AB
rdm[j, 6, 5] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 5] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 5] <- 0 #----------------------------------------------------|OS = ABC
# TS = BC
rdm[j, 1, 6] <- 1 #----------------------------------------------------|OS = U
rdm[j, 2, 6] <- 0 #----------------------------------------------------|OS = A
rdm[j, 3, 6] <- exp(             rho[2, j] ) #-------------------------|OS = B
rdm[j, 4, 6] <- exp(                         rho[3, j] ) #-------------|OS = C
rdm[j, 5, 6] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 6] <- exp(             rho[2, j] + rho[3, j] ) #-------------|OS = BC
rdm[j, 7, 6] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 6] <- 0 #----------------------------------------------------|OS = ABC
# TS = AC
rdm[j, 1, 7] <- 1 #----------------------------------------------------|OS = U
rdm[j, 2, 7] <- exp( rho[1, j] ) #-------------------------------------|OS = A
rdm[j, 3, 7] <- 0 #----------------------------------------------------|OS = B
rdm[j, 4, 7] <- exp(                         rho[3, j] ) #-------------|OS = C
rdm[j, 5, 7] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 7] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 7] <- exp( rho[1, j] +             rho[3, j] ) #-------------|OS = AC
rdm[j, 8, 7] <- 0 #----------------------------------------------------|OS = ABC
# TS = ABC
rdm[j, 1, 8] <- 1 #|---------------------------------------------------|OS = U
rdm[j, 2, 8] <- exp( rho[1, j] ) #-------------------------------------|OS = A
rdm[j, 3, 8] <- exp(             rho[2, j] ) #-------------------------|OS = B
rdm[j, 4, 8] <- exp(                         rho[3, j] ) #-------------|OS = C
rdm[j, 5, 8] <- exp( rho[1, j] + rho[2, j] ) #-------------------------|OS = AB
rdm[j, 6, 8] <- exp(             rho[2, j] + rho[3, j] ) #-------------|OS = BC
rdm[j, 7, 8] <- exp( rho[1, j] +             rho[3, j] ) #-------------|OS = AC
rdm[j, 8, 8] <- exp( rho[1, j] + rho[2, j] + rho[3, j] ) #-------------|OS = ABC
######
# Fill in the linear predictors for the transition matrices
######
for( i in 1:nspec ) {
# base occupancy
psinit[i ,j]  <- inprod( a[i, ], psi_cov[j, ] )
# base colonization
gam[i, j] <-     inprod( b[i, ], gam_cov[j, ] ) 
# base extinction
eps[i, j] <-     inprod( d[i, ], eps_cov[j, ] ) 
# base detection probability
rho[i,j] <-      inprod( f[i, ], rho_cov[j, ] ) 
# 
}
} # closes for loop for j (sites) all the way up at the top of the model
#####
# Priors
######
  for(i in 1:nspec){
    # Initial Occupancy
    for( psip in 1:ncov_psi ){
      a[i, psip] ~ dlogis(0, 1)
    }
    # Colonization
    for( gamp in 1:ncov_gam ){
      b[i, gamp] ~ dlogis(0, 1)
    }  
    # Extinction
    for( epsp in 1:ncov_eps ){
      d[i, epsp] ~ dlogis(0, 1)
    }
    # Detection
    for( rhop in 1:ncov_rho ){
      f[i, rhop] ~ dlogis(0, 1)
    }
  }
}
