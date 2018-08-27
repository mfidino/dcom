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
  fin1[j] <- 1 - (psi_a[j] * (1 - psi_b[j]) * (1 - psi_c[j]) + 
    (1 - psi_a[j]) * psi_b[j] * (1 - psi_c[j]) +
    (1 - psi_a[j]) * (1 - psi_b[j]) * psi_c[j] + 
      psi_a[j] * psi_b[j] * (1 - psi_c[j]) +
      (1 - psi_a[j]) * psi_b[j] * psi_c[j] +
      psi_a[j] * (1 - psi_b[j]) * psi_c[j] +
      psi_a[j] * psi_b[j] * psi_c[j])
#
fsm[j, 1] <- fin1[j]  #------------|U
fsm[j, 2] <- psi_a[j] * (1 - psi_b[j]) * (1 - psi_c[j]) #-------------------|A
fsm[j, 3] <- (1 - psi_a[j]) * psi_b[j] * (1 - psi_c[j]) #-------------------|B
fsm[j, 4] <- (1 - psi_a[j]) * (1 - psi_b[j]) * psi_c[j] #-------------------|C
fsm[j, 5] <- psi_a[j] * psi_b[j] * (1 - psi_c[j]) #-------------------------|AB
fsm[j, 6] <- (1 - psi_a[j]) * psi_b[j] * psi_c[j] #-------------------------|BC
fsm[j, 7] <- psi_a[j] * (1 - psi_b[j]) * psi_c[j] #-------------------------|AC
fsm[j, 8] <- psi_a[j] * psi_b[j] * psi_c[j] #-------------------------------|ABC
######
# Latent state for dynamic part of model
# tpm = transition probability matrix. All columns sum to 1.
# dim(tpm)[1] = site j 
# dim(tpm)[2] = state at time t
# dim(tpm)[3] = state at time t-1
######
# U to ...
#fin2[j] <- 1 - 
#(gam_a[j] * (1 - gam_b[j]) * (1 - gam_c[j]) +
#(1 - gam_a[j]) * gam_b[j] * (1 - gam_c[j]) +
#(1 - gam_a[j]) * (1 - gam_b[j]) * gam_c[j] +
#gam_a[j] * gam_b[j] * (1 - gam_c[j]) +
#gam_a[j] * (1 - gam_b[j]) * gam_c[j] +
#(1 - gam_a[j]) * gam_b[j] * gam_c[j] +
#gam_a[j] * gam_b[j] * gam_c[j] )
#
tpm[j, 1, 1] <- 1 - tpm[j,2,1] - tpm[j,3,1] - tpm[j,4,1] - tpm[j,5,1] - 
                    tpm[j,6,1] - tpm[j,7,1] - tpm[j,8,1] #----------|U
tpm[j, 2, 1] <- gam_a[j] * (1 - gam_b[j]) * (1 - gam_c[j]) #----------------|A
tpm[j, 3, 1] <- (1 - gam_a[j]) * gam_b[j] * (1 - gam_c[j]) #----------------|B
tpm[j, 4, 1] <- (1 - gam_a[j]) * (1 - gam_b[j]) * gam_c[j] #----------------|C
tpm[j, 5, 1] <- gam_a[j] * gam_b[j] * (1 - gam_c[j]) #----------------------|AB
tpm[j, 6, 1] <- gam_a[j] * (1 - gam_b[j]) * gam_c[j] #----------------------|AC
tpm[j, 7, 1] <- (1 - gam_a[j]) * gam_b[j] * gam_c[j] #----------------------|BC
tpm[j, 8, 1] <- gam_a[j] * gam_b[j] * gam_c[j] # ---------------------------|ABC
# A to ...
#fin3[j] <- 1 - 
#(eps_a[j] * (1 - gam_b[j]) * (1 - gam_c[j])+
#eps_a[j] * gam_b[j] * (1 - gam_c[j])+ #---------------
#eps_a[j] * (1 - gam_b[j]) * gam_c[j] +#---------------
#(1 - eps_a[j]) * gam_b_a[j] * (1 - gam_c_a[j])+ #------
#(1 - eps_a[j]) * (1 - gam_b_a[j]) * gam_c_a[j]+ #------
#eps_a[j] * gam_b[j] * gam_c[j]+ #----------------------
#(1 - eps_a[j]) * gam_b_a[j] * gam_c_a[j]) #-------------
#
tpm[j, 1, 2] <- eps_a[j] * (1 - gam_b[j]) * (1 - gam_c[j]) #----------------|U
tpm[j, 2, 2] <- 1 - tpm[j,1,2] - tpm[j,3,2] - tpm[j,4,2] - tpm[j,5,2] - 
                    tpm[j,6,2] - tpm[j,7,2] - tpm[j,8,2] #------|A
tpm[j, 3, 2] <- eps_a[j] * gam_b[j] * (1 - gam_c[j]) #----------------------|B
tpm[j, 4, 2] <- eps_a[j] * (1 - gam_b[j]) * gam_c[j] #----------------------|C
tpm[j, 5, 2] <- (1 - eps_a[j]) * gam_b_a[j] * (1 - gam_c_a[j]) #------------|AB
tpm[j, 6, 2] <- (1 - eps_a[j]) * (1 - gam_b_a[j]) * gam_c_a[j] #------------|AC
tpm[j, 7, 2] <- eps_a[j] * gam_b[j] * gam_c[j] #----------------------------|BC
tpm[j, 8, 2] <- (1 - eps_a[j]) * gam_b_a[j] * gam_c_a[j] #------------------|ABC
# B to ...
#fin4[j] <- 1 - 
#((1 - gam_a[j]) * eps_b[j] * (1 - gam_c[j]) +
#gam_a[j] * eps_b[j] * (1 - gam_c[j]) +
#(1 - gam_a[j]) * eps_b[j] * gam_c[j] +
#gam_a_b[j] * (1 - eps_b[j]) * (1 - gam_c_b[j]) +
#gam_a[j] * eps_b[j] * gam_c[j]  +
#(1 - gam_a_b[j]) * (1 - eps_b[j]) * gam_c_b[j] +
#gam_a_b[j] * (1 - eps_b[j]) * gam_c_b[j] )
#
tpm[j, 1, 3] <- (1 - gam_a[j]) * eps_b[j] * (1 - gam_c[j]) #----------------|U
tpm[j, 2, 3] <- gam_a[j] * eps_b[j] * (1 - gam_c[j]) #----------------------|A
tpm[j, 3, 3] <- 1 - tpm[j,1,3] - tpm[j,2,3] - tpm[j,4,3] - tpm[j,5,3] - 
                    tpm[j,6,3] - tpm[j,7,3] - tpm[j,8,3]
tpm[j, 4, 3] <- (1 - gam_a[j]) * eps_b[j] * gam_c[j] #----------------------|C
tpm[j, 5, 3] <- gam_a_b[j] * (1 - eps_b[j]) * (1 - gam_c_b[j]) #------------|AB
tpm[j, 6, 3] <- gam_a[j] * eps_b[j] * gam_c[j]  #---------------------------|AC
tpm[j, 7, 3] <- (1 - gam_a_b[j]) * (1 - eps_b[j]) * gam_c_b[j] #------------|BC
tpm[j, 8, 3] <- gam_a_b[j] * (1 - eps_b[j]) * gam_c_b[j] #------------------|ABC
# C to ...
#fin5[j] <- 1 - 
#((1 - gam_a[j]) * (1 - gam_b[j]) * eps_c[j] +
#gam_a[j] * (1 - gam_b[j]) * eps_c[j] +
#(1 - gam_a[j]) * gam_b[j] * eps_c[j] +
#gam_a[j] * gam_b[j] * eps_c[j] +
#gam_a_c[j] * (1 - gam_b_c[j]) * (1 - eps_c[j]) +
#(1 - gam_a_c[j]) * gam_b_c[j] * (1 - eps_c[j]) +
#gam_a_c[j] * gam_b_c[j] * (1 - eps_c[j]) )
#
tpm[j, 1, 4] <- (1 - gam_a[j]) * (1 - gam_b[j]) * eps_c[j] #----------------|U
tpm[j, 2, 4] <- gam_a[j] * (1 - gam_b[j]) * eps_c[j] #----------------------|A
tpm[j, 3, 4] <- (1 - gam_a[j]) * gam_b[j] * eps_c[j] #----------------------|B
tpm[j, 4, 4] <- 1 - tpm[j,1,4] - tpm[j,2,4] - tpm[j,3,4] - tpm[j,5,4] - 
                    tpm[j,6,4] - tpm[j,7,4] - tpm[j,8,4] #------|C
tpm[j, 5, 4] <- gam_a[j] * gam_b[j] * eps_c[j] #----------------------------|AB
tpm[j, 6, 4] <- gam_a_c[j] * (1 - gam_b_c[j]) * (1 - eps_c[j])  #-----------|AC
tpm[j, 7, 4] <- (1 - gam_a_c[j]) * gam_b_c[j] * (1 - eps_c[j]) #------------|BC
tpm[j, 8, 4] <- gam_a_c[j] * gam_b_c[j] * (1 - eps_c[j]) #------------------|ABC
# AB to ..
#fin6[j] <- 1 -
#(eps_a[j] * eps_b[j] * (1 - gam_c[j]) +
#(1 - eps_a[j]) * eps_b_a[j] * (1 - gam_c_a[j]) +
#eps_a_b[j] * (1 - eps_b[j]) * (1 - gam_c_b[j]) +
#eps_a[j] * eps_b[j] * gam_c[j] +
#(1 - eps_a[j]) * eps_b_a[j] * gam_c_a[j] +
#eps_a_b[j] * (1 - eps_b[j]) * gam_c_b[j] +
#(1 - eps_a_b[j]) * (1 - eps_b_a[j]) * gam_c_ab[j] )
#
tpm[j, 1, 5] <- eps_a[j] * eps_b[j] * (1 - gam_c[j]) #----------------------|U
tpm[j, 2, 5] <- (1 - eps_a[j]) * eps_b_a[j] * (1 - gam_c_a[j]) #------------|A
tpm[j, 3, 5] <- eps_a_b[j] * (1 - eps_b[j]) * (1 - gam_c_b[j]) #------------|B
tpm[j, 4, 5] <- eps_a[j] * eps_b[j] * gam_c[j] #----------------------------|C
tpm[j, 5, 5] <- 1 - tpm[j,1,5] - tpm[j,2,5] - tpm[j,3,5] - tpm[j,4,5] - 
                    tpm[j,6,5] - tpm[j,7,5] - tpm[j,8,5] #---|AB
tpm[j, 6, 5] <- (1 - eps_a[j]) * eps_b_a[j] * gam_c_a[j] #------------------|AC
tpm[j, 7, 5] <- eps_a_b[j] * (1 - eps_b[j]) * gam_c_b[j] #------------------|BC
tpm[j, 8, 5] <- (1 - eps_a_b[j]) * (1 - eps_b_a[j]) * gam_c_ab[j] #---------|ABC
# AC to ...
#fin7[j] <- 1 - 
#(eps_a[j] * (1 - gam_b[j]) * eps_c[j] +
#(1 - eps_a[j]) * (1 - gam_b_a[j]) * eps_c_a[j] +
#eps_a[j] * gam_b[j] * eps_c[j] +
#eps_a_c[j] * (1 - gam_b_c[j]) * (1 - eps_c[j]) +
#(1 - eps_a[j]) * gam_b_a[j] * eps_c_a[j] +
#eps_a_c[j] * gam_b_c[j] * (1 - eps_c[j]) +
#(1 - eps_a_c[j]) * gam_b_ac[j] * (1 - eps_c_a[j]) )
#
tpm[j, 1, 6] <- eps_a[j] * (1 - gam_b[j]) * eps_c[j] #----------------------|U
tpm[j, 2, 6] <- (1 - eps_a[j]) * (1 - gam_b_a[j]) * eps_c_a[j] #------------|A
tpm[j, 3, 6] <- eps_a[j] * gam_b[j] * eps_c[j] #----------------------------|B
tpm[j, 4, 6] <- eps_a_c[j] * (1 - gam_b_c[j]) * (1 - eps_c[j]) #------------|C
tpm[j, 5, 6] <- (1 - eps_a[j]) * gam_b_a[j] * eps_c_a[j] #------------------|AB
tpm[j, 6, 6] <- 1 - tpm[j,1,6] - tpm[j,2,6] - tpm[j,3,6] - tpm[j,4,6] - 
                    tpm[j,5,6] - tpm[j,7,6] - tpm[j,8,6] #---|AC
tpm[j, 7, 6] <- eps_a_c[j] * gam_b_c[j] * (1 - eps_c[j]) #------------------|BC
tpm[j, 8, 6] <- (1 - eps_a_c[j]) * gam_b_ac[j] * (1 - eps_c_a[j]) #------------|ABC
# BC to ...
#fin8[j] <- 1 - 
#((1 - gam_a[j]) * eps_b[j] * eps_c[j] +
#gam_a[j] * eps_b[j] * eps_c[j] +
#(1 - gam_a_b[j]) * (1 - eps_b[j]) * eps_c_b[j] +
#(1 - gam_a_c[j]) * eps_b_c[j] * (1 - eps_c[j]) +
#gam_a_b[j] * (1 - eps_b[j]) * eps_c_b[j] +
#gam_a_c[j] * eps_b_c[j] * (1 - eps_c[j]) +
#gam_a_bc[j] * (1 - eps_b_c[j]) * (1 - eps_c_b[j]) )
#
tpm[j, 1, 7] <- (1 - gam_a[j]) * eps_b[j] * eps_c[j] #----------------------|U
tpm[j, 2, 7] <- gam_a[j] * eps_b[j] * eps_c[j] #----------------------------|A
tpm[j, 3, 7] <- (1 - gam_a_b[j]) * (1 - eps_b[j]) * eps_c_b[j] #------------|B
tpm[j, 4, 7] <- (1 - gam_a_c[j]) * eps_b_c[j] * (1 - eps_c[j]) #------------|C
tpm[j, 5, 7] <- gam_a_b[j] * (1 - eps_b[j]) * eps_c_b[j] #------------------|AB
tpm[j, 6, 7] <- gam_a_c[j] * eps_b_c[j] * (1 - eps_c[j]) #------------------|AC
tpm[j, 7, 7] <- 1 - tpm[j,1,7] - tpm[j,2,7] - tpm[j,3,7] - tpm[j,4,7] - 
                    tpm[j,5,7] - tpm[j,6,7] - tpm[j,8,7] #---|BC
tpm[j, 8, 7] <- gam_a_bc[j] * (1 - eps_b_c[j]) * (1 - eps_c_b[j]) #---------|ABC 
# ABC to ...
#fin9[j] <- 1 - 
#(eps_a[j] * eps_b[j] * eps_c[j] +
#(1 - eps_a[j]) * eps_b_a[j] * eps_c_a[j] +
#eps_a_b[j] * (1 - eps_b[j]) * eps_c_b[j] +
#eps_a_c[j] * eps_b_c[j] * (1 - eps_c[j]) +
#(1 - eps_a_b[j]) * (1 - eps_b_a[j]) * eps_c_ab[j] +
#(1 - eps_a_c[j]) * eps_b_ac[j] * (1 - eps_c_a[j]) +
#eps_a_bc[j] * (1 - eps_b_c[j]) * (1 - eps_c_b[j]) )
#
tpm[j, 1, 8] <- eps_a[j] * eps_b[j] * eps_c[j] #----------------------------|U
tpm[j, 2, 8] <- (1 - eps_a[j]) * eps_b_a[j] * eps_c_a[j] #------------------|A
tpm[j, 3, 8] <- eps_a_b[j] * (1 - eps_b[j]) * eps_c_b[j] #------------------|B
tpm[j, 4, 8] <- eps_a_c[j] * eps_b_c[j] * (1 - eps_c[j]) #------------------|C
tpm[j, 5, 8] <- (1 - eps_a_b[j]) * (1 - eps_b_a[j]) * eps_c_ab[j] #---------|AB
tpm[j, 6, 8] <- (1 - eps_a_c[j]) * eps_b_ac[j] * (1 - eps_c_a[j]) #---------|AC 
tpm[j, 7, 8] <- eps_a_bc[j] * (1 - eps_b_c[j]) * (1 - eps_c_b[j]) #---------|BC
tpm[j, 8, 8] <- 1 - tpm[j,1,8] - tpm[j,2,8] - tpm[j,3,8] - tpm[j,4,8] - 
                    tpm[j,5,8] - tpm[j,6,8] - tpm[j,7,8]
######
# detection matrix (OS = observed state, TS = true state)
# rdm = rho detection matrix. Each column sums to 1.
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
rdm[j, 1, 2] <- 1 - rho_a[j] #-----------------------------------------|OS = U
rdm[j, 2, 2] <- rho_a[j] #---------------------------------------------|OS = A
rdm[j, 3, 2] <- 0 #----------------------------------------------------|OS = B
rdm[j, 4, 2] <- 0 #----------------------------------------------------|OS = C
rdm[j, 5, 2] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 2] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 2] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 2] <- 0 #----------------------------------------------------|OS = ABC
# TS = B
rdm[j, 1, 3] <- 1 - rho_b[j] #-----------------------------------------|OS = U
rdm[j, 2, 3] <- 0 #----------------------------------------------------|OS = A
rdm[j, 3, 3] <- rho_b[j] #---------------------------------------------|OS = B
rdm[j, 4, 3] <- 0 #----------------------------------------------------|OS = C
rdm[j, 5, 3] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 3] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 3] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 3] <- 0 #----------------------------------------------------|OS = ABC
#TS = C
rdm[j, 1, 4] <- 1 - rho_c[j] #-----------------------------------------|OS = U
rdm[j, 2, 4] <- 0 #----------------------------------------------------|OS = A
rdm[j, 3, 4] <- 0 #----------------------------------------------------|OS = B
rdm[j, 4, 4] <- rho_c[j] #---------------------------------------------|OS = C
rdm[j, 5, 4] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 4] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 4] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 4] <- 0 #----------------------------------------------------|OS = ABC
# TS = AB
rdm[j, 1, 5] <- (1 - rho_a[j]) * (1 - rho_b_a[j]) #--------------------|OS = U
rdm[j, 2, 5] <- rho_a[j] * (1 - rho_b_a[j]) #--------------------------|OS = A
rdm[j, 3, 5] <- (1 - rho_a[j]) * rho_b_a[j] #--------------------------|OS = B
rdm[j, 4, 5] <- 0 #----------------------------------------------------|OS = C
rdm[j, 5, 5] <- rho_a[j] * rho_b_a[j] #--------------------------------|OS = AB
rdm[j, 6, 5] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 5] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 5] <- 0 #----------------------------------------------------|OS = ABC
# TS = BC
rdm[j, 1, 6] <- (1 - rho_b[j]) * (1 - rho_c[j]) #----------------------|OS = U
rdm[j, 2, 6] <- 0 #----------------------------------------------------|OS = A
rdm[j, 3, 6] <- rho_b[j] * (1 - rho_c[j]) #----------------------------|OS = B
rdm[j, 4, 6] <- (1 - rho_b[j]) * rho_c[j] #----------------------------|OS = C
rdm[j, 5, 6] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 6] <- rho_b[j] * rho_c[j] #----------------------------------|OS = BC
rdm[j, 7, 6] <- 0 #----------------------------------------------------|OS = AC
rdm[j, 8, 6] <- 0 #----------------------------------------------------|OS = ABC
# TS = AC
rdm[j, 1, 7] <- (1 - rho_a[j]) * (1 - rho_c_a[j]) #--------------------|OS = U
rdm[j, 2, 7] <- rho_a[j] * (1 - rho_c_a[j]) #--------------------------|OS = A
rdm[j, 3, 7] <- 0 #----------------------------------------------------|OS = B
rdm[j, 4, 7] <- (1 - rho_a[j]) * rho_c_a[j] #--------------------------|OS = C
rdm[j, 5, 7] <- 0 #----------------------------------------------------|OS = AB
rdm[j, 6, 7] <- 0 #----------------------------------------------------|OS = BC
rdm[j, 7, 7] <- rho_a[j] * rho_c_a[j] #--------------------------------|OS = AC
rdm[j, 8, 7] <- 0 #----------------------------------------------------|OS = ABC
# TS = ABC
rdm[j, 1, 8] <- (1 - rho_a[j]) * (1 - rho_b_a[j]) * (1 - rho_c_a[j]) #-|OS = U
rdm[j, 2, 8] <- rho_a[j] * (1 - rho_b_a[j]) * (1 - rho_c_a[j]) #-------|OS = A
rdm[j, 3, 8] <- (1 - rho_a[j]) * rho_b_a[j] * (1 - rho_c_a[j]) #-------|OS = B
rdm[j, 4, 8] <- (1 - rho_a[j]) * (1 - rho_b_a[j]) * rho_c_a[j] #-------|OS = C
rdm[j, 5, 8] <- rho_a[j] * rho_b_a[j] * (1 - rho_c_a[j]) #-------------|OS = AB
rdm[j, 6, 8] <- (1 - rho_a[j]) * rho_b_a[j] * rho_c_a[j] #-------------|OS = BC
rdm[j, 7, 8] <- rho_a[j] * (1 - rho_b_a[j]) * rho_c_a[j] #-------------|OS = AC
rdm[j, 8, 8] <- rho_a[j] * rho_b_a[j] * rho_c_a[j] #-------------------|OS = ABC
######
# Fill in the linear predictors for the transition matrices
######
#
# base occupancy
psi_a[j] <- ilogit( inprod( a[1, ], psi_cov[j, ] ) )
psi_b[j] <- ilogit( inprod( a[2, ], psi_cov[j, ] ) )
psi_c[j] <- ilogit( inprod( a[3, ], psi_cov[j, ] ) )
# base colonization
gam_a[j] <- ilogit( inprod( b[1, ], gam_cov[j, ] ) )
gam_b[j] <- ilogit( inprod( b[2, ], gam_cov[j, ] ) )
gam_c[j] <- ilogit( inprod( b[3, ], gam_cov[j, ] ) )
# colonization species a inxs
gam_a_b[j] <- ilogit( inprod(b[1, ], gam_cov[j, ]) +  
                      inprod(g[1, ], pi_cov[j, ]) ) # a|b
gam_a_c[j] <- ilogit( inprod(b[1, ], gam_cov[j, ]) + 
                      inprod(g[2, ], pi_cov[j, ]) ) # a|c
gam_a_bc[j] <- ilogit( inprod(b[1, ], gam_cov[j, ]) + 
                       inprod(g[1, ], pi_cov[j, ]) + # a|b
                       inprod(g[2, ], pi_cov[j, ]) ) # a|c
# colonization species b inxs
gam_b_a[j] <- ilogit( inprod(b[2, ], gam_cov[j, ]) + 
                      inprod(g[3, ], pi_cov[j, ]) ) # b|a
gam_b_c[j] <- ilogit( inprod(b[2, ], gam_cov[j, ]) + 
                      inprod(g[4, ], pi_cov[j, ]) ) # b|c
gam_b_ac[j] <- ilogit( inprod(b[2, ], gam_cov[j, ]) + 
                       inprod(g[3, ], pi_cov[j, ]) + # b|a
                       inprod(g[4, ], pi_cov[j, ]) ) # b|c
# colonization species c inxs
gam_c_a[j] <- ilogit( inprod(b[3, ], gam_cov[j, ]) + 
                      inprod(g[5, ], pi_cov[j, ]) ) # c|a
gam_c_b[j] <- ilogit( inprod(b[3, ], gam_cov[j, ]) + 
                      inprod(g[6, ], pi_cov[j, ]) ) # c|b
gam_c_ab[j] <- ilogit( inprod(b[3, ], gam_cov[j, ]) + 
                       inprod(g[5, ], pi_cov[j, ]) + # c|a
                       inprod(g[6, ], pi_cov[j, ]) ) # c|b
# base extinction
eps_a[j] <- ilogit( inprod( d[1, ], eps_cov[j, ] ) )
eps_b[j] <- ilogit( inprod( d[2, ], eps_cov[j, ] ) )
eps_c[j] <- ilogit( inprod( d[3, ], eps_cov[j, ] ) )
# extinction species a inxs
eps_a_b[j] <- ilogit( inprod(d[1, ], eps_cov[j, ]) +  
                      inprod(h[1, ], tau_cov[j, ]) ) # a|b
eps_a_c[j] <- ilogit( inprod(d[1, ], eps_cov[j, ]) + 
                      inprod(h[2, ], tau_cov[j, ]) ) # a|c
eps_a_bc[j] <- ilogit( inprod(d[1, ], eps_cov[j, ]) + 
                         inprod(h[1, ], tau_cov[j, ]) + # a|b
                         inprod(h[2, ], tau_cov[j, ]) ) # a|c
# extinction species b inxs
eps_b_a[j] <- ilogit( inprod(d[2, ], eps_cov[j, ]) + 
                        inprod(h[3, ], tau_cov[j, ]) ) # b|a
eps_b_c[j] <- ilogit( inprod(d[2, ], eps_cov[j, ]) + 
                        inprod(h[4, ], tau_cov[j, ]) ) # b|c
eps_b_ac[j] <- ilogit( inprod(d[2, ], eps_cov[j, ]) + 
                         inprod(h[3, ], tau_cov[j, ]) + # b|a
                         inprod(h[4, ], tau_cov[j, ]) ) # b|c
# extinction species c inxs
eps_c_a[j] <- ilogit( inprod(d[3, ], eps_cov[j, ]) + 
                        inprod(h[5, ], tau_cov[j, ]) ) # c|a
eps_c_b[j] <- ilogit( inprod(d[3, ], eps_cov[j, ]) + 
                        inprod(h[6, ], tau_cov[j, ]) ) # c|b
eps_c_ab[j] <- ilogit( inprod(d[3, ], eps_cov[j, ]) + 
                         inprod(h[5, ], tau_cov[j, ]) + # c|a
                         inprod(h[6, ], tau_cov[j, ]) ) # c|b
# base detection probability
rho_a[j] <- ilogit( inprod(f[1, ], rho_cov[j,]) )
rho_b[j] <- ilogit( inprod(f[2, ], rho_cov[j,]) )
rho_c[j] <- ilogit( inprod(f[3, ], rho_cov[j,]) )
# species b inxs
rho_b_a[j] <- ilogit( inprod(f[2, ], rho_cov[j,]) + 
                      inprod(l[1, ], delta_cov[j,]) )
# species c inxs
rho_c_a[j] <- ilogit( inprod(f[3, ], rho_cov[j,]) + 
                        inprod(l[2, ], delta_cov[j,]) )
} # closes for loop for j (sites) all the way up at the top of the model
#####
# Priors
######
for(i in 1:nspec){
# Initial Occupancy
  for( psip in 1:ncov_psi ){
    a[i, psip] ~ dlogis( 0, 1 )
  }
# Colonization
  for( gamp in 1:ncov_gam ){
    b[i, gamp] ~ dlogis( 0, 1 )
  }  
# Extinction
  for( epsp in 1:ncov_eps ){
    d[i, epsp] ~ dlogis( 0, 1 )
  }
# Detection
  for( rhop in 1:ncov_rho ){
    f[i, rhop] ~ dlogis( 0, 1 )
  }
}
for( k in 1:n_inxs ){
# Inxs on colonization
  for( pip in 1:ncov_pi ){
    g[k, pip] ~ dlogis( 0, 1 )
  }
# Inxs on extinction
  for( taup in 1:ncov_tau ){
    h[k, taup] ~ dlogis( 0, 1 )
  }
}
# Inxs on detection
for( delp in 1:ncov_delta ){
  l[1, delp] ~ dlogis( 0, 1 )
  l[2, delp] ~ dlogis( 0, 1 )
}
}