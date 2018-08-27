###############################
#
# calculate_cpo
#
# Written by Mason Fidino
#
###############################

# This function calculates the cpo for each datapoint and creates
#  the summary statsitic -log(sum(cpo)) for each model. For clarity,
#  I've written this to loop through each step in the MCMC chain.
#  This slows down the code considerably. This function could be made faster
#  by applying matrix algebra to the whole mcmc chain instead,
#  but requires some manipulations to the transition matrices to do so.

# This was also set up for our own analysis where we only have coyote
#  influence the detection rates of the smaller species (and not the other way
#  around). To make it for all of the species the rho detection matrix would
#  need to be updated as well as the number of rows in the matrix l.

# Arguments for this function
# --------------------------- #
# mm = model matrix of jags output as a matrix
#  e.g., mm <- as.matrix(jags_output, chains = TRUE) via coda package
#  columns of this matrix represent different parameters while rows are
#  posterior simulations

# data_list = the same data_list supplied to JAGS to fit the DCOM
#   This was generated in fit_softmax_model.R

calculate_cpo <- function(mm = NULL, data_list = NULL){

  # ensure mm is a matrix
  if(!is.matrix(mm)) {
    stop("The object mm must be a matrix")
  }
  
  # Assigning all the objects in data_list so I don't have to
  # index them the whole time.
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  # This is the likelihood matrix, which stores the likelihood of each
  # observations based off the parameters in the model for each step
  # of the mcmc chain.
  lik_mat <- matrix(0, ncol = nyear * nsite, nrow = nrow(mm))
  for(o in 1:nrow(mm)){ # going through each step of the mcmc chain
  
    # initial occupancy
    a <- matrix(mm[o,grep("a", colnames(mm))], ncol = ncov_psi, nrow = nspec)
    # colonization 
    b <- matrix(mm[o,grep("b\\[", colnames(mm))], ncol = ncov_gam, nrow = nspec)
    # extinction
    d <- matrix(mm[o,grep("d\\[", colnames(mm))], ncol = ncov_eps, nrow = nspec)
    # detection
    f <- matrix(mm[o,grep("f", colnames(mm))], ncol = ncov_rho, nrow = nspec)
    # colonization inxs
    g <- matrix(mm[o,grep("g", colnames(mm))], ncol = ncov_pi, nrow = n_inxs) 
    if(sum(is.na(g)>0)) g[is.na(g)] <- 0 # make 0 if not in the model
    # extinciton inxs
    h <- matrix(mm[o,grep("h", colnames(mm))], ncol = ncov_tau, nrow = n_inxs)
    if(sum(is.na(h)>0)) h[is.na(h)] <- 0 # make 0 if not in the model
    # detection inxs (remove -1 if doing all species)
    l <- matrix(mm[o,grep("l", colnames(mm))], ncol = ncov_delta, nrow = nspec - 1)
    if(sum(is.na(l)>0)) l[is.na(l)] <- 0 # make 0 if not in the model
    # estimated community state at time t and site j
    z <- matrix(mm[o,grep("z", colnames(mm))], ncol = nyear, nrow = nsite)

    # occupancy linear predictor
    psinit  <- a %*% t(psi_cov)
    # colonization linear predictor
    gam <-     b %*% t(gam_cov)
    # extinction linear predictor
    eps <-     d %*% t(eps_cov)
    # detection linear predictor
    rho <-     f %*% t(rho_cov)
    # set up arrays for the species interactions.
    # these arrays get filled in this order:
    #######
    #0|4|5#
    #-|-|-#
    #1|0|6# Zeros along the diagonal
    #-|-|-#
    #2|3|0#
    #######
    # We have set it up in the model such that the species in the off-diagonal
    # influence the species along the diagonal. For example, pi[3,2] is the 
    # influence species 2 has on the colonization rate of species 3 while 
    # pi[2, 1] is the influence species 1 has on species 2. More generally it
    # is pi[species of interest, species conditioning on]. 
    pi <- tau <-  gam_one <- eps_one <- array(0, dim = c(nspec, nspec, nsite))
    
    for(k in 1:n_inxs){
      # inxs on colonization 
      pi[rows_vec[k], cols_vec[k], ] <-  pi_cov[j, ] %*% g[k, ]
      # inxs on extinction 
      tau[rows_vec[k], cols_vec[k], ] <- tau_cov[j, ] %*% h[k, ] 
      # linear predictor for colonization | one species present at t-1
      gam_one[rows_vec[k], cols_vec[k], ] <- 
         b[rows_vec[k], ] %*%  t(gam_cov[, ]) + pi[rows_vec[k], cols_vec[k], ]
      # linear predictor for extinction | one species present at t-1
      eps_one[rows_vec[k], cols_vec[k], ] <- 
         d[rows_vec[k], ] %*%t(eps_cov[, ]) + tau[rows_vec[k], cols_vec[k], ]
    }
    gam_two <- eps_two <- matrix(0, nrow = nspec, ncol = nsite)
    
    
    # make linear predictors for eps and gam when 2 species are present
    gam_two[1, ] <- b[1, ] %*%  t(gam_cov[, ]) + pi[1,2, ] + pi[1,3, ]
    gam_two[2, ] <- b[2, ] %*%  t(gam_cov[, ]) + pi[2,1, ] + pi[2,3, ]
    gam_two[3, ] <- b[3, ] %*%  t(gam_cov[, ]) + pi[3,1, ] + pi[3,2, ]
    # make linear predictors for eps and gam when 2 species are present
    eps_two[1, ] <- d[1, ] %*%  t(eps_cov[, ]) + tau[1,2, ] + tau[1,3, ]
    eps_two[2, ] <- d[2, ] %*%  t(eps_cov[, ]) + tau[2,1, ] + tau[2,3, ]
    eps_two[3, ] <- d[3, ] %*%  t(eps_cov[, ]) + tau[3,1, ] + tau[3,2, ]
    
    # coyote influence on detection for opossum and raccoon
    delta  <-  f[2:3,] %*% t(rho_cov) +  l %*% t(delta_cov)
    
    # This is the numerator of the softmax function for each of the 8 community
    # states a site can be in.
    fsm <- matrix(0, ncol = 8, nrow = nsite)
    fsm[, 1] <- 1 #---------------------------------------------------------|U
    fsm[, 2] <- exp( psinit[1, ] ) #----------------------------------------|A
    fsm[, 3] <- exp( psinit[2, ] ) #----------------------------------------|B
    fsm[, 4] <- exp( psinit[3, ] ) #----------------------------------------|C
    fsm[, 5] <- exp( psinit[1, ] + psinit[2, ] ) #--------------------------|AB
    fsm[, 6] <- exp( psinit[2, ] + psinit[3, ] ) #--------------------------|BC
    fsm[, 7] <- exp( psinit[1, ] + psinit[3, ] ) #--------------------------|AC
    fsm[, 8] <- exp( psinit[1, ] + psinit[2, ] + psinit[3, ] ) #------------|ABC
    ######
    # Latent state for dynamic part of model
    # tpm = transition probability matrix. All columns sum to 1.
    # after dividing by their respective column sum.
    # dim(tpm)[1] = site  
    # dim(tpm)[2] = state at time t
    # dim(tpm)[3] = state at time t-1
    ######
    tpm <- array(0, dim = c(nsite, 8, 8))
    ######
    # Latent state for dynamic part of model
    # tpm = transition probability matrix. All columns sum to 1.
    # dim(tpm)[1] = site j 
    # dim(tpm)[2] = state at time t
    # dim(tpm)[3] = state at time t-1
    ######
    # U to ...
    tpm[ , 1, 1] <- 1 #-----------------------------------------------------|U
    tpm[ , 2, 1] <- exp( gam[1, ] ) #---------------------------------------|A
    tpm[ , 3, 1] <- exp(            gam[2, ] ) #----------------------------|B
    tpm[ , 4, 1] <- exp(                       gam[3, ] ) #-----------------|C
    tpm[ , 5, 1] <- exp( gam[1, ] + gam[2, ] ) #----------------------------|AB
    tpm[ , 6, 1] <- exp( gam[1, ] +            gam[3, ] ) #-----------------|AC
    tpm[ , 7, 1] <- exp(            gam[2, ] + gam[3, ] ) #-----------------|BC
    tpm[ , 8, 1] <- exp( gam[1, ] + gam[2, ] + gam[3, ] ) # ----------------|ABC
    # A to ...
    tpm[, 1, 2] <- exp( eps[1, ] ) #----------------------------------------|U
    tpm[, 2, 2] <- 1 #------------------------------------------------------|A
    tpm[, 3, 2] <- exp( eps[1, ] + gam_one[2, 1, ] ) #----------------------|B
    tpm[, 4, 2] <- exp( eps[1, ] +                   gam_one[3, 1, ] ) #----|C
    tpm[, 5, 2] <- exp(            gam_one[2, 1, ] ) #----------------------|AB
    tpm[, 6, 2] <- exp(                              gam_one[3, 1, ] ) #----|AC
    tpm[, 7, 2] <- exp( eps[1, ] + gam_one[2, 1, ] + gam_one[3, 1, ] ) #----|BC
    tpm[, 8, 2] <- exp(            gam_one[2, 1, ] + gam_one[3, 1, ] ) #----|ABC
    # B to ...
    tpm[, 1, 3] <- exp(                   eps[2, ] ) #----------------------|U
    tpm[, 2, 3] <- exp( gam_one[1, 2, ] + eps[2, ] ) #----------------------|A
    tpm[, 3, 3] <- 1 #------------------------------------------------------|B
    tpm[, 4, 3] <- exp(                   eps[2, ] + gam_one[3, 2, ] ) #----|C
    tpm[, 5, 3] <- exp( gam_one[1, 2, ] ) #---------------------------------|AB
    tpm[, 6, 3] <- exp( gam_one[1, 2, ] + eps[2, ] + gam_one[3, 2, ] ) #----|AC
    tpm[, 7, 3] <- exp(                              gam_one[3, 2, ] ) #----|BC
    tpm[, 8, 3] <- exp( gam_one[1, 2, ] +            gam_one[3, 2, ] ) #----|ABC
    # C to ...
    tpm[, 1, 4] <- exp(                                     eps[3, ] ) #----|U
    tpm[, 2, 4] <- exp( gam_one[1, 3, ] +                   eps[3, ] ) #----|A
    tpm[, 3, 4] <- exp(                   gam_one[2, 3, ] + eps[3, ] ) #----|B
    tpm[, 4, 4] <- 1 #------------------------------------------------------|C
    tpm[, 5, 4] <- exp( gam_one[1, 3, ] + gam_one[2, 3, ] + eps[3, ] ) #----|AB
    tpm[, 6, 4] <- exp( gam_one[1, 3, ]) #----------------------------------|AC
    tpm[, 7, 4] <- exp(                   gam_one[2, 3, ] ) #---------------|BC
    tpm[, 8, 4] <- exp( gam_one[1, 3, ] + gam_one[2, 3, ] ) #---------------|ABC
    # AB to ..
    tpm[, 1, 5] <- exp( eps_one[1, 2, ] + eps_one[2, 1, ] ) #---------------|U
    tpm[, 2, 5] <- exp( eps_one[1, 2, ] ) #---------------------------------|A
    tpm[, 3, 5] <- exp(                   eps_one[2, 1, ] ) #---------------|B
    tpm[, 4, 5] <- exp( eps_one[1, 2, ] + eps_one[2, 1, ] + gam_two[3, ] ) #|C
    tpm[, 5, 5] <- 1 #------------------------------------------------------|AB
    tpm[, 6, 5] <- exp(                   eps_one[2, 1, ] + gam_two[3, ] ) #|AC
    tpm[, 7, 5] <- exp( eps_one[1, 2, ] +                   gam_two[3, ] ) #|BC
    tpm[, 8, 5] <- exp(                                     gam_two[3, ] ) #|ABC
    # AC to ...
    tpm[, 1, 6] <- exp( eps_one[1, 3, ] +                eps_one[3, 1, ] ) #|U
    tpm[, 2, 6] <- exp(                                  eps_one[3, 1, ] ) #|A
    tpm[, 3, 6] <- exp( eps_one[1, 3, ] + gam_two[2, ] + eps_one[3, 1, ] ) #|B
    tpm[, 4, 6] <- exp( eps_one[1, 3, ] ) #---------------------------------|C
    tpm[, 5, 6] <- exp(                   gam_two[2, ] + eps_one[3, 1, ] ) #|AB
    tpm[, 6, 6] <- 1 #------------------------------------------------------|AC
    tpm[, 7, 6] <- exp( eps_one[1, 3, ] + gam_two[2, ] ) #------------------|BC
    tpm[, 8, 6] <- exp(                   gam_two[2, ] ) #------------------|ABC
    # BC to ...
    tpm[, 1, 7] <- exp(                eps_one[2, 3, ] + eps_one[3, 2, ] ) #|U
    tpm[, 2, 7] <- exp( gam_two[1, ] + eps_one[2, 3, ] + eps_one[3, 2, ] ) #|A
    tpm[, 3, 7] <- exp(                                  eps_one[3, 2, ] ) #|B
    tpm[, 4, 7] <- exp(                eps_one[2, 3, ] ) #------------------|C
    tpm[, 5, 7] <- exp( gam_two[1, ] +                   eps_one[3, 2, ] ) #|AB
    tpm[, 6, 7] <- exp( gam_two[1, ] + eps_one[2, 3, ] ) #------------------|AC
    tpm[, 7, 7] <- 1 #------------------------------------------------------|BC
    tpm[, 8, 7] <- exp( gam_two[1, ] ) #------------------------------------|ABC
    # ABC to ...
    tpm[, 1, 8] <- exp( eps_two[1, ] + eps_two[2, ] + eps_two[3, ] ) #------|U
    tpm[, 2, 8] <- exp(                eps_two[2, ] + eps_two[3, ] ) #------|A
    tpm[, 3, 8] <- exp( eps_two[1, ] +                eps_two[3, ] ) #------|B
    tpm[, 4, 8] <- exp( eps_two[1, ] + eps_two[2, ] ) #---------------------|C
    tpm[, 5, 8] <- exp(                               eps_two[3, ] ) #------|AB
    tpm[, 6, 8] <- exp(                eps_two[2, ] ) #---------------------|AC 
    tpm[, 7, 8] <- exp( eps_two[1, ] ) #------------------------------------|BC
    tpm[, 8, 8] <- 1 #------------------------------------------------------|ABC
    ######
    # detection matrix (OS = observed state, TS = true state)
    # rdm = rho detection matrix. Each row sums to 1.
    # OS along rows, TS along columns
    ######
    # TS = U
    rdm <- array(0, dim = c(nsite, 8, 8))
    rdm[, 1, 1] <- 1 #-------------------------------------------------|OS = U
    rdm[, 2, 1] <- 0 #-------------------------------------------------|OS = A
    rdm[, 3, 1] <- 0 #-------------------------------------------------|OS = B
    rdm[, 4, 1] <- 0 #-------------------------------------------------|OS = C
    rdm[, 5, 1] <- 0 #-------------------------------------------------|OS = AB
    rdm[, 6, 1] <- 0 #-------------------------------------------------|OS = BC
    rdm[, 7, 1] <- 0 #-------------------------------------------------|OS = AC
    rdm[, 8, 1] <- 0 #-------------------------------------------------|OS = ABC
    # TS = A
    rdm[, 1, 2] <- 1 #-------------------------------------------------|OS = U
    rdm[, 2, 2] <- exp( rho[1, ] ) #-----------------------------------|OS = A
    rdm[, 3, 2] <- 0 #-------------------------------------------------|OS = B
    rdm[, 4, 2] <- 0 #-------------------------------------------------|OS = C
    rdm[, 5, 2] <- 0 #-------------------------------------------------|OS = AB
    rdm[, 6, 2] <- 0 #-------------------------------------------------|OS = BC
    rdm[, 7, 2] <- 0 #-------------------------------------------------|OS = AC
    rdm[, 8, 2] <- 0 #-------------------------------------------------|OS = ABC
    # TS = B
    rdm[, 1, 3] <- 1 #-------------------------------------------------|OS = U
    rdm[, 2, 3] <- 0 #-------------------------------------------------|OS = A
    rdm[, 3, 3] <- exp( rho[2, ] ) #-----------------------------------|OS = B
    rdm[, 4, 3] <- 0 #-------------------------------------------------|OS = C
    rdm[, 5, 3] <- 0 #-------------------------------------------------|OS = AB
    rdm[, 6, 3] <- 0 #-------------------------------------------------|OS = BC
    rdm[, 7, 3] <- 0 #-------------------------------------------------|OS = AC
    rdm[, 8, 3] <- 0 #-------------------------------------------------|OS = ABC
    #TS = C
    rdm[, 1, 4] <- 1 #-------------------------------------------------|OS = U
    rdm[, 2, 4] <- 0 #-------------------------------------------------|OS = A
    rdm[, 3, 4] <- 0 #-------------------------------------------------|OS = B
    rdm[, 4, 4] <- exp( rho[3, ] ) #-----------------------------------|OS = C
    rdm[, 5, 4] <- 0 #-------------------------------------------------|OS = AB
    rdm[, 6, 4] <- 0 #-------------------------------------------------|OS = BC
    rdm[, 7, 4] <- 0 #-------------------------------------------------|OS = AC
    rdm[, 8, 4] <- 0 #-------------------------------------------------|OS = ABC
    # TS = AB
    rdm[, 1, 5] <- 1 #-------------------------------------------------|OS = U
    rdm[, 2, 5] <- exp( rho[1, ] ) #-----------------------------------|OS = A
    rdm[, 3, 5] <- exp(            delta[1, ] ) #----------------------|OS = B
    rdm[, 4, 5] <- 0 #-------------------------------------------------|OS = C
    rdm[, 5, 5] <- exp( rho[1, ] + delta[1, ] ) #----------------------|OS = AB
    rdm[, 6, 5] <- 0 #-------------------------------------------------|OS = BC
    rdm[, 7, 5] <- 0 #-------------------------------------------------|OS = AC
    rdm[, 8, 5] <- 0 #-------------------------------------------------|OS = ABC
    # TS = BC
    rdm[, 1, 6] <- 1 #-------------------------------------------------|OS = U
    rdm[, 2, 6] <- 0 #-------------------------------------------------|OS = A
    rdm[, 3, 6] <- exp( rho[2, ] ) #-----------------------------------|OS = B
    rdm[, 4, 6] <- exp(            rho[3, ] ) #------------------------|OS = C
    rdm[, 5, 6] <- 0 #-------------------------------------------------|OS = AB
    rdm[, 6, 6] <- exp( rho[2, ] + rho[3, ] ) #------------------------|OS = BC
    rdm[, 7, 6] <- 0 #-------------------------------------------------|OS = AC
    rdm[, 8, 6] <- 0 #-------------------------------------------------|OS = ABC
    # TS = AC
    rdm[, 1, 7] <- 1 #-------------------------------------------------|OS = U
    rdm[, 2, 7] <- exp( rho[1, ] ) #-----------------------------------|OS = A
    rdm[, 3, 7] <- 0 #-------------------------------------------------|OS = B
    rdm[, 4, 7] <- exp(            delta[2, ] ) #----------------------|OS = C
    rdm[, 5, 7] <- 0 #-------------------------------------------------|OS = AB
    rdm[, 6, 7] <- 0 #-------------------------------------------------|OS = BC
    rdm[, 7, 7] <- exp( rho[1, ] + delta[2, ] ) #----------------------|OS = AC
    rdm[, 8, 7] <- 0 #-------------------------------------------------|OS = ABC
    # TS = ABC
    rdm[, 1, 8] <- 1 #|------------------------------------------------|OS = U
    rdm[, 2, 8] <- exp( rho[1, ] ) #-----------------------------------|OS = A
    rdm[, 3, 8] <- exp(            delta[1, ] ) #----------------------|OS = B
    rdm[, 4, 8] <- exp(                         delta[2, ] ) #---------|OS = C
    rdm[, 5, 8] <- exp( rho[1, ] + delta[1, ] ) #----------------------|OS = AB
    rdm[, 6, 8] <- exp(            delta[1, ] + delta[2, ] ) #---------|OS = BC
    rdm[, 7, 8] <- exp( rho[1, ] +              delta[2, ] ) #---------|OS = AC
    rdm[, 8, 8] <- exp( rho[1, ] + delta[1, ] + delta[2, ] ) #---------|OS = ABC

    # This is the probability of detecting each community state per site
    pstate <- array(0, dim = c(nsite, 8, nyear))
    # A likelihood matrix to feed back into the full MCMC chain lik_mat
    lik <- matrix(0, ncol = nyear, nrow = nsite)
      for(j in 1:nsite){

        # Product of this times binomial coefficient = 
        # multinomial likelihood for detection
        # This is for the first year, which uses the fsm matrix below
        # divide by the appropriate column sum (indexed by z matrix)
        # to get the probability, which is then raised to the number
        # of times each state was observed.
        pstate[j, , 1] <- (rdm[j, ,z[j, 1] ] /
                       sum(rdm[j, ,z[j, 1] ] ) )^ state_count[j, 1, ]
        
        # fsm indexed to z to get likelihood of latent state
        # at first season.
  
        lik[j,1] <- ( binco[j, 1] * prod( pstate[j, ,1] ) ) * 
                    ( fsm[j, z[j, 1] ] / sum(fsm[j, ] ) )
          for( t in 2:nyear ){
            # Product of this x binomial coefficient = 
            # multinomial likelihood for detection
            pstate[j, , t] <-  (rdm[j, , z[j, t] ]/ 
                           sum( rdm[j, , z[j, t] ]) )^ state_count[j, t, ] 
    # tpm indexed to site, state at times t, state at t-1 to get likelihood
    # for latent state
            lik[j, t] <- ( binco[j, t] * prod( pstate[j, ,t] ) ) * 
                         ( tpm[ j, z[ j, t ] , z[ j, t-1 ] ] / # numerator smax
                      sum( tpm[ j,1:nout,z[j, t-1 ] ] ) ) # denominator smax
             } # close t (year)
           } # close j (site)

    lik_mat[o,] <- as.numeric(lik) # put lik into likelihood matrix
  }

  # these are sites by seasons that sampling did not happen
  # we do not want them to contribute towards the likelihood
  to_zero <- rep(0, ncol(lik_mat))
  
    for(i in 1:length(to_zero)){
      to_zero[i] <- sum(lik_mat[,i]==0)
    }
  # these are all unobserved sites
  to_go <- which(to_zero == nrow(mm))
  # remove them from the likelihood matrix
  lik_mat <- lik_mat[,-to_go]

  # calculate cpo
  CPO <--sum(log(nrow(mm) /(apply(1/lik_mat, 2, sum))))
  # return cpo
  return(CPO)
}
    


