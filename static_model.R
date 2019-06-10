model{
  ######
  # Detection model
  ######
  for(site in 1:nsite) {
      for(week in 1:nsurvey) {
        y[site, week] ~ dcat( rdm[site, ( 1:nout ) , z[site] ] )
      }
  }
  #######
  # Latent state model
  #######
  for(j in 1:nsite) {
    z[j, 1] ~ dcat( fsm[j, ( 1:nout )] )
  }
  for( j in 1:nsite ) {
    ######
    # Fill in all of the transition probabilities
    ######
    ######
    # Latent state
    ######
    #  Probabilities for each state
    fsm[j, 1] <- 1 #--------------------------------------------------------|U
    fsm[j, 2] <- exp( psi_one[1, j] ) #-------------------------------------|A
    fsm[j, 3] <- exp( psi_one[2, j] ) #-------------------------------------|B
    fsm[j, 4] <- exp( psi_one[3, j] ) #-------------------------------------|C
    fsm[j, 5] <- exp( psi_two[1, j] ) #-------------------------------------|AB
    fsm[j, 6] <- exp( psi_two[2, j] ) #-------------------------------------|BC
    fsm[j, 7] <- exp( psi_two[2, j] ) #-------------------------------------|AC
    fsm[j, 8] <- exp( psi_three[j]  ) #-------------------------------------|ABC
    for(k in 1:nsurvey){
    ######
    # detection matrix (OS = observed state, TS = true state)
    # rdm = rho detection matrix. Each row sums to 1.
    # OS along rows, TS along columns
    ######
    # TS = U
    rdm[j, k, 1, 1] <- 1 #----------------------------------------------------|OS = U
    rdm[j, k, 2, 1] <- 0 #----------------------------------------------------|OS = A
    rdm[j, k, 3, 1] <- 0 #----------------------------------------------------|OS = B
    rdm[j, k, 4, 1] <- 0 #----------------------------------------------------|OS = C
    rdm[j, k, 5, 1] <- 0 #----------------------------------------------------|OS = AB
    rdm[j, k, 6, 1] <- 0 #----------------------------------------------------|OS = BC
    rdm[j, k, 7, 1] <- 0 #----------------------------------------------------|OS = AC
    rdm[j, k, 8, 1] <- 0 #----------------------------------------------------|OS = ABC
    # TS = A
    rdm[j, k, 1, 2] <- 1 #----------------------------------------------------|OS = U
    rdm[j, k, 2, 2] <- exp( rho[1, j] ) #-------------------------------------|OS = A
    rdm[j, k, 3, 2] <- 0 #----------------------------------------------------|OS = B
    rdm[j, k, 4, 2] <- 0 #----------------------------------------------------|OS = C
    rdm[j, k, 5, 2] <- 0 #----------------------------------------------------|OS = AB
    rdm[j, k, 6, 2] <- 0 #----------------------------------------------------|OS = BC
    rdm[j, k, 7, 2] <- 0 #----------------------------------------------------|OS = AC
    rdm[j, k, 8, 2] <- 0 #----------------------------------------------------|OS = ABC
    # TS = B
    rdm[j, k, 1, 3] <- 1 #----------------------------------------------------|OS = U
    rdm[j, k, 2, 3] <- 0 #----------------------------------------------------|OS = A
    rdm[j, k, 3, 3] <- exp( rho[2, j] ) #-------------------------------------|OS = B
    rdm[j, k, 4, 3] <- 0 #----------------------------------------------------|OS = C
    rdm[j, k, 5, 3] <- 0 #----------------------------------------------------|OS = AB
    rdm[j, k, 6, 3] <- 0 #----------------------------------------------------|OS = BC
    rdm[j, k, 7, 3] <- 0 #----------------------------------------------------|OS = AC
    rdm[j, k, 8, 3] <- 0 #----------------------------------------------------|OS = ABC
    #TS = C
    rdm[j, k, 1, 4] <- 1 #----------------------------------------------------|OS = U
    rdm[j, k, 2, 4] <- 0 #----------------------------------------------------|OS = A
    rdm[j, k, 3, 4] <- 0 #----------------------------------------------------|OS = B
    rdm[j, k, 4, 4] <- exp( rho[3, j] ) #-------------------------------------|OS = C
    rdm[j, k, 5, 4] <- 0 #----------------------------------------------------|OS = AB
    rdm[j, k, 6, 4] <- 0 #----------------------------------------------------|OS = BC
    rdm[j, k, 7, 4] <- 0 #----------------------------------------------------|OS = AC
    rdm[j, k, 8, 4] <- 0 #----------------------------------------------------|OS = ABC
    # TS = AB
    rdm[j, k, 1, 5] <- 1 #----------------------------------------------------|OS = U
    rdm[j, k, 2, 5] <- exp( rho[1, j] ) #-------------------------------------|OS = A
    rdm[j, k, 3, 5] <- exp(               delta[2, j] ) #---------------------|OS = B
    rdm[j, k, 4, 5] <- 0 #----------------------------------------------------|OS = C
    rdm[j, k, 5, 5] <- exp( delta[1, j] + delta[2, j] ) #---------------------|OS = AB
    rdm[j, k, 6, 5] <- 0 #----------------------------------------------------|OS = BC
    rdm[j, k, 7, 5] <- 0 #----------------------------------------------------|OS = AC
    rdm[j, k, 8, 5] <- 0 #----------------------------------------------------|OS = ABC
    # TS = BC
    rdm[j, k, 1, 6] <- 1 #----------------------------------------------------|OS = U
    rdm[j, k, 2, 6] <- 0 #----------------------------------------------------|OS = A
    rdm[j, k, 3, 6] <- exp( rho[2, j] ) #-------------------------------------|OS = B
    rdm[j, k, 4, 6] <- exp(             rho[3, j] ) #-------------------------|OS = C
    rdm[j, k, 5, 6] <- 0 #----------------------------------------------------|OS = AB
    rdm[j, k, 6, 6] <- exp( rho[2, j] + rho[3, j] ) #-------------------------|OS = BC
    rdm[j, k, 7, 6] <- 0 #----------------------------------------------------|OS = AC
    rdm[j, k, 8, 6] <- 0 #----------------------------------------------------|OS = ABC
    # TS = AC
    rdm[j, k, 1, 7] <- 1 #----------------------------------------------------|OS = U
    rdm[j, k, 2, 7] <- exp( rho[1, j] ) #-------------------------------------|OS = A
    rdm[j, k, 3, 7] <- 0 #----------------------------------------------------|OS = B
    rdm[j, k, 4, 7] <- exp(             delta[2, j] ) #-----------------------|OS = C
    rdm[j, k, 5, 7] <- 0 #----------------------------------------------------|OS = AB
    rdm[j, k, 6, 7] <- 0 #----------------------------------------------------|OS = BC
    rdm[j, k, 7, 7] <- exp( rho[1, j] + delta[2, j] ) #-----------------------|OS = AC
    rdm[j, k, 8, 7] <- 0 #----------------------------------------------------|OS = ABC
    # TS = ABC
    rdm[j, k, 1, 8] <- 1 #|---------------------------------------------------|OS = U
    rdm[j, k, 2, 8] <- exp( rho[1, j] ) #-------------------------------------|OS = A
    rdm[j, k, 3, 8] <- exp(             delta[1, j] ) #-----------------------|OS = B
    rdm[j, k, 4, 8] <- exp(                           delta[2, j] ) #---------|OS = C
    rdm[j, k, 5, 8] <- exp( rho[1, j] + delta[1, j] ) #-----------------------|OS = AB
    rdm[j, k, 6, 8] <- exp(             delta[1, j] + delta[2, j] ) #---------|OS = BC
    rdm[j, k, 7, 8] <- exp( rho[1, j] +               delta[2, j] ) #---------|OS = AC
    rdm[j, k, 8, 8] <- exp( rho[1, j] + delta[1, j] + delta[2, j] ) #---------|OS = ABC
    }
    ######
    # Occupancy linear predictors...
    ######
    for( i in 1:nspec ) {
      # ...for states A, B, and C (in that order)
      psi_one[i, j]  <- inprod( a[i, ], psi_cov[j, ] )
    }
    # ...for states AB, BC, and AC (in that order)
    psi_two[1, j] <- psi_one[1, j] + psi_one[2, j] + inprod( b[1, ], psi_inxs_cov[j, ] )
    psi_two[2, j] <- psi_one[2, j] + psi_one[3, j] + inprod( b[2, ], psi_inxs_cov[j, ] )
    psi_two[3, j] <- psi_one[1, j] + psi_one[3, j] + inprod( b[3, ], psi_inxs_cov[j, ] )
    # ...for state ABC
    psi_three[j] <- sum( psi_one[ , j] ) + inprod( b[1, ], psi_inxs_cov[j, ] ) +
     inprod( b[2, ], psi_inxs_cov[j, ] ) + inprod( b[3, ], psi_inxs_cov[j, ] )
    #######
    # Detection linear predictors...
    #######
    for(i in 1:nspec){
      for(k in 1:)
      # ...for when the TRUE state is either A, B, or C (in that order)
      rho_one[i, j] <- inprod(d[i,], rho_cov[])
    }
    # coyote influence on detection for opossum and raccoon
    delta[1, j]  <- inprod(f[2, ], rho_cov[j, ] ) + inprod( l[1, ], delta_cov[j, ] )
    delta[2, j]  <- inprod(f[3, ], rho_cov[j, ] ) + inprod( l[2, ], delta_cov[j, ] )
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
      #b0[i, gamp] ~ dlogis(0, 1)
    }
    # Extinction
    for( epsp in 1:ncov_eps ){
      d[i, epsp] ~ dlogis(0, 1)
      #d0[i, epsp] ~ dlogis(0, 1)
    }
    # Detection
    for( rhop in 1:ncov_rho ){
      f[i, rhop] ~ dlogis(0, 1)
    }
  }
  for( k in 1:n_inxs ){
    # Inxs on colonization
    for( pip in 1:ncov_pi ){
      g[k, pip] ~ dlogis(0, 1)
    }
    # Inxs on extinction
    for( taup in 1:ncov_tau ){
      h[k, taup] ~ dlogis(0, 1)
    }
  }
  # Inxs on detection
  for( delp in 1:ncov_delta ){
    l[1, delp] ~ dlogis(0, 1)
    l[2, delp] ~ dlogis(0, 1)
  }
}
