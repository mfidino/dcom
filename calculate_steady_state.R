###############################
#
# calculate_steady_state
#
# Written by Mason Fidino
#
###############################

# This function calculates the expected occupancy from the transition matrix
#  in our 3 species coyote, opossum, and raccoon example. As the best fit model
#  indicated that interactions vary as a function of urbanization, we calculate
#  the expected occupancy of each commmunity state across this metric.
#  This function requires:

# Arguments for this function
# --------------------------- #
# mm = model matrix of jags output as a matrix
#  e.g., mm <- as.matrix(jags_output, chains = TRUE) via coda package
#  columns of this matrix represent different parameters while rows are
#  posterior simulations

# data_list = the same data_list supplied to JAGS to fit the DCOM
#   This was generated in fit_softmax_model.R
# ncores = number of cores to run this function on in parallel

# Note that the paramaters are given the same names as in the manuscript.
#  For example, 'a' parameters are for initial occupancy.

# This function returns a list of length MCMC samples (i.e., nrow(mm)). Each
#  element in this list is an 81 x 8 matrix. Rows correspond to the urban
#  PCA that range from -4 to 4 by steps of 0.1 while columns are the 
#  stationary occupancy probability of each state. In our example these states
#  are: No species, coyote, opossum, raccoon, coyote & opossum, coyote & 
#  raccoon, opossum & raccoon, all three species. Each element in the list
#  corresponds to these predictions for that particular step in the MCMC chain
#  from the DCOM.

# Packages required:
#  doParallel
#  markovchain

mm <- mm[,2:51]

calculate_steady_state <- function(mm = NULL, data_list = NULL, ncores = 7){
  
  # load packages
  library(doParallel)
  library(markovchain)
  
  # urban PCA metric
  xp <- seq(-4, 4, 0.1)
  
  # x by 2 matrix. 1 for Intercept value.
  xp <- cbind(1, xp)
  # ensure mm (MCMC model output) is a matrix
  if(!is.matrix(mm)) {
    stop("The object mm must be a matrix")
  }
  
  # Assigning all the objects in data_list so I don't have to
  # index them the whole time.
  for(i in 1:length(data_list)){
    assign(names(data_list)[i], data_list[[i]])
  }
  
  # not really the number of sites, just nrow of urban PCA for predictions
  nsite <- nrow(xp)
  
  # compute steady states 
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  n <- nrow(mm)
 #og_mm2 <- mm

 
 #mm <- t(as.matrix(apply(mm, 2, median)  ))
  
  smat <- foreach(o = 1:nrow(mm), .packages = c('markovchain')) %dopar% { # going through each step of the mcmc chain
                                                  
    steady_mat <- array(0, dim = c(nrow(xp), 8))
    # initial occupancy
    a <- matrix(mm[o,grep("a", colnames(mm))], ncol = ncov_psi, nrow = nspec)
    # colonization 
    b <- matrix(mm[o,grep("b\\[", colnames(mm))], ncol = ncov_gam, nrow = nspec)
    # colonization 
    # extinction
    d <- matrix(mm[o,grep("d\\[", colnames(mm))], ncol = ncov_eps, nrow = nspec)
    # colonization inxs
    g <- matrix(mm[o,grep("g", colnames(mm))], ncol = ncov_pi, nrow = n_inxs) 
    if(sum(is.na(g)>0)) g[is.na(g)] <- 0 # make 0 if not in the model
    # extinciton inxs
    h <- matrix(mm[grep("h", colnames(mm))], ncol = ncov_tau, nrow = n_inxs)
    if(sum(is.na(h)>0)) h[is.na(h)] <- 0 # make 0 if not in the model
    

    # occupancy linear predictor
    psinit  <- a %*% t(xp)
    # colonization linear predictor
    gam <-     b %*% t(xp)
    # extinction linear predictor
    eps <-     d %*% t(xp)
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
    
    gam_one <- eps_one <- array(0, dim = c(nspec, nspec, nsite))
      for(k in 1:n_inxs){
        # inxs on colonization linear predictor
        gam_one[rows_vec[k], cols_vec[k], ] <-  g[k, ]
        # inxs on extinction linear predictor
        eps_one[rows_vec[k], cols_vec[k], ] <-  h[k,]  
      }
    gam_two <- matrix(0, ncol = nspec, nrow = nsite)
    gam_two[,1] <- b[1,] %*% t(xp) + gam_one[1,2, ] + gam_one[1,3,]
    gam_two[,2] <- b[2,] %*% t(xp) + gam_one[2,1, ] + gam_one[2,3,]
    gam_two[,3] <- b[3,] %*% t(xp) + gam_one[3,1, ] + gam_one[3,2,]
    
    eps_two <- matrix(0, ncol = nspec, nrow = nsite)
    eps_two[,1] <- d[1,] %*% t(xp) + eps_one[1,2, ] + eps_one[1,3,]
    eps_two[,2] <- d[2,] %*% t(xp) + eps_one[2,1, ] + eps_one[2,3,]
    eps_two[,3] <- d[3,] %*% t(xp) + eps_one[3,1, ] + eps_one[3,2,]
    
    for(k in 1:n_inxs){
      # inxs on colonization linear predictor
      gam_one[rows_vec[k], cols_vec[k], ] <- 
        t(b[rows_vec[k],] %*% t(xp)) +  g[k, ]
      # inxs on extinction linear predictor
      eps_one[rows_vec[k], cols_vec[k], ] <- 
        t(d[rows_vec[k],] %*% t(xp)) +   h[k,]  
    }

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
    # dim(tpm)[1] = site  
    # dim(tpm)[2] = state at time t
    # dim(tpm)[3] = state at time t-1
    ######
    # U to ...
    tpm[, 1, 1] <- 1 #---------------------------------------------------------|U
    tpm[, 2, 1] <- exp( gam[1, ] ) #------------------------------------------|A
    tpm[, 3, 1] <- exp(             gam[2, ] ) #------------------------------|B
    tpm[, 4, 1] <- exp(                         gam[3, ] ) #------------------|C
    tpm[, 5, 1] <- exp( gam[1, ] + gam[2, ] ) #------------------------------|AB
    tpm[, 6, 1] <- exp( gam[1, ] +             gam[3, ] ) #------------------|AC
    tpm[, 7, 1] <- exp(             gam[2, ] + gam[3, ] ) #------------------|BC
    tpm[, 8, 1] <- exp( gam[1, ] + gam[2, ] + gam[3, ] ) # -----------------|ABC
    # A to ...
    tpm[, 1, 2] <- exp( eps[1, ] ) #------------------------------------------|U
    tpm[, 2, 2] <- 1 #---------------------------------------------------------|A
    tpm[, 3, 2] <- exp( eps[1, ] + gam_one[2, 1, ] ) #-----------------------|B
    tpm[, 4, 2] <- exp( eps[1, ] +                    gam_one[3, 1, ] ) #----|C
    tpm[, 5, 2] <- exp(             gam_one[2, 1, ] ) #-----------------------|AB
    tpm[, 6, 2] <- exp(                                gam_one[3, 1, ] ) #----|AC
    tpm[, 7, 2] <- exp( eps[1, ] + gam_one[2, 1, ] + gam_one[3, 1, ] ) #----|BC
    tpm[, 8, 2] <- exp(             gam_one[2, 1, ] + gam_one[3, 1, ] ) #----|ABC
    # B to ...
    tpm[, 1, 3] <- exp(                    eps[2, ] ) #-----------------------|U
    tpm[, 2, 3] <- exp( gam_one[1, 2, ] + eps[2, ] ) #-----------------------|A
    tpm[, 3, 3] <- 1 #---------------------------------------------------------|B
    tpm[, 4, 3] <- exp(                    eps[2, ] + gam_one[3, 2, ] ) #----|C
    tpm[, 5, 3] <- exp( gam_one[1, 2, ] ) #-----------------------------------|AB
    tpm[, 6, 3] <- exp( gam_one[1, 2, ] + eps[2, ] + gam_one[3, 2, ] ) #----|AC
    tpm[, 7, 3] <- exp(                                gam_one[3, 2, ] ) #----|BC
    tpm[, 8, 3] <- exp( gam_one[1, 2, ] +             gam_one[3, 2, ] ) #----|ABC
    # C to ...
    tpm[, 1, 4] <- exp(                                       eps[3, ] ) #----|U
    tpm[, 2, 4] <- exp( gam_one[1, 3, ] +                    eps[3, ] ) #----|A
    tpm[, 3, 4] <- exp(                    gam_one[2, 3, ] + eps[3, ] ) #----|B
    tpm[, 4, 4] <- 1 #---------------------------------------------------------|C
    tpm[, 5, 4] <- exp( gam_one[1, 3, ] + gam_one[2, 3, ] + eps[3, ] ) #----|AB
    tpm[, 6, 4] <- exp( gam_one[1, 3, ]) #------------------------------------|AC
    tpm[, 7, 4] <- exp(                    gam_one[2, 3, ] ) #----------------|BC
    tpm[, 8, 4] <- exp( gam_one[1, 3, ] + gam_one[2, 3, ] ) #----------------|ABC
    # AB to ..
    tpm[, 1, 5] <- exp( eps_one[1, 2, ] + eps_one[2, 1, ] ) #----------------|U
    tpm[, 2, 5] <- exp( eps_one[1, 2, ] ) #-----------------------------------|A
    tpm[, 3, 5] <- exp(                    eps_one[2, 1, ] ) #----------------|B
    tpm[, 4, 5] <- exp( eps_one[1, 2, ] + eps_one[2, 1, ] + gam_two[,3 ] ) #|C
    tpm[, 5, 5] <- 1 #---------------------------------------------------------|AB
    tpm[, 6, 5] <- exp(                    eps_one[2, 1, ] + gam_two[,3 ] ) #|AC
    tpm[, 7, 5] <- exp( eps_one[1, 2, ] +                    gam_two[,3 ] ) #|BC
    tpm[, 8, 5] <- exp(                                       gam_two[,3 ] ) #|ABC
    # AC to ...
    tpm[, 1, 6] <- exp( eps_one[1, 3, ] +                 eps_one[3, 1, ] ) #|U
    tpm[, 2, 6] <- exp(                                    eps_one[3, 1, ] ) #|A
    tpm[, 3, 6] <- exp( eps_one[1, 3, ] + gam_two[,2 ] + eps_one[3, 1, ] ) #|B
    tpm[, 4, 6] <- exp( eps_one[1, 3, ] ) #-----------------------------------|C
    tpm[, 5, 6] <- exp(                    gam_two[,2 ] + eps_one[3, 1, ] ) #|AB
    tpm[, 6, 6] <- 1 #---------------------------------------------------------|AC
    tpm[, 7, 6] <- exp( eps_one[1, 3, ] + gam_two[,2 ] ) #-------------------|BC
    tpm[, 8, 6] <- exp(                    gam_two[,2 ] ) #-------------------|ABC
    # BC to ...
    tpm[, 1, 7] <- exp(                 eps_one[2, 3, ] + eps_one[3, 2, ] ) #|U
    tpm[, 2, 7] <- exp( gam_two[,1 ] + eps_one[2, 3, ] + eps_one[3, 2, ] ) #|A
    tpm[, 3, 7] <- exp(                                    eps_one[3, 2, ] ) #|B
    tpm[, 4, 7] <- exp(                 eps_one[2, 3, ] ) #-------------------|C
    tpm[, 5, 7] <- exp( gam_two[,1 ] +                    eps_one[3, 2, ] ) #|AB
    tpm[, 6, 7] <- exp( gam_two[,1 ] + eps_one[2, 3, ] ) #-------------------|AC
    tpm[, 7, 7] <- 1 #---------------------------------------------------------|BC
    tpm[, 8, 7] <- exp( gam_two[,1 ] ) #--------------------------------------|ABC
    # ABC to ...
    tpm[, 1, 8] <- exp( eps_two[,1 ] + eps_two[,2 ] + eps_two[,3 ] ) #------|U
    tpm[, 2, 8] <- exp(                 eps_two[,2 ] + eps_two[,3 ] ) #------|A
    tpm[, 3, 8] <- exp( eps_two[,1 ] +                 eps_two[,3 ] ) #------|B
    tpm[, 4, 8] <- exp( eps_two[,1 ] + eps_two[,2 ] ) #----------------------|C
    tpm[, 5, 8] <- exp(                                 eps_two[,3 ] ) #------|AB
    tpm[, 6, 8] <- exp(                 eps_two[,2 ] ) #----------------------|AC 
    tpm[, 7, 8] <- exp( eps_two[,1 ] ) #--------------------------------------|BC
    tpm[, 8, 8] <- 1 #--------------------------------------------------------|ABC

   
   # convert to probability, Divide by column sum of each
    for(i in 1:nrow(xp)){
      tpm[i,,] <- apply(tpm[i,,], 2, function(x) x/ sum(x))
      mark <- 
        new("markovchain", states = c("U", "A", "B", "C", "AB", "AC", "BC", "ABC"),
                  byrow = TRUE, t(tpm[i,,]))
      steady_mat[i,]  <-   steadyStates(mark)
    }
    steady_mat
  }
    stopCluster(cl)
    return(smat)
    

}

