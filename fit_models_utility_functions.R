##############################
#
# fit models utility functions
#
# written by Mason Fidino
#
##############################

# generates initial values for mcmc chains for all parameters
# for inxs model
inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = sp_inits,
      a = matrix(rnorm(data_list$nspec * data_list$ncov_psi),
                 ncol = data_list$ncov_psi, nrow = data_list$nspec),
      b = matrix(rnorm(data_list$nspec * data_list$ncov_gam),
                 ncol = data_list$ncov_gam, nrow = data_list$nspec),
      d = matrix(rnorm(data_list$nspec * data_list$ncov_eps),
                 ncol = data_list$ncov_eps, nrow = data_list$nspec),
      f = matrix(rnorm(data_list$nspec * data_list$ncov_rho),
                 ncol = data_list$ncov_rho, nrow = data_list$nspec),
      g = matrix(rnorm(data_list$n_inxs * data_list$ncov_pi),
                 ncol = data_list$ncov_pi, nrow = data_list$n_inxs),
      h = matrix(rnorm(data_list$n_inxs * data_list$ncov_tau),
                 ncol = data_list$ncov_tau, nrow = data_list$n_inxs),
      l = matrix(rnorm((data_list$nspec - 1) * data_list$ncov_delta),
                 ncol = data_list$ncov_delta, nrow = data_list$nspec-1),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}

# generates initial values for mcmc chains for all parameters
# for softmax model without interactions
inits_no_inxs <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = sp_inits,
      a = matrix(rnorm(data_list$nspec * data_list$ncov_psi),
                 ncol = data_list$ncov_psi, nrow = data_list$nspec),
      b = matrix(rnorm(data_list$nspec * data_list$ncov_gam),
                 ncol = data_list$ncov_gam, nrow = data_list$nspec),
      d = matrix(rnorm(data_list$nspec * data_list$ncov_eps),
                 ncol = data_list$ncov_eps, nrow = data_list$nspec),
      f = matrix(rnorm(data_list$nspec * data_list$ncov_rho),
                 ncol = data_list$ncov_rho, nrow = data_list$nspec),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Marsaglia-Multicarry",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}


# a function used to index interactive parameters in the model.
# In the model we put the inxs in a matrix
# these matrices get filled in this order:
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
make_inxs <- function(nspec = NULL){
  ns <- nspec
  col_l <- row_l <- col_u <- row_u <-numeric((ns * ns - ns)/2)
  for(i in 1:(ns-1)){
    # algorithm to fill lower columns
    col_l[min(which(col_l ==0)):(min(which(col_l ==0)) + 
                                   ns-1-i)] <-  rep(i, ns-i)
    # algorithm to fill lower rows
    row_l[min(which(row_l == 0)):(min(which(row_l==0))+ 
                                    ns-1-i)] <- (i+1):ns
    # algorithm to fill upper rows
    row_u[min(which(row_u == 0)):(min(which(row_u==0))+ 
                                    ns-(ns+1)+i)] <- 1:i
    # algorithm to fill upper columns
    col_u[min(which(col_u == 0)):(min(which(col_u==0))+ 
                                    ns-(ns+1)+i)] <- rep(i + 1, i)
  }
  row_vec <- c(row_l, row_u)
  col_vec <- c(col_l, col_u)
  return(list(rows_vec = row_vec,
              cols_vec = col_vec))
}

