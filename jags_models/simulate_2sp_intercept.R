
logit <- function(x) log(x/ (1 - x))

a <- b <- d <- g <- h <- f <- l <-  rep(0,2)
a[1] <- logit(0.8)
a[2] <- logit(0.4)
b[1] <- logit(0.6)
b[2] <- logit(0.4)
d[1] <- logit(0.8)
d[2] <- logit(0.6)
g[1] <- logit(0.4)
g[2] <- logit(0.5)
h[1] <- logit(0.6)
h[2] <- logit(0.2)
f[1] <- logit(0.65)
f[2] <- logit(0.7)
l[1] <- logit(0.3)
l[2] <- logit(0.4)
set.seed(5)
a <- matrix(rnorm(4, 0, 0.8), 2, 2, TRUE)
b <- matrix(rnorm(4, 0, 0.8), 2, 2, TRUE)
d <- matrix(rnorm(4, 0, 0.8), 2, 2)
g <- matrix(rnorm(4, 0, 0.8), 2, 2)
h <- matrix(rnorm(4, 0, 0.8), 2, 2)
f <- matrix(rnorm(4, 0, 0.8), 2, 2)
l <- matrix(rnorm(4, 0, 0.8), 2, 2)
nsite <- 200
nyear <- 5
x <- rnorm(200)
x <- cbind(1, x)
# fill out the linear predictors
psia <- plogis(x %*% a[1,])
psib <- plogis(x %*% a[2,])
#
gama <- plogis(x %*% b[1,])
gamb <- plogis(x %*% b[2,])
#
epsa <- plogis(x %*% d[1,])
epsb <- plogis(x %*% d[2,])
#
gamab <- plogis(x %*%b[1,] + x %*% g[1,])
gamba <- plogis(x %*%b[2,] + x %*% g[2,])
#
epsab <- plogis(x %*% d[1,] + x %*% h[1,])
epsba <- plogis(x %*% d[2,] + x %*% h[2,])
#
pa <- plogis(x %*% f[1,])
pb <- plogis(x %*% f[2,])
#
pab <- plogis(x %*% f[1,] + x %*% l[1,])
pba <- plogis(x %*% f[2,] + x %*% l[2,])

fsm <- matrix(0, ncol = 4, nrow = nsite)
# Fill transition probs first year
fsm[,1] <- (1 - psia) * (1 - psib)
fsm[,2] <- psia * (1 - psib)
fsm[,3] <- (1 - psia) * psib
fsm[,4] <- psia * psib
tpm <- array(0, dim = c(nsite, 4, 4))
# fill transition probs across years
tpm[,1,1] <- (1 - gama) * (1 - gamb)
tpm[,2,1] <- gama * (1 - gamb)
tpm[,3,1] <- (1 - gama) * gamb
tpm[,4,1] <- gama * gamb
#
tpm[,1,2] <- epsa * (1 - gamb)
tpm[,2,2] <- (1 - epsa) * (1 - gamba)
tpm[,3,2] <- epsa * gamb
tpm[,4,2] <- (1 - epsa) * gamba
#
tpm[,1,3] <- (1 - gama) * epsb
tpm[,2,3] <- gama * epsb
tpm[,3,3] <- (1 - gamab) * (1 - epsb)
tpm[,4,3] <- gamab * (1 - epsb)
#
tpm[,1,4] <- epsab * epsba
tpm[,2,4] <- (1 - epsab) * epsba
tpm[,3,4] <- epsab * (1 - epsba)
tpm[,4,4] <- (1 - epsab) * (1 - epsba)
#
rdm <- array(0, dim = c(nsite, 4, 4))
# rho detection matrix
rdm[,1,1] <- 1
rdm[,2,1] <- 0
rdm[,3,1] <- 0
rdm[,4,1] <- 0
#
rdm[,1,2] <- (1 - pa)
rdm[,2,2] <- pa
rdm[,3,2] <- 0
rdm[,4,2] <- 0
#
rdm[,1,3] <- (1 - pb)
rdm[,2,3] <- 0
rdm[,3,3] <- pb
rdm[,4,3] <- 0
# 
rdm[,1,4] <- (1 - pab) * (1 - pba)
rdm[,2,4] <- pab * (1 - pba)
rdm[,3,4] <- (1 - pab) * pba
rdm[,4,4] <- pab * pba


nsite <- 200
nyear <- 5
nsurvey <- 4
nout <- 4
library(LaplacesDemon)

#######
# Latent state model
#######
z <- matrix(0, ncol = nyear, nrow = nsite)
for(j in 1:nsite) {
  # for first season
  z[j, 1] <- rcat(1, fsm[j,( 1:nout )] )
  for(t in 2:nyear){
    z[j, t] <- rcat(1, tpm[j,( 1:nout ) , z[ j, t-1]] )
  }
}


######
# Detection model
######
y <- array(0, dim= c(nsite, nyear, nsurvey))
for(j in 1:nsite) {
  for(ti in 1:nyear) {
    for(week in 1:nsurvey) {
      y[j, ti, week] <- rcat(1, rdm[j,( 1:nout ) , z[j, ti] ] )
    }
  }
}


library(runjags)
library(rjags)


load.module('glm')
a <- b <- d <- g <- h <- f <- l <-  rep(0,2)

inits <- function(chain){
  gen_list <- function(chain = chain){
    list( 
      z = z,
      a = matrix(rnorm(4), 2, 2),
      b = matrix(rnorm(4), 2, 2),
      d = matrix(rnorm(4), 2, 2),
      h = matrix(rnorm(4), 2, 2),
      l = matrix(rnorm(4), 2, 2),
      f = matrix(rnorm(4), 2, 2),
      g = matrix(rnorm(4), 2, 2),
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

data_list <- list(y = y, nsite = 200, nyear = 5, nout = 4, nsurvey = 4,
                  nspec = 2, x = x)

mout <- run.jags(model = "./jags_models/test_jags_state_covariate.R",
                             monitor = c("a", "b", "d", "f", "h", "g", "l"), 
                             data = data_list,
                             n.chains = 7,
                             inits = inits,
                             adapt = 100,
                             burnin = 10000,
                             sample = ceiling(10000/7),
                             thin = 5,
                             summarise = FALSE,
                             plots = FALSE,
                             method = "parallel")
a
summary(mout)[,1:4]

m2 <- as.mcmc.list(mout)
library(mcmcplots)
caterplot(m2, "g", reorder = FALSE)
points(rev(c(1:4)) ~ as.numeric(g))

