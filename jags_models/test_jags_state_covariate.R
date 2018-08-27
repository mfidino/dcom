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
for(j in 1:nsite){
# Fill transition probs first year
fsm[j,1] <- (1 - psia[j]) * (1 - psib[j])
fsm[j,2] <- psia[j] * (1 - psib[j])
fsm[j,3] <- (1 - psia[j]) * psib[j]
fsm[j,4] <- psia[j] * psib[j]
# fill transition probs across years
tpm[j,1,1] <- (1 - gama[j]) * (1 - gamb[j])
tpm[j,2,1] <- gama[j] * (1 - gamb[j])
tpm[j,3,1] <- (1 - gama[j]) * gamb[j]
tpm[j,4,1] <- gama[j] * gamb[j]
#
tpm[j,1,2] <- epsa[j] * (1 - gamb[j])
tpm[j,2,2] <- (1 - epsa[j]) * (1 - gamba[j])
tpm[j,3,2] <- epsa[j] * gamb[j]
tpm[j,4,2] <- (1 - epsa[j]) * gamba[j]
#
tpm[j,1,3] <- (1 - gama[j]) * epsb[j]
tpm[j,2,3] <- gama[j] * epsb[j]
tpm[j,3,3] <- (1 - gamab[j]) * (1 - epsb[j])
tpm[j,4,3] <- gamab[j] * (1 - epsb[j])
#
tpm[j,1,4] <- epsab[j] * epsba[j]
tpm[j,2,4] <- (1 - epsab[j]) * epsba[j]
tpm[j,3,4] <- epsab[j] * (1 - epsba[j])
tpm[j,4,4] <- (1 - epsab[j]) * (1 - epsba[j])
#
# rho detection matrix
rdm[j,1,1] <- 1
rdm[j,2,1] <- 0
rdm[j,3,1] <- 0
rdm[j,4,1] <- 0
#
rdm[j,1,2] <- (1 - pa[j])
rdm[j,2,2] <- pa[j]
rdm[j,3,2] <- 0
rdm[j,4,2] <- 0
#
rdm[j,1,3] <- (1 - pb[j])
rdm[j,2,3] <- 0
rdm[j,3,3] <- pb[j]
rdm[j,4,3] <- 0
# 
rdm[j,1,4] <- (1 - pab[j]) * (1 - pba[j])
rdm[j,2,4] <- pab[j] * (1 - pba[j])
rdm[j,3,4] <- (1 - pab[j]) * pba[j]
rdm[j,4,4] <- pab[j] * pba[j]
# fill out the linear predictors
psia[j] <- ilogit( inprod(a[1,], x[j,]) )
psib[j] <- ilogit( inprod(a[2,], x[j,]) )
#
gama[j] <- ilogit( inprod(b[1,], x[j,]) )
gamb[j] <- ilogit( inprod(b[2,], x[j,]) )
#
epsa[j] <- ilogit( inprod(d[1,], x[j,]) )
epsb[j] <- ilogit( inprod(d[2,], x[j,]) )
#
gamab[j] <- ilogit( inprod(b[1,], x[j,] ) + inprod(g[1,], x[j,]) )
gamba[j] <- ilogit( inprod(b[2,], x[j,] ) + inprod(g[2,], x[j,]) )
#
epsab[j] <- ilogit( inprod(d[1,], x[j,] ) + inprod(h[1,], x[j,]) )
epsba[j] <- ilogit( inprod(d[2,], x[j,] ) + inprod(h[2,], x[j,]) )
#
pa[j] <- ilogit( inprod(f[1,], x[j,]) )
pb[j] <- ilogit( inprod(f[2,], x[j,]) )
#
pab[j] <- ilogit( inprod(f[1,], x[j,] ) + inprod(l[1,], x[j,]) )
pba[j] <- ilogit( inprod(f[2,], x[j,] ) + inprod(l[2,], x[j,]) )
}
# priors
for(i in 1:nspec){
  for(par in 1:2){
  a[i,par] ~ dlogis(0,1)
  b[i,par] ~ dlogis(0,1)
  d[i,par] ~ dlogis(0,1)
  g[i,par] ~ dlogis(0,1)
  h[i,par] ~ dlogis(0,1)
  f[i,par] ~ dlogis(0,1)
  l[i,par] ~ dlogis(0,1)
  }
}
}