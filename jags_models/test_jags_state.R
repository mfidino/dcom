model{
######
# Detection model
######
for(j in 1:nsite) {
  for(ti in 1:nyear) {
    for(week in 1:nsurvey) {
      y[j, ti, week] ~ dcat( rdm[( 1:nout ) , z[j, ti] ] )
    }
  }
}
#######
# Latent state model
#######
for(j in 1:nsite) {
  # for first season
  z[j, 1] ~ dcat( fsm[ ( 1:nout )] )
  for(t in 2:nyear){
    z[j, t] ~ dcat( tpm[( 1:nout ) , z[ j, t-1]] )
  }
}
# Fill transition probs first year
fsm[1] <- (1 - psia) * (1 - psib)
fsm[2] <- psia * (1 - psib)
fsm[3] <- (1 - psia) * psib
fsm[4] <- psia * psib
# fill transition probs across years
tpm[1,1] <- (1 - gama) * (1 - gamb)
tpm[2,1] <- gama * (1 - gamb)
tpm[3,1] <- (1 - gama) * gamb
tpm[4,1] <- gama * gamb
#
tpm[1,2] <- epsa * (1 - gamb)
tpm[2,2] <- (1 - epsa) * (1 - gamba)
tpm[3,2] <- epsa * gamb
tpm[4,2] <- (1 - epsa) * gamba
#
tpm[1,3] <- (1 - gama) * epsb
tpm[2,3] <- gama * epsb
tpm[3,3] <- (1 - gamab) * (1 - epsb)
tpm[4,3] <- gamab * (1 - epsb)
#
tpm[1,4] <- epsab * epsba
tpm[2,4] <- (1 - epsab) * epsba
tpm[3,4] <- epsab * (1 - epsba)
tpm[4,4] <- (1 - epsab) * (1 - epsba)
#
# rho detection matrix
rdm[1,1] <- 1
rdm[2,1] <- 0
rdm[3,1] <- 0
rdm[4,1] <- 0
#
rdm[1,2] <- (1 - pa)
rdm[2,2] <- pa
rdm[3,2] <- 0
rdm[4,2] <- 0
#
rdm[1,3] <- (1 - pb)
rdm[2,3] <- 0
rdm[3,3] <- pb
rdm[4,3] <- 0
# 
rdm[1,4] <- (1 - pab) * (1 - pba)
rdm[2,4] <- pab * (1 - pba)
rdm[3,4] <- (1 - pab) * pba
rdm[4,4] <- pab * pba
# fill out the linear predictors
psia <- ilogit(a[1])
psib <- ilogit(a[2])
#
gama <- ilogit(b[1])
gamb <- ilogit(b[2])
#
epsa <- ilogit(d[1])
epsb <- ilogit(d[2])
#
gamab <- ilogit(b[1] + g[1])
gamba <- ilogit(b[2] + g[2])
#
epsab <- ilogit(d[1] + h[1])
epsba <- ilogit(d[2] + h[2])
#
pa <- ilogit(f[1])
pb <- ilogit(f[2])
#
pab <- ilogit(f[1] + l[1])
pba <- ilogit(f[2] + l[2])
# priors
for(i in 1:nspec){
  a[i] ~ dlogis(0,1)
  b[i] ~ dlogis(0,1)
  d[i] ~ dlogis(0,1)
  g[i] ~ dlogis(0,1)
  h[i] ~ dlogis(0,1)
  f[i] ~ dlogis(0,1)
  l[i] ~ dlogis(0,1)
}
}