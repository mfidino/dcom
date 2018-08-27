

# This function generates the predictions on the probability scale
# for colonization and extinction.

# mat:       The model matrix from the jags model
# sp:        Vector of the two species the plot should be about. The first number
#              is the species of interest while the second number is the species
#              we are conditioning on being at the site.
# p1:        Parameter name for non-interactive effects. character
# p2:        Parameter name for interaction effects. character
# xp:        The covariates used in the analyis

# returns a named list. 'w' is the rate of interest with the other species 
# while 'wo' is the rate without the other species.

make_preds <- function(mat, sp = NULL, p1 ,p2, xp = NULL){
  # get mcmc matrix for non inxs
  m1 <- mat[, grep(paste0(p1,"\\[",sp[1]), colnames(mat))]
  # make rows and col vec
  j1 <- matrix(c(make_inxs(3)$cols_vec, make_inxs(3)$rows_vec),
    ncol = 2, nrow = 6)
  # paste together each row and collapse
  j1 <- apply(j1, 1, paste0, collapse = "")
  # # the reverse of sp is the parameter we are looking for
  m2_ind <- which(j1 == paste0(rev(sp), collapse = ""))
  # get mcmc matrix for inxs
  m2 <- mat[,grep(paste0(p2, "\\[", m2_ind), colnames(mat))]
  # calculating the occupancy rate at equilibrium
  without <- plogis(m1 %*% t(xp))
  without <- apply(without, 2, quantile, probs = c(0.025,0.5,0.975))
  #
  with <- plogis(m1 %*% t(xp) + m2 %*% t(xp[,1]))
  with <- apply(with, 2, quantile, probs = c(0.025,0.5,0.975))
  #
  return(list(wo = without, w = with))
}


# this makes the axes for the sub-plots of the colonization and extinction
#  figures.
make_axes_ce <- function(){
  axis(1, at=seq(-4,4, 1), labels=F, tck=-.025)
  
  mtext(text = sprintf("%.0f",seq(-4,4, 2)), 1, 
    line = 1.6, at = seq(-4,4,2), cex = 1.3)
  
  axis(2, at = seq(0,1,0.2), labels = F, tck = -0.025)
  
  mtext(text = sprintf("%.1f", seq(0,1,0.2)), 2,
    line = 0.95, las = 2, at = seq(0,1,0.2), cex = 1.3)
  
}


# makes the axes for each subplot
# x: Put numbers on x axis? requires Boolean
# y: Put numbers on y axis? requires Boolean
axes_occ_plot <- function(x = NULL, y = NULL){
  axis(1, at=seq(-4,4, 1), labels=F, tck=-.035/2)
  axis(1, at=seq(-4,4, 2), labels=F, tck=-.035)
  
  if(x){
    mtext(text = sprintf("%.0f",seq(-4,4, 2)), 1, 
      line = 0.3, at = seq(-4,4,2), cex = 0.7)
  }
  axis(2, at = seq(0,1,0.25), labels = F, tck = -0.035)
  axis(2, at = seq(0,1,0.125), labels = F, tck = -0.035/2)
  
  if(y){
    mtext(text = sprintf("%.2f", seq(0,1,0.25)), 2,
      line = 0.4, las = 2, at = seq(0,1,0.25), cex = 0.7)
  }
}


# A function to generate each sub plot for the occupancy plot.
# This requires a good number of arguments. 

# steady:    This is the object 'steady_CI' in plotting_script.R. It is a 
#              81 x 3 x 8 array (URB covariate x 0.025, 0.5, 0.975 quantile x 
#              eight community states.) and contains the expected occupancy rate
#              of each community state along a gradient of urbanization.
# pred_covs: A vector of the slope covariate (URB) used to predict the occupancy
#              rates for 'steady_CI'. Here, it is seq(-4,4,0.1).
# lab_x:     A Boolean that designates if the X axis of a subfigure should
#              be labeled with numbers. If TRUE, numbers are added. If FALSE,
#              only the axes are created.
# lab_y:     A Boolean that designates if the y axis of a subfigure should
#              be labeled with numbers. If TRUE, numbers are added. If FALSE,
#              only the axes are created.
# lab:       A character that designates what label should be put in the top
#              left corner of the sub-figure.
# d_list:    The 'data_list' object from the best fit model in
#              the 'fit_softmax_model.R' script.
# sp_init:   A matrix of observed community states at each site and season. This
#              is hard coded to be the 'sp_inits' object, which was created in
#              'fit_softmax_model.R' and supplied to the initial value functions
#              for JAGS. Look to 'fit_softmax_model.R' to determine how this
#              object is made.
# com:       A numeric element that indexes which community state is meant to 
#              be plotted. 1 = U, 2 = A, 3 = B, 4 = C, 5 = AB, 6 = AC, 7 = BC,
#              8 = ABC.

occ_plot <- function(steady = NULL, pred_covs = NULL, lab_x = NULL, lab_y = NULL,
  lab = NULL, d_list = data_list, sp_init = sp_inits, com = NULL){
  
  # make base plot
  plot(steady[,2,com] ~ pred_covs, type = 'n', ylim = c(0,1), bty = 'l', 
    xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  # these are for the 95% credible interval, used at the end
  x1 <- pred_covs
  x2 <- rev(pred_covs)
  y1 <- steady_CI[, 1, com]
  y2 <- rev(steady_CI[, 3, com])
  
  # function to add axis labels
  axes_occ_plot(x = lab_x, y = lab_y)
  # for sub-figures, textx and texty ensure that sub-figure names are
  # all in the same location.
  u <- par("usr")
  textx <- u[1] + 0.01 * (u[2] - u[1])
  texty <- u[4] - 0.07 * (u[4] - u[3])
  text(textx, texty, lab, pos = 4)
  
  # discretizing the urbanization covariates
  x <- d_list$gam_cov[, 2]
  cov_site_by_season <-
    matrix(x, ncol = ncol(sp_inits), nrow = length(x))
  # remove the one site with data > 3.5
  cov_site_by_season[cov_site_by_season > 3.5] <- 0
  is_com <- sp_init == com
  com_covs <- is_com * cov_site_by_season
  com_covs[com_covs == 0] <- NA
  not_com <- sp_init != com
  not_com_covs <- not_com * cov_site_by_season
  not_com_covs [not_com_covs == 0] <- NA
  tmp1 <- cut(com_covs, seq(-4, 4, 0.5))
  tmp2 <- cut(not_com_covs, seq(-4, 4, 0.5))
  
  myplot <- table(tmp1) / (table(tmp1) + table(tmp2))
  myplot[is.na(myplot)] <- NA
  #myplot[myplot == 0] <- NA
  myplot <- as.numeric(myplot)
  k <- table(tmp1)
  kn <-  table(tmp1) + table(tmp2)
  myses <- sqrt((myplot * (1 - myplot) * k) / kn)
  my_y <- seq(-3.75, 3.75, 0.5)
  points(myplot ~ my_y,
    pch = 16,
    cex = 0.8,
    col = "#606575")
  for (i in 1:length(my_y)) {
    if (is.na(myplot[i]))
      next
    y_se <- c(myplot[i] + myses[i], myplot[i] - myses[i])
    
    y_se[y_se > 1] <- 1
    y_se[y_se < 0] <- 0
    lines(y = y_se,
      x = rep(my_y[i], 2),
      col = "#606575")
  }
  
  polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .60), border = NA)
  lines(steady_CI[,2,com] ~ xp[,2], lwd = 3, col = alpha("#0D205F", 0.6))
  
}