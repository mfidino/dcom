###########################
#
# plotting script
#
# Written by Mason Fidino
#
###########################

# load packages
library(runjags)
library(coda)
library(scales)

# need make_inxs function
source("fit_models_utility_functions.R")
source("plotting_script_utility_functions.R")

# read in the best fit model
mout <- readRDS("./model_output/model_3_inxs_intcov.rds")

# convert it to a matrix
mm <- as.matrix(as.mcmc.list(mout), chains = TRUE)

# remove the 'z' parameters
mm <- mm[, -grep("z", colnames(mm))]

# create covariates for predictions
# the urbanization pca went from roughly -4 to 4
xp <- seq(-4, 4, by = 0.1)
xp <- cbind(1, xp) # add intercept term


# create the colonization and extinction figures. I've left the code
#  seperate for each subfigure as there are a number of unique actions
#  that need to be taken for each specific sub-figure. Following this,
#  vector graphics of the species were added via Inkscape and the image
#  was resized via ImageMagick



png("./plots/colonization.png", height = 12, width = 12, units = "in",
    res = 300)

# set layout manually
m <- matrix(10, ncol = 14, nrow = 14)
m[2:5,2:13] <- rep(1:3, each = 16)
m[6:9, 2:13] <- rep(4:6, each = 16)
m[10:13, 2:13] <- rep(7:9, each = 16)
layout(m)
par(mar = c(2,2,2,2), xpd = NA)
# first
plot(1~1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(-4,4), ylim = c(0,1))


mtext("Colonization of\ncoyote...", cex = 1.7)
# opossum|coyote
tp <- make_preds(mm, sp = c(2,1), "b", "g", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", las = 1,
     xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
mtext("Colonization of\nopossum...", cex = 1.7)

# raccoon|Coyote
tp <- make_preds(mm, sp = c(3,1), "b", "g", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1",
     xaxt= "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
mtext("Colonization of\nraccoon...", cex = 1.7)
text("conditional on\ncoyote", x = 5.5, y = 0.5, las = 1, srt = 270,
     cex = 3)

# coyote|opossum
tp <- make_preds(mm, sp = c(1,2), "b", "g", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
text("Colonization probability", x = -6.7, y = 0.5, las = 1, srt = 90,
     cex = 4)

# middle
plot(1~1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(-4, 4), ylim = c(0,1))
par(fg = "white", lend = 2)
legend(x = -4, y = 0.9, c("Absent", "Present"), lwd = 20, 
       col = c(alpha("#3854B1", 0.3), alpha("#2CC44E", 0.3)), lty = c(1,1), cex = 1.9, seg.len = 1, bty = "n" , horiz = TRUE, title = "Conditional species is..." )
par(fg = "black", lend = 0)
legend(x = -4, y = 0.9, c("Absent", "Present"), lwd = 4, 
       col = c("#3854B1", "#2CC44E"), lty = c(1,3), cex = 1.9, bg = "NA", seg.len = 1,
       bty = "n", horiz = TRUE, title = "Conditional species is...")

# raccoon|opossum
tp <- make_preds(mm, sp = c(3,2), "b", "g", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
text("conditional on\nopossum", x = 5.5, y = 0.5, las = 1, srt = 270,
     cex = 3)

# coyote|raccoon
tp <- make_preds(mm, sp = c(1,3), "b", "g", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()

# opossum|raccoon
tp <- make_preds(mm, sp = c(2,3), "b", "g", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
text("Urbanization", x = 0, y = -0.27, las = 0,
     cex = 4)
arrows(x0 = 3.8, y0 = -0.35, x1 = 14, lwd = 2, length = 0.15)
arrows(x0 = -3.8, y0 = -0.35, x1 = -14, lwd = 2, length = 0.15)

text("Higher housing density", x = 8.7, y = -0.41, cex = 2)
text("Higher tree cover", x = -8.9, y = -0.41, cex = 2)

# end
plot(1~1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(-4,4), ylim = c(0,1))
#legend("center", c("Absent", "Present"), lwd = 8, col = c("#3854B1", "#2CC44E"),
#       bty = "n")
text("conditional on\nraccoon", x = 5.5, y = 0.5, las = 1, srt = 270,
     cex = 3)


dev.off()


png("./plots/extinction.png", height = 12, width = 12, units = "in",
    res = 300)

m <- matrix(10, ncol = 14, nrow = 14)
m[2:5,2:13] <- rep(1:3, each = 16)
m[6:9, 2:13] <- rep(4:6, each = 16)
m[10:13, 2:13] <- rep(7:9, each = 16)
layout(m)
par(mar = c(2,2,2,2), xpd = NA)
# first
plot(1~1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "")

mtext("Extinction of\ncoyote...", cex = 1.7)
# opossum|coyote
tp <- make_preds(mm, sp = c(2,1), "d", "h", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", las = 1,
     xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
mtext("Extinction of\nopossum...", cex = 1.7)

# raccoon|Coyote
tp <- make_preds(mm, sp = c(3,1), "d", "h", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1",
     xaxt= "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
mtext("Extinction of\nraccoon...", cex = 1.7)
text("conditional on\ncoyote", x = 5.5, y = 0.5, las = 1, srt = 270,
     cex = 3)

# coyote|opossum
tp <- make_preds(mm, sp = c(1,2), "d", "h", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
text("Extinction probability", x = -6.7, y = 0.5, las = 1, srt = 90,
     cex = 4)

# middle
plot(1~1, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(-4,4), ylim = c(0,1))
par(fg = "white", lend = 2)
legend(x = -4, y = 0.9, c("Absent", "Present"), lwd = 20, 
       col = c(alpha("#3854B1", 0.3), alpha("#2CC44E", 0.3)), lty = c(1,1), cex = 1.9, seg.len = 1, bty = "n" , horiz = TRUE, title = "Conditional species is..." )
par(fg = "black", lend = 0)
legend(x = -4, y = 0.9, c("Absent", "Present"), lwd = 4, 
       col = c("#3854B1", "#2CC44E"), lty = c(1,3), cex = 1.9, bg = "NA", seg.len = 1,
       bty = "n", horiz = TRUE, title = "Conditional species is...")

# raccoon|opossum
tp <- make_preds(mm, sp = c(3,2), "d", "h", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
text("conditional on\nopossum", x = 5.5, y = 0.5, las = 1, srt = 270,
     cex = 3)

# coyote|raccoon
tp <- make_preds(mm, sp = c(1,3), "d", "h", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()

# opossum|raccoon
tp <- make_preds(mm, sp = c(2,3), "d", "h", xp)
x1 <- xp[,2]
x2 <- rev(xp[,2])
y1 <- tp$wo[1,]
y2 <- rev(tp$wo[3,])
plot(tp$wo[2,] ~ xp[,2], lwd = 2, type = "l", ylim = c(0,1), bty = "l",
     xlab = "", ylab = "",
     main = "", col = "#3854B1", xaxt = "n", yaxt = "n")
polygon(c(x1, x2), c(y1, y2), col = alpha("#3854B1", .30), border = NA)
polygon(c(x1, x2), c(tp$w[1,], rev(tp$w[3,])), 
        col = alpha("#2CC44E", 0.3), border = NA)
lines(tp$wo[2,] ~ xp[,2], lwd = 3, col = "#3854B1")
lines(tp$w[2, ] ~ xp[,2], lwd = 3, lty = 3, col = "#2CC44E")
make_axes_ce()
text("Urbanization", x = 0, y = -0.27, las = 0,
     cex = 4)
arrows(x0 = 3.8, y0 = -0.35, x1 = 14, lwd = 2, length = 0.15)
arrows(x0 = -3.8, y0 = -0.35, x1 = -14, lwd = 2, length = 0.15)
#text("4", x = 14.7, y = -0.35, cex = 3)
#text("-4", x = -14.7, y = -0.35, cex = 3)
text("Higher housing density", x = 8.7, y = -0.41, cex = 2)
text("Higher tree cover", x = -8.9, y = -0.41, cex = 2)

# end
plot(1~1, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(-4,4), ylim = c(0,1))
text("conditional on\nraccoon", x = 5.5, y = 0.5, las = 1, srt = 270,
     cex = 3)

dev.off()


# plot out the steady state as well (figure 3). see calculate_steady_state.R
# for details.

smat <- calculate_steady_state(mm = mm, data_list = data_list, ncores = 7)

# put into array
steady_mat <- array(0, dim = c(nrow(mm), nrow(xp), 8))
for(i in 1:nrow(mm)){
  steady_mat[i,,] <- smat[[i]]
}

# this is the median and 95% CI across the urban PCA for each community state.
# covariate x low, med, high x community state
steady_CI <- array(0, dim = c(81, 3, 8))
for (i in 1:81) {
  steady_CI[i, , ] <- apply(steady_mat[, i, ], 2,
    quantile, probs = c(0.025, 0.5, 0.975))
}

#########################################
# plot figure 3
#########################################




png("./plots/fidino_3.png", height = 6, width = 4, units = "in", res = 300)

# create the layout for the plot
m <- matrix(9, ncol = 15, nrow = 14)
m[2:4, c(2:7)  +1 ] <- 1
m[2:4, c(8:13) +1 ] <- 2
m[5:7, c(2:7)  +1 ] <- 3
m[5:7, c(8:13) +1 ] <- 4
m[8:10, c(2:7)  +1  ] <- 5
m[8:10,  c(8:13) +1 ] <- 6
m[11:13,c(2:7)  +1 ] <- 7
m[11:13, c(8:13) +1 ] <- 8

# set layout for plot
layout(m)

# set margins for better looking plot and xpd = NA so we can add
#  text to the margin of the plots
par(mar = c(1.5,0.8,0.8,1.5), xpd = NA)

# first
occ_plot(steady_CI, pred_covs = xp[,2], lab_x = FALSE, lab_y = TRUE, com = 1,
  lab = "A) No species")
# second
occ_plot(steady_CI, pred_covs = xp[,2], lab_x = FALSE, lab_y = FALSE, com = 2,
  lab = "B) Coyote")
# third
occ_plot(steady_CI, pred_covs = xp[,2], lab_x = FALSE, lab_y = TRUE, com = 3,
  lab = "C) Opossum")
# fourth
occ_plot(steady_CI, pred_covs = xp[,2], lab_x = FALSE, lab_y = FALSE, com = 4,
  lab = "D) Raccoon")
# fifth
occ_plot(steady_CI, pred_covs = xp[,2], lab_x = FALSE, lab_y = TRUE, com = 5,
  lab = "E) Coyote & opossum")
# sixth
occ_plot(steady_CI, pred_covs = xp[,2], lab_x = FALSE, lab_y = FALSE, com = 6,
  lab = "F) Coyote & raccoon")
# seventh
occ_plot(steady_CI, pred_covs = xp[,2], lab_x = TRUE, lab_y = TRUE, com = 7,
  lab = "G) Opossum & raccoon")
# eigth
occ_plot(steady_CI, pred_covs = xp[,2], lab_x = TRUE, lab_y = FALSE, com = 8,
  lab = "H) No species")

text(-5.25, -0.4, "Urbanization", cex = 1.5, xpd = NA)
text(-18, 2.6, srt = 90, "Occupancy probability", cex = 1.5)

# the following code was run via imagemagick to generate the final figure

convert fidino_3.png -trim -bordercolor white -border 10 -resize 900x1350 fidino_3.eps

