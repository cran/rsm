### R code from vignette source 'rsm-plots.Rtex'

###################################################
### code chunk number 1: rsm-plots.Rtex:62-63
###################################################
swiss2.lm <- lm(Fertility ~ poly(Agriculture, Education, degree=2), data=swiss)


###################################################
### code chunk number 2: basics
###################################################
library(rsm)
par(mfrow=c(1,3))
image(swiss2.lm, Education ~ Agriculture)
contour(swiss2.lm, Education ~ Agriculture)
persp(swiss2.lm, Education ~ Agriculture, zlab = "Fertility")


###################################################
### code chunk number 3: lessbas
###################################################
persp(swiss2.lm, Education ~ Agriculture, col = "blue", 
  bounds = list(Agriculture=c(20,70), Education=c(0,30)),
  zlab = "Predicted Fertility", contours = list(z="top", col="orange"), 
  theta = -145, phi = 35, shade = 1)


###################################################
### code chunk number 4: rsm-plots.Rtex:100-101
###################################################
heli.rsm <- rsm(ave ~ block + SO(x1,x2,x3,x4), data = heli)


###################################################
### code chunk number 5: cont1
###################################################
par(mfrow = c(2,3))
contour (heli.rsm, ~ x1 + x2 + x3 + x4)


###################################################
### code chunk number 6: rsm-plots.Rtex:133-139
###################################################
xs <- canonical(heli.rsm)$xs
myhook <- list()
myhook$post.plot <- function(lab) {
  idx <- sapply(lab[3:4], grep, names(xs))
  points (xs[idx[1]], xs[idx[2]], pch=2, col="red")
}


###################################################
### code chunk number 7: cont2
###################################################
par(mfrow = c(2,3))
contour (heli.rsm, ~ x1 + x2 + x3 + x4, image = TRUE,
  at = xs, hook = myhook)


###################################################
### code chunk number 8: rsm-plots.Rtex:154-157 (eval = FALSE)
###################################################
## pdf(file = "heli-cps.pdf")
## contour (heli.rsm, ~ x1 + x2 + x3 + x4, image = TRUE, at = xs, hook = myhook)
## dev.off()


###################################################
### code chunk number 9: rsm-plots.Rtex:167-173 (eval = FALSE)
###################################################
## png.hook <- list()
## png.hook$pre.plot <- function(lab) 
##     png(file = paste(lab[3], lab[4], ".png", sep = ""))
## png.hook$post.plot = function(lab)
##     dev.off()
## contour (heli.rsm, ~ x1 + x2 + x3 + x4, image = TRUE, at = xs, hook = png.hook)


###################################################
### code chunk number 10: persp
###################################################
persp (heli.rsm, ~ x1 + x2 + x3 + x4, at = xs, 
  col = rainbow(50), contours = "colors")


