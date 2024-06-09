## ----echo = FALSE, results = "hide", message = FALSE--------------------------------------------------------
require("rsm")
knitr::opts_chunk$set(fig.width = 5.5, class.output = "ro")

## -----------------------------------------------------------------------------------------------------------
swiss2.lm <- lm(Fertility ~ poly(Agriculture, Education, degree=2), data=swiss)

## ----fig=TRUE, fig.width=8.5, fig.height=3------------------------------------------------------------------
library(rsm)
par(mfrow=c(1,3))
image(swiss2.lm, Education ~ Agriculture)
contour(swiss2.lm, Education ~ Agriculture)
persp(swiss2.lm, Education ~ Agriculture, zlab = "Fertility")

## ----fig=TRUE, fig.height=5.5-------------------------------------------------------------------------------
persp(swiss2.lm, Education ~ Agriculture, col = "blue", 
  bounds = list(Agriculture=c(20,70), Education=c(0,30)),
  zlab = "Predicted Fertility", 
  contours = list(z="top", col="orange", shade = 1), 
  theta = -135, phi = 35)

## -----------------------------------------------------------------------------------------------------------
heli.rsm <- rsm(ave ~ block + SO(x1,x2,x3,x4), data = heli)

## ----fig=TRUE, fig.height=6, fig.width=8--------------------------------------------------------------------
par(mfrow = c(2,3))
contour (heli.rsm, ~ x1 + x2 + x3 + x4)

## -----------------------------------------------------------------------------------------------------------
xs <- canonical(heli.rsm)$xs          # stat.pt. in coded units
SP <- code2val(xs, codings(heli.rsm)) # in decoded units
myhook <- list()
myhook$post.plot <- function(lab) {
  idx <- sapply(lab[3:4], grep, names(xs))
  points (SP[idx[1]], SP[idx[2]], pch = 2, col = "red")
}

## ----fig=TRUE, fig.height=6, fig.width=8--------------------------------------------------------------------
par(mfrow = c(2,3))
contour (heli.rsm, ~ x1 + x2 + x3 + x4, image = TRUE,
  at = xs, hook = myhook)

## ----fig=TRUE, fig.height = 5.5-----------------------------------------------------------------------------
persp (heli.rsm, x2 ~ x1, at = xs, col = rainbow(50), contours = "colors")
persp (heli.rsm, x4 ~ x1, at = xs, col = rainbow(50), contours = "colors")

