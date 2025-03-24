## ----echo = FALSE, results = "hide", message = FALSE--------------------------
require("rsm")
knitr::opts_chunk$set(fig.width = 4.5, class.output = "ro")

suppressWarnings(RNGversion("2.15.3")) ### This saved my bacon reproducing a very old vignette!
set.seed(19482012)

# basic generator
simb = function(flour, sugar, butter, sd=0.58) {
  x1 = (flour - 1.2)/.1
  x2 = (sugar - .25)/.1
  x3 = (butter - .4)/.15
  y = 32 - x1^2 - .5*x3^2 - sqrt(abs(x2 + x3))^2
  rnorm(length(y), y, sd)
}

# requires DECODED data!
simBake = function(data) {
  blkeff = rnorm(1, 0, 2.3)
  round(blkeff + simb(data$flour, data$sugar, data$butter), 1)
}

## -----------------------------------------------------------------------------
library(rsm)
expt1 = cube(~ x1 + x2,  x3 ~ x1 * x2, n0 = 4,
            coding = c(x1 ~ (flour - 1)/.1, x2 ~ (sugar - .5)/.1, x3 ~ (butter - .25)/.1))

## -----------------------------------------------------------------------------
expt1

## -----------------------------------------------------------------------------
as.data.frame(expt1)

## ----echo = FALSE-------------------------------------------------------------
SAVESEED = .Random.seed

## ----fig=TRUE, fig.width=6, fig.height=3.5------------------------------------
par(mfrow=c(1,2))
varfcn(expt1, ~ FO(x1,x2,x3))
varfcn(expt1, ~ FO(x1,x2,x3), contour = TRUE)

## ----error = TRUE-------------------------------------------------------------
try({
varfcn(expt1, ~ SO(x1,x2,x3))
})

## -----------------------------------------------------------------------------
try(
  djoin(expt1, star(n0 = 2, alpha = "rotatable"))
)

## ----fig=TRUE, fig.height=3.5, fig.width=6------------------------------------
par(mfrow=c(1,2))
followup = djoin(expt1, star(n0 = 2, alpha = 1.5))
varfcn(followup, ~ Block + SO(x1,x2,x3), main = "Followup")
varfcn(followup, ~ Block + SO(x1,x2,x3), contour = TRUE, main = "Block + SO(x1,x2,x3)")

## ----echo = FALSE-------------------------------------------------------------
.Random.seed = SAVESEED
expt1$rating = simBake(decode.data(expt1))

## -----------------------------------------------------------------------------
expt1

## -----------------------------------------------------------------------------
anal1 = rsm(rating ~ FO(x1,x2,x3), data=expt1)
summary(anal1)

## -----------------------------------------------------------------------------
( sa1 = steepest(anal1) )

## -----------------------------------------------------------------------------
expt2 = dupe(sa1[2:9, ])

## ----echo = FALSE-------------------------------------------------------------------------------------------
expt2$rating = simBake(expt2)
options(width=110)

## -----------------------------------------------------------------------------------------------------------
expt2

## ----fig = TRUE, scale=.46, fig.height=4.5------------------------------------------------------------------
plot(rating ~ dist, data = expt2)
anal2 = lm(rating ~ poly(dist, 2),  data = expt2)
with(expt2, {
    ord = order(dist)
    lines(dist[ord], predict(anal2)[ord])
})

## -----------------------------------------------------------------------------------------------------------
expt3 = dupe(expt1)
codings(expt3) = c(x1 ~ (flour - 1.25)/.1,  x2 ~ (sugar - .45)/.1,  x3 ~ (butter - .25)/.1)

## ----echo = FALSE-------------------------------------------------------------------------------------------
expt3$rating = simBake(decode.data(expt3))

## -----------------------------------------------------------------------------------------------------------
expt3

## -----------------------------------------------------------------------------------------------------------
anal3 = rsm(rating ~ FO(x1,x2,x3), data=expt3)
summary(anal3)

## -----------------------------------------------------------------------------------------------------------
expt4 = foldover(expt3, variable = "x1")
expt4$rating = NULL  ### discard previous rating data
expt4   # Here's the new protocol

## ----echo = FALSE-------------------------------------------------------------------------------------------
expt4$rating = simBake(decode.data(expt4))

## -----------------------------------------------------------------------------------------------------------
expt4

## -----------------------------------------------------------------------------------------------------------
names( djoin(expt3, expt4) )

## -----------------------------------------------------------------------------------------------------------
anal4 = rsm(rating ~ Block + FO(x1,x2,x3), data = djoin(expt3, expt4))
summary(anal4)

## ----fig = TRUE, fig.height=3.5, fig.width=6----------------------------------------------------------------
expt5 = star(expt4, n0 = 2, alpha = "orthogonal")
par(mfrow=c(1,2))
comb = djoin(expt3, expt4, expt5)
varfcn(comb, ~ Block + SO(x1,x2,x3), main = "Further augmented")
varfcn(comb, ~ Block + SO(x1,x2,x3), contour = TRUE, main = "2nd order")

## ----echo = FALSE-------------------------------------------------------------------------------------------
expt5$rating = simBake(decode.data(expt5))

## -----------------------------------------------------------------------------------------------------------
expt5

## -----------------------------------------------------------------------------------------------------------
anal5 = rsm(rating ~ Block + SO(x1,x2,x3), data = djoin(expt3, expt4, expt5))
summary(anal5)

## -----------------------------------------------------------------------------------------------------------
steepest(anal5)

## -----------------------------------------------------------------------------------------------------------
expt6 = dupe(steepest(anal5, dist = (2:9)/3))

## ----echo = FALSE-------------------------------------------------------------------------------------------
expt6$rating = simBake(expt6)

## -----------------------------------------------------------------------------------------------------------
expt6

## ----fig = TRUE, fig.height=3.5-----------------------------------------------------------------------------
par(mar=c(4,4,0,0)+.1)
plot(rating ~ dist, data = expt6)
anal6 = lm(rating ~ poly(dist, 2),  data = expt6)
with(expt6, {
    ord = order(dist)
    lines(dist[ord], predict(anal6)[ord])
})

## -----------------------------------------------------------------------------------------------------------
expt7  =  ccd( ~ x1 + x2 + x3,  n0 = c(0, 2),  alpha = "orth",  coding  =  c(
            x1 ~ (flour - 1.25)/.1,  x2 ~ (sugar - .3)/.1,  x3 ~ (butter - .3)/.1))

## ----echo = FALSE-------------------------------------------------------------------------------------------
expt7$rating = simBake(decode.data(expt7))

## -----------------------------------------------------------------------------------------------------------
expt7

## -----------------------------------------------------------------------------------------------------------
anal7 = rsm(rating ~ Block + SO(x1,x2,x3), data = expt7)
summary(anal7)

## ----fig = TRUE, fig.width=8, fig.height=2.5----------------------------------------------------------------
par(cex.lab=1.25, cex.axis=1, cex.sub=1.5, mar=.1+c(4.5,7,0,0))
par(mfrow=c(1,3))
contour(anal7, ~ x1 + x2 + x3, at = xs(anal7), image = TRUE)

## ----fig=TRUE, fig.width=7, fig.height=2.3------------------------------------------------------------------
fits = predict(anal7)
resids = resid(anal7)
boot.raw = suppressMessages(
    replicate(200, xs(update(anal7, fits + sample(resids, replace=TRUE) ~ .))))
boot = code2val(as.data.frame(t(boot.raw)), codings=codings(anal7))
par(mar=.1+c(4,5,0,0), cex.lab=1.5)
par(mfrow = c(1,3))
plot(sugar ~ flour, data = boot, col = "gray");   points(1.215, .282, col = "red", pch = 7)
plot(butter ~ flour, data = boot, col = "gray");  points(1.215, .364, col = "red", pch = 7)
plot(butter ~ sugar, data = boot, col = "gray");  points(.282, .364, col = "red", pch = 7)

