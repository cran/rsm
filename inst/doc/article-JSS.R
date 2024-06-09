## ----echo = FALSE, results = "hide", message = FALSE--------------------------
require("rsm")
knitr::opts_chunk$set(fig.width = 4.5, class.output = "ro")

## -----------------------------------------------------------------------------
library(`rsm`)
ChemReact

## -----------------------------------------------------------------------------
CR1 <- coded.data(ChemReact1, x1 ~ (Time - 85)/5, x2 ~ (Temp - 175)/5)
CR1

## -----------------------------------------------------------------------------
as.data.frame(CR1)

## -----------------------------------------------------------------------------
code2val(data.frame(x1 = c(0.25, 0.5), x2 = c(-1.5, -0.5)), codings(CR1))

## -----------------------------------------------------------------------------
bbd(3, n0 = 2, coding =
  list(x1 ~ (Force - 20)/3, x2 ~ (Rate - 50)/10, x3 ~ Polish - 4))

## -----------------------------------------------------------------------------
ccd.pick(5, n.c = c(8, 16), blks.c = c(1, 2, 4),
  wbr.s = 1:2, restrict = "N<=65")

## -----------------------------------------------------------------------------
des1 <- ccd (y1 + y2 ~ A + B + C + D,
  generators = E ~ - A * B * C * D, n0 = c(6, 1))

## -----------------------------------------------------------------------------
des10 <- ccd( ~ A + B + C + D + E,
  blocks = Blk ~ c(A * B * C, C * D * E), n0 = c(2, 4))

## ----fig=TRUE, fig.height=5.5, fig.width=5------------------------------------
varfcn(des10, ~ Blk + SO(A,B,C,D,E), dist = seq(0, 3, by = .1))
varfcn(des10, ~ Blk + SO(A,B,C,D,E), dist = seq(0, 3, by = .1), contour = TRUE)

## -----------------------------------------------------------------------------
ccd(2, n0 = c(1, 1), inscribed = TRUE, randomize = FALSE)

## -----------------------------------------------------------------------------
CR1.rsm <- rsm(Yield ~ FO(x1, x2), data = CR1)
summary(CR1.rsm)

## ----results="hide"-----------------------------------------------------------
CR1.rsmi <- update(CR1.rsm, . ~ . + TWI(x1, x2))
summary(CR1.rsmi)

## -----------------------------------------------------------------------------
( CR2 <- djoin(CR1, ChemReact2) )

## -----------------------------------------------------------------------------
CR2.rsm <- rsm(Yield ~ Block + SO(x1, x2), data = CR2)
summary(CR2.rsm)

## -----------------------------------------------------------------------------
heli.rsm <- rsm(ave ~ block + SO(x1, x2, x3, x4), data = heli)
summary(heli.rsm)

## ----fig=TRUE, fig.height=6, fig.width=8--------------------------------------
par(mfrow = c(2, 3))
contour(heli.rsm, ~ x1 + x2 + x3 + x4, image = TRUE,
  at = summary(heli.rsm)$canonical$xs)

## -----------------------------------------------------------------------------
steepest(CR1.rsm, dist = c(0, 0.5, 1))

## -----------------------------------------------------------------------------
canonical.path(heli.rsm, dist = seq(-5, 5, by = 0.5))

## -----------------------------------------------------------------------------
CO = as.coded.data(codata,  x1 ~ (Ethanol - 0.2)/0.1,  x2 ~ A.F.ratio - 15)
names(CO)[3] = "CO.conc"
head(CO)

## -----------------------------------------------------------------------------
CO.rsm = rsm(CO.conc ~ SO(x1,x2), data = CO)
canonical(CO.rsm)

## -----------------------------------------------------------------------------
canonical(CO.rsm, threshold = 0)

## ----fig=TRUE, fig.height=5, message=FALSE------------------------------------
contour(CO.rsm, x2 ~ x1, bounds = list(x1 = c(-16, 2), x2 = c(-2, 16)), 
        zlim = c(-100, 100), col = "gray", decode = FALSE)
lines(c(-1,1,1,-1,-1), c(-1,-1,1,1,-1), col = "green") # design region
points(x2 ~ x1, data = canonical.path(CO.rsm), 
        col = "blue", pch = 1 + 6*(dist == 0))
points(x2 ~ x1, data = canonical.path(CO.rsm, threshold = 0), 
        col = "red", pch = 1 + 6*(dist == 0))
points(x2 ~ x1, data=steepest(CO.rsm), 
        col = "magenta", pch = 1 + 6*(dist == 0))

