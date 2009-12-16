\name{ccd}
\alias{ccd}
\title{Generate a central-composite design}
\description{
This function can generate a broad range of central-composite designs with flexible choices of replications, aliasing of predictors and fractional blocks, and choices of axis or 'star' points.
}
\usage{
ccd(basis, generators, blocks = "Block", n0 = 4, alpha = "orthogonal", 
    wbreps = 1, bbreps = 1, randomize = TRUE, coding)
}
\arguments{
  \item{basis}{A formula, or an integer giving the number of variables.  If the formula has a left-hand side, the variables named there are appended to the design and initialized to \code{NA}.}
  \item{generators}{Optional formula or list of formulas to generate aliased variables}
  \item{blocks}{A string or a formula.  If a character string, it is the name of the blcoking factor; if a formula, the left-hand side is used as the name of the blcoking factor, and the formula(s) on the right-hand side are used to generate separate fractional blocks.}
  \item{n0}{Integer or vector of two integers giving the number of center points.
    If only one value, this is the number of center pouints in each block.
    If two values are given, \code{n0[1]} specifies the number of center points in 
    'cube' blocks, and \code{n0[2]} specifies the number of center points in 'star' blocks.}
  \item{alpha}{If numeric, the position of the 'star' points.  May also be a character string; if it matches \code{"orthogonal"},the star points are positioned to block the design orthogonally; if it matches \code{"rotatable"}, the star points are chosen to make the design rotatable.  The default is to generate an orthogonally blocked design.}
  \item{wbreps}{Number(s) of within-block replications.  If this is a vector of length 2, then separate numbers are used for the 'cube' and the 'axis' blocks respectively.}
  \item{bbreps}{Number(s) of between-block replications (i.e., number of repeats of each block).  f this is a vector of length 2, then separate numbers are used for the 'cube' and the 'axis' blocks respectively.}
  \item{randomize}{Boolean determinig whether or not to randomize the design.  Each block is randomized separately.}
  \item{coding}{Optional list of formulas.  If this is provided, a \code{\link{coded.data}} object is returned with the specified coding.}
}
\details{
Central-composite designs (CCDs) are popular designs for use in response-surface exploration.  They are blocked designs consisting of at least one 'cube' block (two-level factorial or fractional factorial, plus center points) , and at least one 'star' block (points along each axis at positions \code{-alpha} and \code{+alpha}), plus center points.  Everything is assumed to be on a coded scale, where the cube portion of the design has values of -1 and 1 for each variable, and the center points are 0.  The codings may be provided, if desired, using the \code{coding} argument.

The \code{basis} argument determines a basic design used to create cube blocks.  For example, \code{basis = ~A+B+C} would generate a basic design of 8 runs.     Use \code{generators} if you want additional variables; for example, \code{generators = c(D~-A*B, E~B*C)} added to the above would generate a 5-factor design with defining relation \code{I = -ABD = BCE = -ACDE}.  

If you want the cube points divided into fractional blocks, give the formula(s) in the \code{blocks} argument.  For instance, suppose we have \code{basis = A+B+C+D+E}, \code{generators = F~-A*C*D}, and \code{blocks = Day ~ c(A*B*C, C*D*E)}.  Then the blocking variable will be named \code{"Day"}, and the basic 32-run design (half-fraction of 6 factors) will be divided into 4 blocks of 8 runs each, based on the combination of signs of \code{A*B*C} and \code{C*D*E}.  Notice that \code{Day} will be confounded with its generators, the interaction thereof, and all aliases of these: \code{Day = ABC = CDE = ABDE = -BDF = -ADCF = -BCEF = -AEF}.  To each of these blocks, we will add \code{n0} (or \code{n0[1]}) center points.

A star block contains \code{n0} (or \code{n0[2]}) center points, plus points at \code{+/- alpha} on each coordinate axis.  You may specify the \code{alpha} you want, or have it chosen to achieve orthogonality of blocks (the default) or rotatability.  Conditions for these criteria to hold are described in detail in references such as Myers and Montgomery (2005).

\code{wbreps} values greater than 1 cause cube points or star points to be replicated in the same block.  \code{bbreps} values greater than 1 cause additional blocks to be added to the design.  By default, the design is randomized so that the run order in each block is random.  The order of appearance of the blocks is not randomized.

A couple of convenience features are added.  If \code{basis} is an integer instead of a formula, default variable names of \code{x1, x2, ...} are used; for example, \code{basis = 3} is equivalent to \code{basis = ~x1+x2+x3}.  You may specify generators or fractional blocks using the same notation.  You may also give variables on the left-hand side when \code{basis} is a formula, and these variables are added to the returned design.  For example, \code{Yield+ProcessTime ~ A+B+C} would generate a design in factors \code{A, B, C} (plus others if there are generators), and additional variables \code{Yield} and \code{ProcessTime} initialized at \code{NA}.
}
\value{
A \code{\link{data.frame}} with the generated design; or if \code{coding} is specified, a \code{\link{coded.data}} object.  The blocking variable will be a \code{\link{factor}}; all other variables will be numeric.
}
\references{
Lenth RV (2009). ``Response-Surface Methods in R, Using rsm'', 
  \emph{Journal of Statistical Software}, 32(7), 1--17. 
  \url{http://www.jstatsoft.org/v32/i07/.}

Myers, RH, Montgomery, DC, and Anderson-Cook, CM (2009), 
\emph{Response Surface Methodology} (3rd ed.), Wiley.
}
\author{Russell V. Lenth}
\note{
  Poor choices of \code{generators} and/or \code{blocks} can alias or partially alias some effects needed to estimate a second-order response surface.  \code{ccd} runs a trial model on just the cube portion of the experiment, and issues a warning if there are aliased effects.  The presence of center points and axis points may help mitigate this problem and make the surface estimable, but it is best to avoid designs where this aliasing occurs.
  
  The function \code{\link{ccd.pick}} is available to help determine good choices for arguments to this function.
  
  In  practice, you may generate the whole design, but initially collect data only on one or more 'cube' blocks, which is good enough for estimating a first-order model.  Depending on the results of analysis, we either may follow a path of steepest ascent, or continue with data from the 'star' block(s) (and possibly more 'cube' blocks) so that we can estimate a second-order model.
  
  An alternative to a CCD when you want to go straight to second-order modeling is a Box-Behnken design, generated by \code{\link{bbd}}.  These designs are not as various or flexible as CCDs, but they can require fewer runs.
}
\seealso{\code{\link{ccd.pick}}, \code{\link{coded.data}}, \code{\link{bbd}}}
\examples{
library(rsm)

### Generate a standard 3-variable CCD with 12 runs in each block
des = ccd(3, n0=c(4,6), coding = list(x1 ~ (Temp - 150)/10, 
          x2 ~ (Pres - 50)/5, x3 ~ Feedrate - 4))
decode.data(des)

### Generate a 5-variable design in 2 blocks.  The cube block has 16 runs
### This design will have alpha = 2; it is both orthogonal and rotatable
ccd(~ x1 + x2 + x3 + x4,  x5 ~ x1 * x2 * x3 * x4, n0 = c(6,1))

### Generate a 5-variable design with 5 blocks:
### 4 blocks of 8 cube points each, and 1 block with star points 
### You'll get alpha = 2.366; if you add alpha = "rot", you'll get 2.378
des = ccd(y1 + y2 ~ A + B + C + D + E, , Shift ~ c(-A*B*C, C*D*E), n0=c(2,4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}
