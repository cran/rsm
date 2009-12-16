\name{bbd}
\alias{bbd}
\title{Generate a Box-Behnken design}
\description{
This function can generate a Box-Behnken design in 3 to 7 factors, and optionally will block it orthogonally if there are 4 or 5 factors.  It can also randomize the design and create a \code{coded.data} object.
}
\usage{
bbd(k, n0 = 4, block = (k == 4 | k == 5), randomize = TRUE, coding)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{A formula, or an integer giving the number of variables.  If the formula has a left-hand side, the variables named there are appended to the design and initialized to \code{NA}.}
  \item{n0}{Number of center points in each block.}
  \item{block}{Boolean specifying whether or not to block the design; or a character string (taken as \code{TRUE}) giving the desired name for the blocking factor.
  Only BBDs with 4 or 5 factors can be blocked.  A 4-factor BBD has three orthogonal blocks, and a 5-factor BBD has two.}
  \item{randomize}{Boolean determining whether or not to randomize the design.  If \code{block} is \code{TRUE}, each block is randomized separately.}
  \item{coding}{Optional list of formulas.  If this is provided, a \code{\link{coded.data}} object is returned with the specified coding.}
}
\details{
Box-Behnken designs (BBDs) are useful designs for fitting second-order response-surface models.  They use only three levels of each factor (compared with 5 for central-composite designs) and sometimes fewer runs are required than a CCD.
This function uses an internal table of BBDs; it only works for 3 to 7 factors.

If \code{k} is specified as a formula, the names in the formula determine the names of the factors in the generated design.  Otherwise, the names will be \code{x1, x2, ...}.
}
\value{
A \code{\link{data.frame}} with the generated design; or if \code{coding} is specified, a \code{\link{coded.data}} object.  The blocking variable, if present, will be a \code{\link{factor}}; all other variables will be numeric.
}
\note{
  To avoid aliasing the pure-quadratic terms, you must use a positive value \code{n0}.
}
\references{
Lenth RV (2009). ``Response-Surface Methods in R, Using rsm'', 
  \emph{Journal of Statistical Software}, 32(7), 1--17. 
  \url{http://www.jstatsoft.org/v32/i07/.}

Myers, RH, Montgomery, DC, and Anderson-Cook, CM (2009), 
\emph{Response Surface Methodology} (3rd ed.), Wiley.
}
\author{Russell V. Lenth}
\seealso{\code{\link{ccd}}, \code{\link{coded.data}}}
\examples{
library(rsm)

### Simple 3-factor case, not randomized so structure is evident
bbd(3, randomize=FALSE)

### 5-factor BBD, divided between two plants
bbd(y1 + y2 ~ A + B + C + D + E,  n0 = 5,  block = "Plant")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}
