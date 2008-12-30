\name{contour.lm}
\alias{contour.lm}
\title{Contour plot(s) of a fitted linear model}
\description{
  \code{contour} method that will work with any \code{lm} object involving two or
  more numerical predictors.
}
\usage{
\method{contour}{lm}(x, form, at, bounds, zlim, image = TRUE, img.col = terrain.colors(50), ...)
}
\arguments{
  \item{x}{A \code{lm} object.}
  \item{form}{A formula, or a list of formulas.}
  \item{at}{Optional named list of fixed values to use for variables other
    than the axis variables.  If not provided, the mean is used for numeric variables,
    or the first level if it is a factor.}
  \item{bounds}{Optional \emph{named} list of bounds or grid values to use for
    the variables having the same names.  See details.}
  \item{zlim}{Optional \code{zlim} setting passed to \code{\link{image}}.
    If not provided, the range of values across all plotted surfaces is used.}
  \item{image}{Set to \code{TRUE} if you want an image plot underneath the contours.}
  \item{img.col}{Color map to use when \code{image=TRUE}.}
  \item{\dots}{Additional arguments passed to \code{\link{contour}}.}
}
\details{
  \code{form} may be a single formula or a list of formulas.  A simple formula like
  \code{x2 ~ x1} will produce a contour plot of the fitted regression surface
  for combinations of \code{x2} (vertical axis) and \code{x1} (horizontal axis).
  A list of several such simple formulas will produce a contour plot for each formula.
  A two-sided formula produces contour plots for each left-hand variable versus each
  right-hand variable (except when they are the same); for example, 
  \code{x1+x3 ~ x2+x3} is equivalent to 
  \code{list(x1~x2, x3~x2, x1~x3)}.  
  A one-sided formula produces contour plots for each pair of variables.  For example,
  \code{~ x1+x2+x3} is equivalent to 
  \code{list(x2~x1, x3~x1, x3~x2)}. 
  
  For any variables not in the \code{bounds} argument, a grid of 26 equally-spaced 
  values in the observed range of that variable is used.  If you specify a vector of 
  length 2, it is interpreted as the desired range for that variable and a grid of 26
  equally-spaced points is generated.  If it is a vector of length 3, the first two elements are used
  as the range, and the third as the number of grid points.
  If it is a vector of length 4 or more, those
  values are used directly as the grid values.
}
\author{Russell V. Lenth}
\seealso{\code{\link{contour}}}
\examples{
library (rsm)
heli.rsm = rsm (ave ~ block + SO(x1, x2, x3, x4), data = heli)

par (mfrow = c (2,3))
contour (heli.rsm, ~x1+x2+x3+x4, at=summary(heli.rsm)$canonical$xs)
}
\keyword{ regression }
\keyword{ hplot }
