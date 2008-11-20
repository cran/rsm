\name{contour.lm}
\alias{contour.lm}
\title{Contour plot of a linear model}
\description{
  \code{contour} method that will work with any \code{lm} object involving two or
  more numerical predictors.
}
\usage{
\method{contour}{lm}(x, varlist, image = TRUE, img.col = terrain.colors(50), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{lm} object.}
  \item{varlist}{A named list with at least two elements.  
  The names should match those of predictors in the model.
  The first element defines the horizontal variable, and the second
  one defines the vertical variable.  Subsequent elements, if specified,
  specify values to hold fixed.}
  \item{image}{Set to \code{TRUE} if you want an image plot underneath the contours.}
  \item{img.col}{Colr mapn to use when \code{image=TRUE}.}
  \item{\dots}{Additional arguments passed to \code{\link{contour}}.}
}
\details{
  The first two elements of \code{varlist} define the primary variables.  
  In each of the primary variables, you may specify:
%%%%  \arguments{
    \item{NULL}{A default grid of 26 values is used}
    \item{A single value}{Used as the number of points on the grid}
    \item{Vector of 2 values}{Used as minimum and maximum}
    \item{Vector of 3 values}{Used as minimum, maximum, and number of points}
    \item{Longer vector}{Used as the grid of values for that coordinate}
%%%%  }
  In most cases, you will want to specify single values for all subsequent entries,
  to specify values to hold fixed.  By default, the mean of each other numerical
  variable, and the first level of each factor variable is used.
}
\value{
  Nothing is returned
}
\author{Russell V. Lenth}
\seealso{\code{\link{contour}}}
\examples{
library(rsm)
CR = coded.data (ChemReact, x1~(Time-85)/5, x2~(Temp-175)/5)
CR.rs2 = rsm (Yield ~ Block + SO(x1,x2), data=CR) 
contour(CR.rs2, list(x1=NULL, x2=NULL))
}
\keyword{ regression }
\keyword{ hplot }
