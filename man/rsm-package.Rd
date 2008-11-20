\name{rsm-package}
\alias{rsm-package}
\docType{package}
\title{
Response-surface analysis
}
\description{
The \code{rsm} package provides functions useful for analyzing
experiments that are done sequentially in hopes of optimizing a response surface.
The function \code{\link{rsm}} is an enhancement of \code{\link{lm}} that provides 
for additional analyses peculiar to response surfaces.  It requires a model formula 
that contains a call to \code{\link{FO}} or \code{\link{SO}} to specify a first- or 
second-order model.  Once the model is fitted, the \code{\link{steepest}} 
function may be used to obtain the direction of steepest ascent (or descent).  
\code{\link{canonical.path}} is an alternative to \code{steepest} for second-order 
response surfaces.

In RSM methods, appropriate coding of data is
important not only for nyumerical stability, but for proper scaling
of results; the function \code{\link{coded.data}} and its relatives facilitate
this coding requirement.  

Finally, two more functions are provided that may be useful beyond response-surface applications.  
\code{\link{contour.lm}} aids in visualizing a response surface, 
or of any other \code{lm} object where a surface is fitted.  \code{\link{model.data}}
recovers the data used in a \code{lm} call, but unlike \code{model.frame}, no
polynomials, factors, etc. are expanded.
}
\details{
\tabular{ll}{
Package: \tab rsm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-11-20\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
Russell V. Lenth

Maintainer: Russ Lenth <russell-lenth@uiowa.edu>
}
\references{
Box, GEP, Hunter, JS, and Hunter, WG (2005), 
\emph{Statistics for Experimenters} (2nd ed.), Wiley-Interscience.

Meyers, RH and Montgomery, DC (2002), 
\emph{Response Surface Methodology} (2nd ed.), Wiley-Interscience.
}
\keyword{ package }
\keyword{ regression }

%%%\seealso{
%%%~~ Optional links to other man pages, e.g. ~~
%%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%%}
