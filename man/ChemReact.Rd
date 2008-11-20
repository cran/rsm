\name{ChemReact}
\alias{ChemReact}
%%%\docType{data}
\title{Chemical Reaction Data}
\description{
  These are uncoded data (in their original units) from a central
  composite design with 2 factors in 2 blocks.
}
\usage{ChemReact}
\format{
  A data frame with 14 observations on the following 4 variables.
  \describe{
    \item{\code{Time}}{a numeric vector; predictor variable}
    \item{\code{Temp}}{a numeric vector; predictor variable}
    \item{\code{Block}}{a factor with levels \code{B1} \code{B2}.
    Block \code{B1} is a first-order design with 3 center points.  
    Block \code{B2} consists of axis points and 3 more center points.}
    \item{\code{Yield}}{a numeric vector; response variable}
  }
}
\source{
  Table 7.7 of Myers, RH and Montgomery, DC (2002), \emph{Response Surfaces} (2nd ed),
  Wiley-Interscience.
}
\keyword{datasets}
