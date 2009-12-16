\name{coded.data}
\alias{coded.data}
\alias{as.coded.data}
\alias{decode.data}
\alias{val2code}
\alias{code2val}
\alias{print.coded.data}
\alias{parse.coding}
\alias{codings}
\alias{codings.coded.data}
\title{Functions for coded data}
\description{
  These functions facilitate the use of coded data in response-surface analysis.
}
\usage{
coded.data(data, ..., formulas = list(...))
as.coded.data(data, ..., formulas = list(...))

decode.data(data)
val2code(X, codings)
code2val(X, codings)

\method{print}{coded.data}(x, ...)

codings(object)
\method{codings}{coded.data}(object)
parse.coding(form)
}
\arguments{
  \item{data}{A data.frame}
  \item{\dots}{Formulas for producing coded variables.}
  \item{formulas}{Provides an alternative way to provide the coding formulas.
  This is useful if you want to retrieve the codings from one dataset and use them for another.}
  \item{X}{A vector, matrix, or data.frame to be coded or decoded.}
  \item{codings}{a list of formulas (see \code{form})}
  \item{x}{Coded dataset to be printed.}
  \item{form}{A formula as described above for \code{\dots}}
  \item{object}{An object that contains coding formulas}
}
\details{
  Typically, coding formulas are of the form \code{x ~ (var - center) / div} 
  where \code{x} and \code{var} are variable names, and \code{center} and \code{div} are numbers.  
  The left-hand side gives the name of the coded variable, and the right-hand side 
  should be a linear expression in  the uncoded variable.  This expression is evaluated at
  variable values of 0 and 1, then the results are used to
  solve for the scale center and divisor.  These results are rounded to 4 digits
  to help ensure that zeros come out exactly.
}
\value{
  \code{coded.data} and \code{as.coded.data} return an object of class \code{coded.data},
  which inherits from \code{\link{data.frame}}.  \code{print.coded.data} is the print 
  method for this class; it simply prints the data.frame and then the formulas.
  A named \code{\link{list}} of the coding formulas is saved in 
  \code{attr(, "codings")}.
  
  Use \code{coded.data} to convert a data.frame in which the variables
  are on their original scales.  The variables named in the formulas are 
  coded and replaced with their coded versions (and also renamed).
  In contrast, \code{as.coded.data} does not modify any of the data; it assumes the variables
  are already coded, and the coding information is simply added.
  
  \code{decode.data} converts a dataset of class \code{coded.data} and returns \code{data.frame} 
  containing the original variables.
  
  \code{code2val} converts coded values to the original scale using the codings provided, 
  and returns an object of the same class as \code{X}.  
  \code{val2code} converts the other direction.  When using these functions, it is 
  essential that the names (or column names in the case of matrices) match those of the corresponding 
  coded or uncoded variables.
  
  \code{parse.coding} is mostly for internal use; it extracts and returns a \code{list} with two vectors:
  a character vector \code{names} with the two variable names, 
  and a numeric vector \code{const} with the center and divisor.
  
  \code{codings} is a generic function for accessing codings.  \code{codings.coded.data} returns the 
  list of coding formulas from a \code{coded.data} object.  See also \code{\link{codings.rsm}}.
}

\seealso{\code{\link{data.frame}}}
\references{
Lenth RV (2009). ``Response-Surface Methods in R, Using rsm'', 
  \emph{Journal of Statistical Software}, 32(7), 1--17. 
  \url{http://www.jstatsoft.org/v32/i07/.}
}
\author{Russell V. Lenth}
\examples{
library(rsm)

CR = coded.data (ChemReact, x1~(Time-85)/5, x2~(Temp-175)/5)
CR
decode.data (CR)  ### will be same as ChemReact

code2val (c(x1=.5, x2=-1), codings = codings(CR))
}
\keyword{regression}
