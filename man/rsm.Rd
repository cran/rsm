\name{rsm}
\alias{rsm}
\alias{summary.rsm}
\alias{print.summary.rsm}
\alias{loftest}
\title{Response-surface regression}
\description{
  Fit a linear model with a response-surface component,
  and produce appropriate analyses and summaries.
}
\usage{
rsm (..., data)
\method{summary}{rsm} (object, ...)
\method{print}{summary.rsm}(x, ...)
loftest (object)
}
\arguments{
  \item{\dots}{In \code{rsm}, arguments that are passed to \code{\link{lm}}.
  The model must include an \code{FO()} or \code{SO()} term.to
  define the response-surface portion of the model.
  
  In \code{summary}, and \code{print}, additional arguments 
  are passed to their generic methods.
  }
  \item{data}{Optional \code{data.frame} as in \code{lm}.  This could also be
  a \code{\link{coded.data}} object, in which case additional information
  is included in subsequent summaries.}
  \item{object}{An object of class \code{rsm}}
  \item{x}{An object produced by \code{summary}}
}
\details{
  In \code{rsm}, the model formula must contain at least an \code{FO} term; optionally, you can add
  a \code{TWI()} term and/or a \code{PQ()} term as well (use the same variables in each!).
  For convenience, specifying \code{SO()} is the same as including \code{FO()}, \code{TWI()}, and \code{PQ()},
  and is the safe, preferred way of specifying a full second-order model.
}
\value{
  \code{rsm} returns an \code{rsm} object, which is a \code{lm} object with 
  additional members as follows:
    \item{order}{The order of the model: 1 for first-order, 1.5 for first-order plus interactions, 
     or 2 for a model that contains square terms.}
    \item{b}{The first-order response-surface coefficients.}
    \item{B}{The matrix of second-order response-surface coefficients, if present.}
    \item{labels}{Labels for the response-surface terms.  These make the summary much more readable.}
    \item{coding}{Coding formulas, if \code{data} is a \code{\link{coded.data}} object.}
  \code{summary} 
  is the summary method for \code{rsm} objects.  It returns an object of class
  \code{summary.rsm}, which is an extension of the \code{summary.lm}
  class with additional list elements:
    \item{sa}{Unit-length vector of the path of steepest ascent
      (first-order models only).}
    \item{canonical}{Canonical analysis (second-order models only).  This is a list with elements \code{xs}, the stationary point, and \code{eigen}, the eigenanalysis
      of \code{B} (see above).}
    \item{lof}{ANOVA table including lack-of-fit test.}
    \item{coding}{Coding formulas in parent \code{rsm} object.}
  Its \code{print} method shows the regression summary,
  followed by an ANOVA and lack-of-fit test.  For first-order models, it shows the direction of 
  steepest ascent, and for second-order models, it shows the canonical analysis of the
  response surface.  
  
  \code{loftest} returns an \code{\link{anova}} object that tests the fitted model against a model 
  that interpolates the means of the response-surface-variable combinations.
}
%%%\references{}
\author{Russell V. Lenth}
%%%\note{ ~~further notes~~ 
%%%
%%% ~Make other sections like Warning with \section{Warning }{....} ~
%%%}
\seealso{\code{\link{FO}}, \code{\link{SO}},
  \code{\link{lm}}, \code{\link[=lm-class]{summary}}, \code{\link{coded.data}}}
\examples{
library(rsm)
CR = coded.data (ChemReact, x1~(Time-85)/5, x2~(Temp-175)/5)

### 1st-order model, using only the first block
CR.rs1 = rsm (Yield ~ FO(x1,x2), data=CR, subset=1:7) 
summary(CR.rs1)

### 2nd-order model, using both blocks
CR.rs2 = rsm (Yield ~ Block + SO(x1,x2), data=CR) 
summary(CR.rs2)
}
\keyword{regression}
