\name{pcg}
\alias{pcg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Poisson-compound Gamma estimator for the species richness }

\description{
Function to calculate the Poisson-compound Gamma estimators of the species number by Wang 2010. This method
is essentially a conditional NPMLE method. The species abundance here is assumed to follow a compound 
Gamma model. The confidence interval is obtained based on a bootstrap procedure.
A Fortran function is called to for the computing. This function requires Fortran compiler installed.
}
\usage{
pcg(n,t=35,C=0,alpha=c(1:10),b=200,seed=NULL,conf=0.95,dis=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{a matrix or a numerical data frame of two columns. It is also called the \dQuote{frequency of frequencies} data in literature. The first column is the frequency \eqn{j=1, 2\ldots}; and the second column is \eqn{n_j}, the  number of species  observed with \eqn{j}
 individuals in the sample.  }
 \item{t}{a positive integer. \code{t} is the cutoff value defining the relatively less abundant species to be used in estimation. The default value for \code{t}=35.
 The  estimator is more sensitive to t compared with pnpmle or unpmle estimators. We recommend to use \eqn{t \ge 20} if the maximum frequency (\eqn{j}) is greater than 20. Otherwise use the maximum frequency of \eqn{j} for \code{t}.}
 \item{C}{integer either 0 or 1. It specifies whether bootstrap confidence interval should be calculated. \dQuote{\code{C}=1} for YES and \dQuote{\code{C}=0} for NO.The default of \code{C} is set as 0.}
  \item{b}{integer. \code{b} specifies the number of bootstrap samples to be generated for confidence interval. It is ignored if \dQuote{\code{C}=0}.}
  \item{alpha}{a positive grid for Gamma shape parameter. \code{alpha} must be a numerical vector for positive numbers. A cross-validation will be used to select a unified shape parameter value for the compound Gamma from the specified \dQuote{\code{alpha}} grid. The default \dQuote{\code{alpha}} grid is \eqn{1,2,\ldots,10}.}
    \item{conf}{a positive number \eqn{\le 1}. \code{conf} specifies the confidence level for confidence interval. The default is 0.95.}
  \item{seed}{a single value, interpreted as an integer. Seed for random number generation}
  \item{dis}{0 or 1. 1 for on-screen display of the mixture output, and 0 for none.}

}
\details{
 The \code{pcg} estimator is computing intensive. The computing of bootstrap confidence interval may take up to a few hours.
}
\value{
 The function \code{pcg} returns a list of:  \code{Nhat}, \code{CI} (if \dQuote{\code{C}=1}) and \code{AlphaModel}.
  \item{Nhat}{point estimate of \code{N}.}
  \item{CI}{bootstrap confidence interval.}
  \item{AlphaModel}{unified shape parameter of compound Gamma selected from cross-validation.}
}
\references{Wang, J.-P. (2010), Estimating the species richness by a Poisson-Compound Gamma model, 97(3): 727-740}
\author{Ji-Ping Wang, Department of Statistics, Northwestern University}
%\note{}
%\seealso{}
\examples{
library(SPECIES)
##load data from the package, 
## \dQuote{butterfly} is the famous butterfly data by Fisher 1943.

data(butterfly)

##output estimate without confidence interval using cutoff t=15
##pcg(butterfly,t=20,C=0,alpha=c(1:10))

##output estimate with confidence interval using cutoff t=15
#pcg(butterfly,t=20,C=1,alpha=c(1:10),b=200)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

