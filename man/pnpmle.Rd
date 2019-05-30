\name{pnpmle}
\alias{pnpmle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Penalized conditional NPML estimator for species richness}
\description{
This function calculate the penalized conditional NPML estimator of the species number by Wang and Lindsay 2005.
This estimator was based on the conditional likelihood of a Poisson mixture model. A penalty term was 
introduced into the  model to prevent the boundary problem discussed in Wang and Lindsay 2008. The confidence interval is 
calculated based on a bootstrap procedure. A Fortran function is called to for the computing.
}
\usage{
pnpmle(n,t=15,C=0,b=200,seed=NULL,conf=0.95,dis=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a matrix or a numerical data frame of two columns. It is also called the \dQuote{frequency of frequencies} data in literature. The first column is the frequency \eqn{j=1, 2\ldots}; and the second column is \eqn{n_j}, the  number of species  observed with \eqn{j}
 individuals in the sample. }
  \item{t}{a positive integer. \code{t} is the cutoff value to define the relatively less abundant species to be used in estimation of the Poisson mixture. The default value is \code{t}=15. The recommendation is to use \eqn{\code{t} \ge 10}.}
  \item{C}{integer either 0 or 1. It specifies whether bootstrap confidence interval should be calculated. \dQuote{\code{C}=1} for YES and \dQuote{\code{C}=0} for NO.The default of \code{C} is set as 0.}
 \item{b}{integer. \code{b} specifies the number of bootstrap samples to be generated for confidence interval. It is ignored if \dQuote{\code{C}=0}.}
    \item{conf}{a positive number \eqn{\le 1}. \code{conf} specifies the confidence level for confidence interval. The default is 0.95.}
  \item{seed}{a single value, interpreted as an integer. Seed for random number generation}
  \item{dis}{0 or 1. 1 for on-screen display of the mixture output, and 0 for none.}
}
%\details{}
\value{
  The function \code{pnpmle} returns a list of:  \code{Nhat}, \code{CI} (if \dQuote{\code{C}=1}).
  \item{Nhat}{Point estimate of \code{N}}
  \item{CI}{bootstrap confidence interval}
}
\references{ Wang, J.-P. Z. and Lindsay, B. G. ,2005, A penalized nonparametric maximum likelihood approach to species richness estimation. Journal of American Statistical Association, 2005,100(471):942-959

Wang, J.-P.,  and Lindsay, B.G., 2008, An exponential partial prior for improving NPML estimation for mixtures, Statistical Methodology, 2008,5:30-45  
}
\author{Ji-Ping Wang,Department of Statistics, Northwestern University}
%\note{ }
%\seealso{ }
\examples{
library(SPECIES)

##load data from the package, 
## \dQuote{butterfly} is the famous butterfly data by Fisher 1943.
#data(butterfly)

##output estimate without confidence interval using cutoff t=15
#pnpmle(butterfly,t=15,C=0)

##output estimate with confidence interval using cutoff t=15
#pnpmle(butterfly,t=15,C=1, b=200)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

