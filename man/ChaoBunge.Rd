\name{ChaoBunge}
\alias{ChaoBunge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coverage-duplication estimator for species richness}

\description{
  This function calculates coverage-duplication based estimator from a Poisson-Gamma model by Chao and Bunge 2002.
}
\usage{
       ChaoBunge(n, t = 10,conf = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a matrix or a numerical data frame of two columns. It is also called the \dQuote{frequency of frequencies} data in literature. The first column is the frequency \eqn{j=1, 2\ldots}; and the second column is \eqn{n_j}, the  number of species  observed with \eqn{j}
 individuals in the sample. }
  \item{t}{a positive integer. \code{t} is the cutoff value to define the relatively less abundant species to be used in estimation. The frequencies \eqn{n_j} of \eqn{j>t} will not be used in estimating the sample coverage. The default value is \code{t}=10.}
    \item{conf}{a positive number \eqn{\le 1}. \code{conf} specifies the confidence level for confidence interval. The default is 0.95.}

}
%\details{}
\value{
 
  The function \code{ChaoBunge} returns a list of:  \code{Nhat}, \code{SE} and \code{CI}.
 
  \item{Nhat}{point estimate.}
  \item{SE}{standard error(s) of the point estimate.}
  \item{CI}{confidence interval using a log transformation explained in Chao 1987.}
}

\references{Chao, A. (1984), Nonparametric Estimation of the Number of Classes in a Population, Scandinavian Journal of Statistics, 11, 265-270.

Chao, A., and Bunge, J. (2002), Estimating the Number of Species in a Stochastic Abundance Model, Biometrics, 58, 531-539.
}
\author{Ji-Ping Wang, Department of Statistics, Northwestern University }
%\note{}
%\seealso{}

\examples{
library(SPECIES)

##load data from the package, 
##"butterfly" is the famous butterfly data by Fisher 1943.

data(butterfly)

##output estimates from all 4 methods using cutoff t=10
ChaoBunge(butterfly,t=10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

