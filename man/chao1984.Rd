\name{chao1984}
\alias{chao1984}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lower-bound estimator for species richness}
\description{
  This function calculates the lower-bound estimator by Chao 1984.}
\usage{
       chao1984(n,conf=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a matrix or a numerical data frame of two columns. It is also called the \dQuote{frequency of frequencies} data in literature.
 The first column is the frequency \eqn{j=1, 2\ldots}; and the second column is \eqn{n_j}, the  number of species  observed with \eqn{j}
 individuals in the sample. }
    \item{conf}{a positive number \eqn{\le 1}. \code{conf} specifies the confidence level for confidence interval. The default is 0.95.}

}
%\details{}
\value{
 
  The function \code{chao1984} returns a list of:  \code{Nhat}, \code{SE} and \code{CI}.
 
  \item{Nhat}{point estimate.}
  \item{SE}{standard error of the point estimate.}
  \item{CI}{confidence interval using a log transformation explained in Chao 1987.}
}

\references{Chao, A. (1984), Nonparametric Estimation of the Number of Classes in a Population, Scandinavian Journal of Statistics, 11, 265-270.

Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics
43, 783-791.
}
\author{Ji-Ping Wang, Department of Statistics, Northwestern University }
%\note{}
%\seealso{}

\examples{
library(SPECIES)

##load data from the package, 
## \dQuote{butterfly}is the famous butterfly data by Fisher 1943.

data(butterfly)
chao1984(butterfly)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

