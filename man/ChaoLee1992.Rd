\name{ChaoLee1992}
\alias{ChaoLee1992}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coverage-based estimators for species richness}
\description{
  This function calculates \code{ACE} and \code{ACE-1} estimators by Chao and Lee 1992 (\code{ACE-1} provides
further bias correction based on \code{ACE}).
}
\usage{
       ChaoLee1992(n, t = 10, method = "all",conf = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a matrix or a numerical data frame of two columns. It is also called the \dQuote{frequency of frequencies} data in literature. The first column is the frequency \eqn{j=1, 2\ldots}; and the second column is \eqn{n_j}, the  number of species  observed with \eqn{j}
 individuals in the sample. }
  \item{t}{a positive integer. \code{t} is the cutoff value to define the relatively less abundant species to be used in estimation. The frequencies \eqn{n_j} of \eqn{j>t} will not be used in estimating the sample coverage. The default value is \code{t}=10.}
  \item{method}{a string. It can be any one of  \dQuote{ACE}, \dQuote{ACE-1}, or \dQuote{all}. The default is \dQuote{all}.}
    \item{conf}{a positive number \eqn{\le 1}. \code{conf} specifies the confidence level for confidence interval. The default is 0.95.}

}
%\details{}
\value{
 
  The function \code{ChaoLee1992} returns a list of:  \code{Nhat}, \code{SE} and \code{CI}.
 
  \item{Nhat}{point estimate of the specified method. If the default \code{method}=\dQuote{all} is used, the function returns an estimate vector including \code{ACE, ACE-1}.}
  \item{SE}{standard error(s) of the point estimate(s).}
  \item{CI}{confidence interval using a log transformation explained in Chao 1987.}
}

\references{
Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics
43, 783-791.

Chao, A., and Lee, S.-M. (1992), Estimating the Number of Classes via Sample Coverage, Journal of the American Statistical Association, 87, 210-217.

}
\author{Ji-Ping Wang, Department of Statistics, Northwestern University }
%\note{}
%\seealso{}

\examples{
library(SPECIES)

##load data from the package, 
## "butterfly" is the famous butterfly data by Fisher 1943.

data(butterfly)

##output estimates from all 4 methods using cutoff t=10
ChaoLee1992(butterfly,t=10,method="all")

##output estimates from ACE method using cutoff t=10
ChaoLee1992(butterfly,t=10,method="ACE")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

