\name{unpmle}
\alias{unpmle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Unconditional NPML estimator for the SPECIES number}
\description{
 This function calculate the unconditional NPML estimator of the species number by Norris and Pollock 1996, 1998.
This estimator was obtained from the full likelihood based on a Poisson mixture model. The confidence interval is 
calculated based on a bootstrap procedure.
}
\usage{
unpmle(n,t=15,C=0,method="W-L",b=200,conf=.95,seed=NULL,dis=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a matrix or a numerical data frame of two columns. It is also called the \dQuote{frequency of frequencies} data in literature. The first column is the frequency \eqn{j=1, 2\ldots}; and the second column is \eqn{n_j}, the  number of species  observed with \eqn{j}
 individuals in the sample.   }
  \item{t}{a positive integer. \code{t} specifies the cutoff value to define the relatively less abundant species to be used in estimation.  The default value for t=15. The estimator is fairly insensitive to the choice of \code{t}. The recommendation is
to use \eqn{t \ge 10}.}
  \item{C}{integer either 0 or 1. It specifies whether bootstrap confidence interval should be calculated. \dQuote{\code{C}=1} for YES and \dQuote{\code{C}=0} for NO.The default of \code{C} is set as 0.}
  \item{method}{string either \dQuote{N-P} or \dQuote{W-L}(default). If \code{method}=\dQuote{N-P}, unconditional NPMLE will be used using an algorithm by  Bonhing and Schon (2005). 
   Sometimes this method can be extremely slow. Alternatively one can use method \dQuote{W-L}, an approximate method (but with high precision and much faster) by Wang and Lindsay 2005.}
   \item{b}{integer. \code{b} specifies the number of bootstrap samples for confidence interval. It is ignored if \dQuote{\code{C}=0}.}
    \item{conf}{a positive number \eqn{\le 1}. \code{conf} specifies the confidence level for confidence interval. The default is 0.95.}
  \item{seed}{a single value, interpreted as an integer. Seed for random number generation}
  
  \item{dis}{0 or 1. 1 for on-screen display of the mixture output, and 0 for none.}
}
\details{The computing is intensive if \code{method}=\dQuote{N-P} is used particularly when extrapolation is large. 
It may takes hours to compute the bootstrap confidence interval. If \code{method}=\dQuote{W-L} is used, computing usually 
is much much faster. Estimates from both methods are often identical.}
\value{
 The function \code{unpmle} returns a list of:  \code{Nhat}, \code{CI} (if \dQuote{\code{C}=1})
  \item{Nhat}{point estimate of N}
  \item{CI}{bootstrap confidence interval.}
}

\references{ 
Norris, J. L. I., and Pollock, K. H. (1996), Nonparametric MLE Under Two Closed Capture-Recapture Models With Heterogeneity, Biometrics, 52,639-649. 

Norris, J. L. I., and Pollock, K. H.(1998), Non-Parametric MLE for Poisson Species Abundance Models Allowing for Heterogeneity Between Species, Environmental and Ecological Statistics, 5, 391-402.

Bonhing, D. and Schon, D.,  (2005), Nonparametric maximum likelihood estimation of population size based on the counting distribution, Journal of the Royal Statistical Society,
Series C: Applied Statistics, 54, 721-737.

Wang, J.-P. Z. and Lindsay, B. G. ,(2005), A penalized nonparametric maximum likelihood approach to species richness estimation. Journal of American Statistical Association, 2005,100(471):942-959
}
\note{The unconditional NPML estimator is unstable from either \code{method='N-P'} or \code{method='W-L'}. Extremely large estimates may occur. 
This is also reflected in that the upper confidence bound often greatly vary from different runs of bootstrap procedure. In contrast the penalized NPMLE by \code{pnpmle} function is much more stable.}
%\seealso{}
\author{Ji-Ping Wang, Department of Statistics, Northwestern University}


\examples{
library(SPECIES)

##load data from the package, 
## "butterfly" is the famous butterfly data by Fisher 1943.

data(butterfly)


##output estimate without confidence interval using cutoff t=15
#unpmle(butterfly,t=15,C=0)

##output estimate with confidence interval using cutoff t=15
#unpmle(butterfly,t=15,C=1,b=200)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

