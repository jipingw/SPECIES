\name{jackknife}
\alias{jackknife}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Jackknife estimator for the species richness}
\description{
 A function implementing the jackknife estimator of the species 
number by Burnham and Overton 1978 and 1979.
}
\usage{
jackknife(n, k = 5, conf = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a matrix or a numerical data frame of two columns. It is also called the \dQuote{frequency of frequencies} data in literature. The first column is the frequency \eqn{j=1, 2\ldots}; and the second column is \eqn{n_j}, the  number of species  observed with \eqn{j}
 individuals in the sample.  }
  \item{k}{a positive integer. \code{k} is the specified Jackknife order. The default is \code{k}=5. Burnham and Overton 1978 and 1979 provided a testing procedure for the maximum order to be used in this estimator. If the specified order \code{k} or default is greater than the 
           order obtained from the testing procedure, the function will automatically use the determined order rather than \code{k}. 
	   Currently this function only provide jackknife estimate up to order 10.}
    \item{conf}{a positive number \eqn{\le 1}. \code{conf} specifies the confidence level for confidence interval. The default is 0.95.
                \code{conf} also specifies the critical value in the sequential test for jackknife order.
}

}
%\details{}
\value{
   The function \code{jackknife} returns a list of:  \code{JackknifeOrder}, \code{Nhat}, \code{SE} and \code{CI}.
  \item{JackknifeOrder}{the jackknife estimator order specified order by the user or determined by the testing procedure.}
  \item{Nhat}{jackknife estimate.}
  \item{SE}{standard error of the jackknife estimate.}
  \item{CI}{confidence interval of the jackknife estimate.}
}
\references{Burnham, K. P., and Overton,W. S. (1978), Estimation of the Size of a Closed Population When Capture Probabilities Vary Among Animals, Biometrika, 65, 625-633.

Burnham, K. P., and Overton,W. S. (1979), Robust Estimation of Population Size When Capture Probabilities Vary Among Animals, Ecology, 60, 927-936.}

\author{Ji-Ping Wang, Department of Statistics, Northwestern University}
%\note{}
%\seealso{ }
\examples{

library(SPECIES)

##load data from the package, 
## "butterfly" is the famous tterfly data by Fisher 1943.

data(butterfly)
jackknife(butterfly,k=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}

