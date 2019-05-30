\name{traffic}
\alias{traffic}
\docType{data}
\title{Traffic data}
\description{The traffice data originally appeared in Simar 1976 where the total number of N is knowns as 9461. Recently
it was re-analyzed by Bohning and Schon 2005.}
%\author{Ji-Ping Wang, Department of Statistics, Northwestern University, Maintainer: \email{jzwang@northwestern.edu}}
\references{
Simar, L. (1976), Maximum likelihood estimation of a compound Poisson process, Annals of Statistics, 4, 1200-1209.
Bohning, D., and Schon, D. (2005), Nonparametric maximum likelihood estimation of population size based on the counting distribution, Journal of the Royal Statistical Society, Series C: Applied Statistics, 54, 721-737.
}

\examples{
##load library
library(SPECIES)

##load data that coming with the package.
data(traffic)
chao1984(traffic)
}
