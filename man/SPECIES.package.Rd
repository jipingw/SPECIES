\name{SPECIES-package}
\alias{SPECIES-package}
\alias{SPECIES}
\docType{package}
\title{An R pacakge for species richness estimation}
\description{
\code{SPECIES} provides multiple functions to compute popular estimators for species richness.These estimators include:
(1) jackknife estimator by Burnham and Overton 1978, 1979; (2) lower-bound estimator by Chao 1984; (3) coverage-base estimators ACE, ACE-1 by Chao and Lee 1992; (4) coverage-duplication estimator from Poisson-Gamma model by Chao and Bunge 2002; (5) unconditional nonparametric maximum likelihood estimator by Norris and Pollock 1996, 1998; (6) penalized nonparametric maximum likelihood estimator by Wang and Lindsay 2005; and (7) Poisson-compound Gamma model with smooth nonparametric maximum likelihood estimation by Wang 2010.
}
\details{
\tabular{ll}{
Package: \tab SPECIES\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-1-22\cr
License: \tab GPL-2\cr
}
functions: chao1984, ChaoBunge, ChaoLee1992, jackknife, pcg ,pnpmle, unpmle;
data: butterfly, cottontail,EST, insect, microbial, traffic
}
\author{
Ji-Ping Wang, Department of Statistics, Northwestern University

Maintainer: jzwang@northwestern.edu
}
\references{

Acinas, S., Klepac-Ceraj, V., Hunt, D., Pharino, C., Ceraj, I., Distel, D., and Polz, M. (2004), Fine-scale phylogenetic architecture of a complex bacterial community. Nature, 430, 551-554.

Bohning, D. and Schon, D., Nonparametric maximum likelihood estimation of population size based on the counting distribution, Journal of the Royal Statistical Society,
Series C: Applied Statistics, 54, 721-737.

Burnham, K. P., and Overton,W. S. (1978), Estimation of the Size of a Closed Population When Capture Probabilities Vary Among Animals, Biometrika, 65, 625-633.

Burnham, K. P., and Overton,W. S. (1979), Robust Estimation of Population Size When Capture Probabilities Vary Among Animals, Ecology, 60, 927-936.

Chao, A. (1984), Nonparametric Estimation of the Number of Classes in a Population, Scandinavian Journal of Statistics, 11, 265-270.

Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.

Chao, A., and Lee, S.-M. (1992), Estimating the Number of Classes via Sample Coverage, Journal of the American Statistical Association, 87, 210-217.

Chao, A., and Bunge, J. (2002), Estimating the Number of Species in a Stochastic Abundance Model, Biometrics, 58, 531-539.

Fisher, R. A., Corbet, A. S., and Williams, C. B. ,(1943), The Relation Between the Number of Species and the Number of Individuals in a Random Sample of an Animal Population, Journal of Animal Ecology, 12, 42-58.

Hong, S. H., and Bunge, J. and Jeon, S.O. and Epstein, S. (2006), Predicting microbial species richness, Proc. Natl. Acad. Sci, 103, 117-122.

Norris, J. L. I., and Pollock, K. H. (1996), Nonparametric MLE Under Two Closed Capture-Recapture Models With Heterogeneity, Biometrics, 52,639-649. 

Norris, J. L. I., and Pollock, K. H.(1998), Non-Parametric MLE for Poisson Species Abundance Models Allowing for Heterogeneity Between Species, Environmental and Ecological Statistics, 5, 391-402.

Simar, L. (1976), Maximum likelihood estimation of a compound Poisson process, Annals of Statistics, 4, 1200-1209.


Wang, J.-P. Z. and Lindsay, B. G. (2005), A penalized nonparametric maximum likelihood approach to species richness estimation. Journal of American Statistical Association, 100(471):942-959.

Wang, J.-P.,  and Lindsay, B.G. (2008), An exponential partial prior for improving NPML estimation for mixtures, Statistical Methodology, 5:30-45. 

Wang, J.-P. (2010), Estimating the species richness by a Poisson-Compound Gamma model, Biometrika, 97(3): 727-740.

Wang, J.-P. (2011), SPECIES: An R Package for Species Richness Estimation, Journal of Statistical Software, 40(9), 1-15, URL: http://www.jstatsoft.org/v40/i09/.
}

\examples{
##load library
library(SPECIES)

## "butterfly" is the famous butterfly data by Fisher 1943.

data(butterfly)

##jackknife method
jackknife(butterfly,k=5)

##using only 'ACE' coverage method
ChaoLee1992(butterfly,t=10, method="all")

##using chao1984 lower bound estimator
chao1984(butterfly)

##using Chao and Bunge coverage-duplication method
ChaoBunge(butterfly,t=10)

##penalized NPMLE method
#pnpmle(butterfly,t=15,C=1,b=200)

##unconditonal NPMLE method
#unpmle(butterfly,t=10,C=1,b=200)

##Poisson-compound Gamma method
#pcg(butterfly,t=20,C=1,b=200)
}
