\name{CARBayesST-package}
\alias{CARBayesST-package}
\alias{CARBayesST}
\docType{package}

\title{
Spatio-Temporal Generalised Linear Mixed Models For Areal Unit Data
}

\description{
Implements a class of spatio-temporal generalised linear mixed models for areal unit data, 
with inference in a Bayesian setting using Markov chain Monte Carlo (McMC) simulation. The 
response variable can be binomial, Gaussian or Poisson. The spatio-temporal autocorrelation
is modelled by  random effects, which are assigned conditional autoregressive (CAR) prior 
distributions. A number of different random effects structures are available, and full details
are given in the vignette accompanying this package and the references below. The initial 
creation the Engineering and Physical Science Research Council (EPSRC) grant EP/J017442/1.
}

\details{
\tabular{ll}{
Package: \tab CARBayesST\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-07-02\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
Author: Duncan Lee and Alastair Rushworth

Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}

\references{
Bernardinelli, L., D. Clayton, C.Pascuto, C.Montomoli, M.Ghislandi, and M. Songini
(1995). Bayesian analysis of space-time variation in disease risk. Statistics in 
Medicine, 14, 2433-2443.

Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time variation in disease risk. Statistics in Medicine, 19, 2555-2567.

Lee, D and Sarran, C  (2015). Controlling for unmeasured confounding and spatial 
misalignment in long-term air pollution and health studies, Environmetrics, to appear.
    
Rushworth, A., Lee, D., and Sarran, C.  An adaptive spatio-temporal smoothing model for estimating trends and step changes in disease risk (2014a) arXiv:1411.0924

Rushworth, A., D. Lee, and R. Mitchell (2014b). A spatio-temporal model for estimating the long-term effects of air pollution on respiratory hospital admissions in Greater London. Spatial and Spatio-temporal Epidemiology 10, 29-38.
}


\examples{
## See the examples in the function specific help files and in the vignette
## accompanying this package.
}
