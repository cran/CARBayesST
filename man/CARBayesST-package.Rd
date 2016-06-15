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
is modelled by  random effects,  which are assigned conditional autoregressive (CAR) style prior 
distributions. A number of different random effects structures are available, and full details 
are given in the vignette accompanying this package and the references below. The creation of 
this package was supported by the Engineering and Physical Sciences Research Council (EPSRC) 
grant EP/J017442/1 and the Medical Research Council (MRC) grant MR/L022184/1.

The main change between version 2.3 and 2.2 is a vastly enhanced vignette with worked examples.
Additionally, a few bugs in the code have been fixed. Version 2.2 itself has the following 
changes from version 2.1. Firstly, missing values are allowed  in the response variable for 
the models ST.CARanova(), ST.CARar() and ST.CARlinear(). Secondly, the spatial and temporal 
dependence parameters (denoted rho) in ST.CARanova(), ST.CARar(), ST.CARlinear() and 
ST.CARsepspatial() are allowed to be fixed at values in the unit interval [0,1] or estimated 
in the model. Thirdly, the modelfit component of the  fitted model object now additionally 
returns the  Watanabe-Akaike Information Criterion (WAIC) and an estimate of the effective 
number of effective parameters (p.w).
}

\details{
\tabular{ll}{
Package: \tab CARBayesST\cr
Type: \tab Package\cr
Version: \tab 2.3\cr
Date: \tab 2016-06-15\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
Author: Duncan Lee, Alastair Rushworth and Gary Napier

Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}

\references{
Bernardinelli, L., D. Clayton, C.Pascuto, C.Montomoli, M.Ghislandi, and M. Songini
(1995). Bayesian analysis of space-time variation in disease risk. Statistics in 
Medicine, 14, 2433-2443.

Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time variation in disease risk. 
Statistics in Medicine, 19, 2555-2567.

Lee, D and Lawson, C  (2016). Quantifying the spatial inequality and temporal trends in
maternal smoking rates in Glasgow, Annals of Applied Statistics, to appear.

Napier, G, D. Lee, C. Robertson, A. Lawson, and K. Pollock (2016). A model to estimate the 
impact of changes in MMR vaccination uptake on inequalities in measles susceptibility in 
Scotland, Statistical Methods in Medical Research,to appear.
        
Rushworth, A., D. Lee, and R. Mitchell (2014). A spatio-temporal model for estimating the 
long-term effects of air pollution on respiratory hospital admissions in Greater London. 
Spatial and Spatio-temporal Epidemiology 10, 29-38.

Rushworth, A., Lee, D., and Sarran, C (2016).  An adaptive spatio-temporal smoothing model for 
estimating trends and step changes in disease risk. Journal of the Royal Statistical Society 
Series C, DOI: 10.1111/rssc.12155.
}


\examples{
## See the examples in the function specific help files and in the vignette
## accompanying this package.
}
