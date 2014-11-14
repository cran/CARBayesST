\name{CARBayesST-package}
\alias{CARBayesST-package}
\alias{CARBayesST}
\docType{package}

\title{
Poisson Log-linear Models with Spatio-temporal Random Effects 
}

\description{
Implements a series of Poisson log-linear models with spatio-temporal random effects. The spatial and temporal autocorrelation in the random effects are induced using the class of Conditional AutoRegressive (CAR) priors, which are a special case of Gaussian Markov Random Fields (GMRF). All models are fitted in a Bayesian context using Markov chain Monte Carlo (McMC) simulation.
}

\details{
\tabular{ll}{
Package: \tab CARBayesST\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-11-13\cr
License: \tab GPL (>= 2)\cr
}


Implements Bayesian Poisson log-linear models for spatio-temporal areal unit count data. The models assume the response is Poisson, and model the mean as a combination of covariates, random effects and a piecewise constant cluster or intercept term (the latter is in ST.cluster and ST.ARCARcluster only). Six different models can be implemented: those proposed by Knorr-Held and Besag (1998), Knorr-Held (2000), Lee and Lawson (2014), Rushworth et al. (2014a) and Rushworth et al. (2014b).
}


\author{
Maintainer: Duncan Lee <Duncan.Lee@glasgow.ac.uk>
}


\references{
Knorr-Held, L. and J. Besag (1998). Modelling Risk from a Disease in Time and Space. Statistics in Medicine 17, 2045-2060.

Knorr-Held, L. (2000). Bayesian modelling of inseparable space-time variation in disease risk. Statistics in Medicine, 19, 2555-2567.

Lee, D. and A. Lawson (2014). Cluster detection and risk estimation for spatio-temporal health data. arXiv:1408.1191v2.

Rushworth, A., Lee, D., and Sarran, C.  An adaptive spatio-temporal smoothing model for estimating trends and step changes in disease risk (2014a) arXiv:1411.0924

Rushworth, A., D. Lee, and R. Mitchell (2014b). A spatio-temporal model for estimating the long-term effects of air pollution on respiratory hospital admissions in Greater London. Spatial and Spatio-temporal Epidemiology 10, 29-38.
}



\examples{
#### Artificial data generated on a square

#### Set up a square lattice region
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
n <- nrow(Grid)
t <- 10


#### set up distance and neighbourhood (W, based on sharing a common border) matrices
distance <-array(0, c(n,n))
W <-array(0, c(n,n))
     for(i in 1:n)
     {
		for(j in 1:n)
		{
		temp <- (Grid[i,1] - Grid[j,1])^2 + (Grid[i,2] - Grid[j,2])^2
		distance[i,j] <- sqrt(temp)
			if(temp==1)  W[i,j] <- 1 
		}	
	}
	
	
#### Generate data
n.all <- n * t
E <- rep(100, n.all)
log.risk <- log(rep(c(rep(1, 70), rep(2, 30)),t))
x <- rnorm(n.all)
risk <- exp(log.risk + 0.1 * x)
mean <- E * risk
Y <- rpois(n=n.all, lambda=mean)
formula <- Y~ offset(log(E)) + x
     

#### Run the models     
\dontrun{model1 <- ST.KHmain(formula, data=NULL, W=W, burnin=5000, n.sample=10000)}
\dontrun{model2 <- ST.KHinteraction(formula, data=NULL, W=W, burnin=5000, n.sample=10000)}
\dontrun{model3 <- ST.cluster(formula, data=NULL, G=3, K=100, burnin=5000, n.sample=10000)}
\dontrun{model4 <- ST.ARCAR(formula, data=NULL, W=W, burnin=5000, n.sample=10000)}
\dontrun{model5 <- ST.ARCARcluster(formula, data=NULL, G=3, W=W, burnin=5000, 
n.sample=10000)}
\dontrun{model6 <- ST.ARCARadaptive(formula, data=NULL, W=W, burnin=5000, 
n.sample=10000)}
\dontrun{print.CARBayesST(model6)}
}