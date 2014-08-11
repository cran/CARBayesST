\name{ST.clustconv}
\alias{ST.clustconv}


\title{
Fit the cluster conv Poisson log-linear model proposed by Lee and Lawson (2014).
}

\description{
The function fits the Bayesian spatio-temporal clustering model with a spatial convolution of Gaussian noise proposed by Lee and Lawson (2014) to Poisson count data. The natural log of the linear predictor is made up of a clustering component and a spatial convolution of zero-mean Gaussian noise with a constant variance, the former of which is piecewise constant with G risk levels or classes. The mean in each risk level changes over time, and each area is allowed to change risk class over time. The mean risk levels and the allocation of an area to a risk class are both temporally correlated. The convolution component is based on a spatial Gaussian kernel function. Inference is based on Markov chain Monte Carlo (McMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}


\usage{
ST.clustconv(formula, data=NULL, coordinates=NULL, G, burnin=0, n.sample=1000, 
thin=1,  blocksize.beta=5, prior.mean.beta=NULL, prior.var.beta=NULL, 
prior.tau2=NULL, prior.sigma2=NULL, prior.alpha=NULL, prior.rho=NULL, 
verbose=TRUE)
}

\arguments{
  \item{formula}{
A formula for the covariate part of the model, using the same notation as for the lm() function. The offsets should also be included here using the offset() function. The response and each covariate should be vectors of length (KN)*1, where K is the number of spatial units and N is the number of time periods. Each vector should be ordered so that the first K data points are the set of all K spatial locations at time 1, the next K are the set of spatial points for time 2 and so on.
}
\item{data}{
A data.frame containing the  variables in the formula.
}
 \item{coordinates}{
A K by 2 matrix giving the coordinates centroids of the K areal units.
}
 \item{G}{
The maximum number of risk classes. Note, some classes can be estimated as being empty.
}
  \item{burnin}{
The number of MCMC samples to discard as the burnin period. Defaults to 0.
}
  \item{n.sample}{
The number of MCMC samples to generate. Defaults to 1,000.
}
  \item{thin}{
The level of thinning to apply to the MCMC samples to reduce their temporal autocorrelation. Defaults to 1.
}
  \item{blocksize.beta}{
The size of the blocks in which to update the regression parameters in the MCMC algorithm. Defaults to 5.
}
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector with values 1000.
}
  \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for the variance of the random Gaussian noise. Defaults to c(0.001, 0.001).  
}
  \item{prior.sigma2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for the random walk prior for the risk levels. Defaults to c(0.001, 0.001).  
}
  \item{prior.alpha}{
The prior maximum for a Unifrom(0, prior.alpha) prior for the temporal correlation model on the allocation variables of an area to a risk class. Defaults to 10.  
}
  \item{prior.rho}{
The prior maximum for a Unifrom(0, prior.rho) prior for the bandwidth of the Gaussian kernel function. Defaults to 10.  
}
     \item{verbose}{
Logical, should the function update the user on its progress.  
}
}

\details{
For further details about how to apply the function see the examples below and in the vignette.
}


\value{
\item{formula }{The formula for the covariate and offset part of the model.}

\item{samples }{A list containing the MCMC samples from the model.}

\item{fitted.values }{A vector containing the fitted value for each area and time point. The vector is ordered so that all spatial units at time one come first and then at time two and so on.}

\item{residuals }{A vector containing the residuals for each area and time point. The vector is ordered so that all spatial units at time one come first and then at time two and so on.}

\item{posterior.Z }{A matrix containing the discrete posterior distribution for each allocation variable. Each row corresponds to one variable, and the columns are the posterior probabilities of being in each of the G risk classes. The matrix is ordered so that all spatial units at time one come first and then at time two and so on.}

\item{median.Z }{A K*N matrix giving the posterior median of the allocation variable assigning a data point to a risk class. Each row corresponds to an areal unit and each column to a time period.}

\item{modelfit }{Model fit criteria including the Deviance Information Criterion (DIC), the effective number of parameters in the model (p.d), DIC3, Watanabe-Akaike Information Criterion (WAIC), and the log Marginal Predictive Likelihood (LMPL).}

\item{summary.results }{A table summarising some of the parameters in the model.}

\item{model }{A text string describing the model.}

\item{accept }{The acceptance probabilities for the parameters.}
}



\references{
Lee, D. and A. Lawson (2014). Cluster detection and risk estimation for spatio-temporal health data. arXiv:1408.1191.
}

\author{
Duncan Lee
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
     

#### Run the model
\dontrun{model5 <- ST.clustconv(formula, data=NULL, G=3, coordinates=Grid, burnin=5000, 
n.sample=10000)}
}