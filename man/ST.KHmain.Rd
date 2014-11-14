\name{ST.KHmain}
\alias{ST.KHmain}


\title{
Fit the space-time Poisson log-linear model proposed by Knorr-Held and Besag (1998) with separable spatio-temporal main effects.
}

\description{
The function fits the Bayesian spatio-temporal model proposed by Knorr-Held and Besag (1998) to Poisson count data. The natural log of the linear predictor is made up of 5 terms, a covariate component, two spatial random effects, correlated and independent, and two temporal random effects, correlated and independent. The independent random effects are assigned zero-mean Gaussian priors with a common variance, while the correlated effects are modelled by intrinsic conditional  autoregressive (ICAR) and first order random walk (RW) models respectively. Inference is based on Markov chain Monte Carlo (McMC) simulation, using a combination of Gibbs sampling and Metropolis steps.
}


\usage{
ST.KHmain(formula, data=NULL, W, burnin=0, n.sample=1000, thin=1,  
prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, verbose=TRUE)
}


\arguments{
  \item{formula}{
A formula for the covariate part of the model, using the same notation as for the lm() function. The offsets should also be included here using the offset() function. The response and each covariate should be vectors of length (KN)*1, where K is the number of spatial units and N is the number of time periods. Each vector should be ordered so that the first K data points are the set of all K spatial locations at time 1, the next K are the set of spatial points for time 2 and so on.
}
\item{data}{
A data.frame containing the  variables in the formula.
}
 \item{W}{
A binary K by K neighbourhood matrix (where K is the number of spatial units). The jkth element equals one if areas (j, k) are spatially close (e.g. share a common border) and is zero otherwise. 
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
  \item{prior.mean.beta}{
A vector of prior means for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector of zeros.
}
  \item{prior.var.beta}{
A vector of prior variances for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector with values 1000.
}
  \item{prior.tau2}{
The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for each random effect variance term tau2. Defaults to c(0.001, 0.001).  
}
     \item{verbose}{
Logical, should the function update the user on its progress.  
}
}

\details{
For further details about how to apply the function see the examples below.
}


\value{
\item{formula }{The formula for the covariate and offset part of the model.}

\item{samples }{A list containing the MCMC samples from the model.}

\item{fitted.values }{A vector containing the fitted value for each area and time point. The vector is ordered so that all spatial units for time period one come first and then time period two and so on.}

\item{residuals }{A vector containing the residuals for each area and time point. The vector is ordered so that all spatial units for time period one come first and then time period two and so on.}

\item{stepchange }{This argument is NULL and is included for comparibility in the output object with the other models in this package.}

\item{modelfit }{Model fit criteria including the Deviance Information Criterion (DIC), the effective number of parameters in the model (p.d), and the log Marginal Predictive Likelihood (LMPL).}

\item{summary.results }{A table summarising some of the parameters in the model.}

\item{model }{A text string describing the model.}

\item{accept }{The acceptance probabilities for the parameters.}
}



\references{
Knorr-Held, L. and J. Besag (1998). Modelling Risk from a Disease in Time and Space. Statistics in Medicine 17, 2045-2060.
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
\dontrun{model1 <- ST.KHmain(formula, data=NULL, W=W, burnin=5000, n.sample=10000)}
}