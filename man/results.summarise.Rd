\name{results.summarise}
\alias{results.summarise}


\title{
Summarise a matrix of Markov chain Monte Carlo samples.
}

\description{
This function takes in a matrix of Markov chain Monte Carlo (McMC) samples from a set of parameters or quantities of interest, such as the fitted values, and calculates posterior quantiles and exceeedence probabilities. The latter are probabilities of the form P(quantity > c|data), where c is a threshold chosen by the user.
}


\usage{
results.summarise(samples, exceedences=NULL, quantiles=0.5)
}

\arguments{
  \item{samples}{
A matrix of McMC samples resulting from fitting one of the models. The object must be of class mcmc from the coda package.
}
  \item{exceedences}{
A vector of threshold levels, c, that you wish to calculate exceedence probabilities for.
}
 \item{quantiles}{
A vector of posterior quantiles required for McMC samples.
}
}

\details{
For further details about how to apply the function see the example below.
}


\value{
\item{exceedence }{A 2 dimensional array containing the requied exceedence probabilities. Each row relates to a parameter, and each column to a different requested exceedence probability. If the argument exceedences is missing this object is NULL.
}

\item{quantile }{A 2 dimensional array containing the requied posterior quantiles. Each row relates to a parameter, and each column to a different requested quantile. If the argument quantiles is missing this object is NULL.
}
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
     

#### Run the models     
\dontrun{model1 <- ST.knorrheld.main(formula, data=NULL, W=W, burnin=5000, n.sample=10000)}
\dontrun{samples <- model1$samples$fitted / E}
\dontrun{results <- results.summarise(samples, exceedences=c(1,1.5,2), 
quantiles=c(0.5, 0.025, 0.975))}
}