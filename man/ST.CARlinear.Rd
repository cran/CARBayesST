\name{ST.CARlinear}
\alias{ST.CARlinear}


\title{
    Fit a spatio-temporal generalised linear mixed model to data, where the spatial
    units have linear time trends with spatially varying intercepts and slopes.
}

\description{
    Fit a spatio-temporal generalised linear mixed model to areal unit data, where the response
    variable can be binomial, Gaussian or Poisson. The linear predictor is modelled by known 
    covariates and an area specific linear time trend. The intercepts and slopes are spatially
    autocorrelated and modelled by the conditional autoregressive (CAR) prior proposed by Leroux
    et al. (1999). The model is similar to that proposed by Bernardinelli et al. (1995) and 
    further details are given in the vignette accompanying this package. Inference is conducted 
    in a Bayesian setting using Markov chain Monte Carlo (McMC) simulation. 
}

\usage{
    ST.CARlinear(formula, family, data=NULL,  trials=NULL, W, burnin, n.sample, thin=1, 
    prior.mean.beta=NULL, prior.var.beta=NULL, prior.mean.alpha=NULL, 
    prior.var.alpha=NULL, prior.nu2=NULL, prior.tau2=NULL, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{formula}{
        A formula for the covariate part of the model using the syntax of the
        lm() function. Offsets can be included here using the offset() function. 
        The response and each covariate should be vectors of length (KT)*1, where
        K is the number of spatial units and T is the number of time periods. Each
        vector should be ordered so that the first K data points are the set of all
        K spatial locations at time 1, the next K are the set of spatial points for
        time 2 and so on.
    }
    \item{family}{
        One of either `binomial', `gaussian' or `poisson', which respectively specify a 
        binomial likelihood model with a logistic link function, a Gaussian likelihood 
        model with an identity link function, or a Poisson likelihood model with a 
        log link function. 
    }
    \item{data}{
        An optional data.frame containing the  variables in the formula.
    }
    \item{trials}{
        A vector the same length as the response containing the total number of trials 
        for each area and time period. Only used if family=`binomial'. 
    }
    \item{W}{
        A K by K neighbourhood matrix (where K is the number of spatial units). 
        Typically a binary specification is used, where the jkth element equals one
        if areas (j, k) are spatially close (e.g. share a common border) and is zero 
        otherwise. The matrix can be non-binary, but each row must contain at least 
        one non-zero entry.
    }   
    \item{burnin}{
        The number of McMC samples to discard as the burnin period.
    }
    \item{n.sample}{
        The number of McMC samples to generate.
    }
    \item{thin}{
        The level of thinning to apply to the McMC samples to reduce their temporal 
        autocorrelation. Defaults to 1.
    }
    \item{prior.mean.beta}{
        A vector of prior means for the regression parameters beta (Gaussian priors are 
        assumed). Defaults to a vector of zeros.
    }
    \item{prior.var.beta}{
        A vector of prior variances for the regression parameters beta (Gaussian priors 
        are assumed). Defaults to a vector with values 1000.
    } 
    \item{prior.mean.alpha}{
        The prior mean for the average slope of the linear time trend alpha (Gaussian priors are 
        assumed). Defaults to  zero.
    }
    \item{prior.var.alpha}{
        The prior variance for the average slope of the linear time trend alpha (Gaussian priors 
        are assumed). Defaults to 1000.
    }  
    \item{prior.nu2}{
        The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
        prior for the Gaussian error variance nu2. Defaults to c(0.001, 0.001) and only used if
        family=`Gaussian'.   
    }
    \item{prior.tau2}{
        The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) 
        prior for the random effect variances tau2. Defaults to c(0.001, 0.001).  
    }
    \item{verbose}{
        Logical, should the function update the user on its progress.  
    }
        }


\value{
    \item{summary.results }{A summary table of the parameters.}
    \item{samples }{A list containing the McMC samples from the model. The `tau2' element
                    of this list has columns (tau2.phi, tau2.delta. Similarly, the `rho' 
                    element of this list has columns (rho.phi, rho.delta).}
    \item{fitted.values }{A vector of fitted values for each area and time period.}
    \item{residuals }{A vector of residuals for each area and time period.}
    \item{modelfit }{Model fit criteria including the Deviance Information Criterion 
                     (DIC), the effective number of parameters in the model (p.d), and the Log 
                    Marginal Predictive Likelihood (LMPL).}
    \item{accept }{The acceptance probabilities for the parameters.}
    \item{localised.structure }{NULL, for compatability with the other models.}
    \item{formula }{The formula for the covariate and offset part of the model.}
    \item{model }{A text string describing the model fit.}
    \item{X }{The design matrix of covariates.}
}

\references{
Bernardinelli, L., D. Clayton, C.Pascuto, C.Montomoli, M.Ghislandi, and M. Songini
(1995). Bayesian analysis of space-time variation in disease risk. Statistics in 
Medicine, 14, 2433-2443.

Leroux, B., X. Lei, and N. Breslow (1999). Estimation of disease rates in small 
areas: A new mixed model for spatial dependence, Chapter Statistical Models in 
Epidemiology, the Environment and Clinical Trials, Halloran, M and Berry, D (eds),
pp. 135-178. Springer-Verlag, New York.
}

\author{
    Duncan Lee
}




\examples{
##################################################
#### Run the model on simulated data on a lattice
##################################################
#### set up the regular lattice    
x.easting <- 1:10
x.northing <- 1:10
Grid <- expand.grid(x.easting, x.northing)
K <- nrow(Grid)
T <- 10
N.all <- K * T

#### set up spatial neighbourhood matrix W
W <-array(0, c(K,K))
     for(i in 1:K)
     {
        for(j in 1:K)
		{
		temp <- (Grid[i,1] - Grid[j,1])^2 + (Grid[i,2] - Grid[j,2])^2
			if(temp==1)  W[i,j] <- 1 
		}	
	}


#### Simulate the elements in the linear predictor and the data
x <- rnorm(n=N.all, mean=0, sd=1)
beta <- 0.1
Q.W <- 0.99 * (diag(apply(W, 2, sum)) - W) + 0.01 * diag(rep(1,K))
Q.W.inv <- solve(Q.W)
phi <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.1 * Q.W.inv))
delta <- mvrnorm(n=1, mu=rep(0,K), Sigma=(0.1 * Q.W.inv))
trend <- array(NA, c(K, T))
time <-(1:T - mean(1:T))/T
    for(i in 1:K)
    {
    trend[i, ] <- phi[i] + delta[i] * time        
    }
trend.vec <- as.numeric(trend)
LP <- 4 + x * beta + trend.vec
mean <- exp(LP)
Y <- rpois(n=N.all, lambda=mean)


#### Run the model
\dontrun{model <- ST.CARlinear(formula=Y~x, family="poisson", W=W, burnin=10000, 
n.sample=50000)}
}
