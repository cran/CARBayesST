gaussian.MVCARar2 <- function(formula, data=NULL,  W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, rho.S=NULL, rho.T=NULL, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
  
  
  
#### Frame object
frame.results <- common.frame.MVST(formula, data, "gaussian")
NK <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
N.all <- length(Y)
J <- ncol(Y)
which.miss <- frame.results$which.miss
n.miss <- N.all - sum(which.miss)
Y.DA <- Y   
  
  
#### Create a missing list
  if(n.miss>0)
  {
  miss.locator <- array(NA, c(n.miss, 2))
  colnames(miss.locator) <- c("row", "column")
  locations <- which(which.miss==0)
  miss.locator[ ,1] <- ceiling(locations/J)
  miss.locator[ ,2] <- locations - (miss.locator[ ,1]-1) * J
  }else
  {}
  
  
#### W matrix
  if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
K <- nrow(W)
N <- NK / K
  if(ceiling(N)!= floor(N)) stop("The number of data points in Y divided by the number of rows in W is not a whole number.", call.=FALSE)
  
  
#### Check on the rho arguments
  if(is.null(rho.S))
  {
  rho <- runif(1)
  fix.rho.S <- FALSE   
  }else
  {
  rho <- rho.S
  fix.rho.S <- TRUE
  }
  if(!is.numeric(rho)) stop("rho.S is fixed but is not numeric.", call.=FALSE)  
  if(rho<0 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)  
  if(rho>1 ) stop("rho.S is outside the range [0, 1].", call.=FALSE)    

  if(is.null(rho.T))
  {
  alpha <- c(runif(1), runif(1))
  fix.rho.T <- FALSE   
  }else
  {
  alpha <- rho.T
  fix.rho.T <- TRUE
  }
  if(!is.numeric(alpha)) stop("rho.T is fixed but is not numeric.", call.=FALSE)  
  if(length(alpha)!=2) stop("rho.T is fixed but is not of length 2.", call.=FALSE)  


#### Priors
  if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
  if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
  if(is.null(prior.Sigma.df)) prior.Sigma.df <- J+1
  if(is.null(prior.Sigma.scale)) prior.Sigma.scale <- diag(rep(1/1000,J))
  if(is.null(prior.nu2)) prior.nu2 <- c(1, 0.01)
prior.beta.check(prior.mean.beta, prior.var.beta, p)
common.prior.varmat.check(prior.Sigma.scale, J)  
prior.var.check(prior.nu2)     
  
  
#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)  
  
  
  
#############################
#### Initial parameter values
#############################
beta <- array(NA, c(p, J))
nu2 <- rep(NA, J)
  for(i in 1:J)
  {
  mod.glm <- lm(Y[ ,i]~X.standardised-1, offset=offset[ ,i])
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.unscaled)) * summary(mod.glm)$sigma
  beta[ ,i] <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
  nu2[i] <- runif(1, var(mod.glm$residuals)*0.5, var(mod.glm$residuals))
  }
  
res.temp <- Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi.vec <- rnorm(n=N.all, mean=0, sd=res.sd)
phi <- matrix(phi.vec, ncol=J, byrow=TRUE)
Sigma <- cov(phi)
Sigma.inv <- solve(Sigma)
regression <- X.standardised %*% beta
fitted <- regression + phi + offset
  
  
  
###############################    
#### Set up the MCMC quantities    
###############################
#### Matrices to store samples    
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, J*p))
samples.nu2 <- array(NA, c(n.keep, J))
samples.phi <- array(NA, c(n.keep, N.all))
samples.Sigma <- array(NA, c(n.keep, J, J))
  if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
  if(!fix.rho.T) samples.alpha <- array(NA, c(n.keep, 2))
samples.loglike <- array(NA, c(n.keep, N.all))
samples.fitted <- array(NA, c(n.keep, N.all))
  if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
  
  
#### Metropolis quantities
accept <- rep(0,4)
accept.all <- rep(0,4)
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.02
Sigma.post.df <- prior.Sigma.df + K * N  
nu2.posterior.shape <- prior.nu2[1] + 0.5 * K * N
  
  
  
##################################
#### Set up the spatial quantities
##################################
#### CAR quantities
W.quants <- common.Wcheckformat.leroux(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin
Wstar <- diag(apply(W,1,sum)) - W
Q <- rho * Wstar + diag(rep(1-rho,K))
  
  
#### Create the determinant     
  if(!fix.rho.S)
  {
  Wstar.eigen <- eigen(Wstar)
  Wstar.val <- Wstar.eigen$values
  det.Q <- sum(log((rho * Wstar.val + (1-rho))))    
  }else
  {} 
  
  
#### Check for islands
W.list<- mat2listw(W)
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
  if(rho==1 & alpha[1]==2 & alpha[2]==-1) 
  {
  Sigma.post.df <- prior.Sigma.df + ((N-2) * (K-n.islands))/2
  }else if(rho==1)
  {
  Sigma.post.df <- prior.Sigma.df + (N * (K-n.islands))/2        
  }else if(alpha[1]==2 & alpha[2]==-1)
  {
  Sigma.post.df <- prior.Sigma.df + ((N-2) * K)/2          
  }else
  {}


#### Beta update quantities
data.precision <- t(X.standardised) %*% X.standardised 
  if(length(prior.var.beta)==1)
  {
  prior.precision.beta <- 1 / prior.var.beta
  }else
  {
  prior.precision.beta <- solve(diag(prior.var.beta))
  }
  
  
  
  ###########################
  #### Run the Bayesian model
  ###########################
  #### Start timer
  if(verbose)
  {
  cat("Generating", n.keep, "post burnin and thinned (if requested) samples.\n", sep = " ")
  progressBar <- txtProgressBar(style = 3)
  percentage.points<-round((1:100/100)*n.sample)
  }else
  {
  percentage.points<-round((1:100/100)*n.sample)     
  }
  
  
  #### Create the MCMC samples      
  for(j in 1:n.sample)
  {
  ####################################
  ## Sample from Y - data augmentation
  ####################################
    if(n.miss>0)
    {
    Y.DA[miss.locator] <- rnorm(n=n.miss, mean=fitted[miss.locator], sd=sqrt(nu2[miss.locator[ ,2]]))    
    }else
    {}
    
    
    
  ##################
  ## Sample from nu2
  ##################
  fitted.current <- regression + phi + offset
  nu2.posterior.scale <- prior.nu2[2] + 0.5 * apply((Y.DA - fitted.current)^2, 2, sum)
  nu2 <- 1 / rgamma(J, nu2.posterior.shape, scale=(1/nu2.posterior.scale))
    
    
    
  ###################
  ## Sample from beta
  ###################
    for(r in 1:J)
    {
    fc.precision <- prior.precision.beta + data.precision / nu2[r]
    fc.var <- solve(fc.precision)    
    fc.temp1 <- t(((Y.DA[, r] - phi[ , r] - offset[ , r]) %*%  X.standardised) / nu2[r]) + prior.precision.beta %*% prior.mean.beta  
    fc.mean <- fc.var %*% fc.temp1
    chol.var <- t(chol(fc.var))
    beta[ ,r] <- fc.mean + chol.var %*% rnorm(p)   
    }
  regression <- X.standardised %*% beta    
    
    
    
  ##################
  ## Sample from phi
  ##################
  #### Create the offset elements
  den.offset <- rho * W.triplet.sum + 1 - rho
  phi.offset <- Y.DA - regression - offset
    
  #### Create the random draws to create the proposal distribution
  Chol.Sigma <- t(chol(proposal.sd.phi*Sigma))
  z.mat <- matrix(rnorm(n=N.all, mean=0, sd=1), nrow=J, ncol=NK)
  innovations <- t(Chol.Sigma %*% z.mat)
    
  #### Update the elements of phi
  temp1 <- gaussianmvar2carupdateRW(W.triplet, W.begfin, W.triplet.sum, K, N, J, phi, alpha[1], alpha[2],  rho, Sigma.inv, nu2, innovations,  phi.offset, den.offset)      
  phi <- temp1[[1]]
    for(r in 1:J)
    {
    phi[ ,r] <- phi[ ,r] - mean(phi[ ,r])    
    }
  accept[1] <- accept[1] + temp1[[2]]
  accept[2] <- accept[2] + NK
    
    
  
  ####################
  ## Sample from Sigma
  ####################
  Sigma.post.scale <- prior.Sigma.scale + t(phi[1:K, ]) %*% Q %*% phi[1:K, ] + t(phi[(K+1):(2*K), ]) %*% Q %*% phi[(K+1):(2*K), ]
    for(t in 3:N)
    {
    phit <- phi[((t-1)*K+1):(t*K), ]
    phitminus1 <- phi[((t-2)*K+1):((t-1)*K), ]
    phitminus2 <- phi[((t-3)*K+1):((t-2)*K), ]
    temp1 <- phit - alpha[1] * phitminus1 - alpha[2] * phitminus2
    Sigma.post.scale <- Sigma.post.scale + t(temp1) %*%  Q %*% temp1
    }
  Sigma <- riwish(Sigma.post.df, Sigma.post.scale)
  Sigma.inv <- solve(Sigma)

    
    
  ######################
  #### Sample from alpha
  ######################
    if(!fix.rho.T)
    {
    temp  <- MVSTrhoTAR2compute(W.triplet, W.triplet.sum, n.triplet, den.offset, K, N, J, phi, rho, Sigma.inv)
    alpha.precision <- matrix(c(temp[[1]], temp[[2]], temp[[2]], temp[[3]]), nrow=2, ncol=2)
    alpha.var <- solve(alpha.precision)
    alpha.mean <- rep(NA, 2)
    alpha.mean[2] <- (temp[[1]] * temp[[5]] - temp[[2]] * temp[[4]]) / (temp[[1]] * temp[[3]] - temp[[2]]^2)
    alpha.mean[1] <-  (temp[[5]] - temp[[3]] * alpha.mean[2]) / temp[[2]]
    alpha <- mvrnorm(n=1, mu=alpha.mean, Sigma=alpha.var)
    }else
    {}
  
  

  ##################
  ## Sample from rho
  ##################
    if(!fix.rho.S)
    {
    ## Propose a new value
    proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)
    proposal.Q <- proposal.rho * Wstar + diag(rep(1-proposal.rho), K)
    proposal.det.Q <-  sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))  
    proposal.den.offset <- proposal.rho * W.triplet.sum + 1 - proposal.rho
    
    ## Compute the quadratic forms based on current and proposed values of rho
    temp1.QF  <- MVSTrhoSAR2compute(W.triplet, W.triplet.sum, n.triplet, den.offset, K, N, J, phi, rho, alpha[1], alpha[2], Sigma.inv)
    temp2.QF  <- MVSTrhoSAR2compute(W.triplet, W.triplet.sum, n.triplet, proposal.den.offset, K, N, J, phi, proposal.rho, alpha[1], alpha[2], Sigma.inv)        
    
    ## Compute the acceptance rate
    logprob.current <- 0.5 * J * N * det.Q - 0.5 * temp1.QF
    logprob.proposal <- 0.5 * J * N * proposal.det.Q - 0.5 * temp2.QF
    hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
    prob <- exp(logprob.proposal - logprob.current + hastings)
      if(prob > runif(1))
      {
      rho <- proposal.rho
      det.Q <- proposal.det.Q
      Q <- proposal.Q
      accept[3] <- accept[3] + 1           
      }else
      {}              
    accept[4] <- accept[4] + 1       
    }else
    {}
  

    
  #########################
  ## Calculate the deviance
  #########################
  fitted <- regression + phi + offset
  loglike <- dnorm(x=as.numeric(t(Y)), mean=as.numeric(t(fitted)), sd=rep(sqrt(nu2), K*N), log=TRUE)
    
  
  
  ###################
  ## Save the results
  ###################
    if(j > burnin & (j-burnin)%%thin==0)
    {
    ele <- (j - burnin) / thin
    samples.beta[ele, ] <- as.numeric(beta)
    samples.nu2[ele, ] <- nu2
    samples.phi[ele, ] <- as.numeric(t(phi))
    samples.Sigma[ele, , ] <- Sigma
      if(!fix.rho.S) samples.rho[ele, ] <- rho
      if(!fix.rho.T) samples.alpha[ele, ] <- alpha
    samples.loglike[ele, ] <- loglike
    samples.fitted[ele, ] <- as.numeric(t(fitted))
      if(n.miss>0) samples.Y[ele, ] <- Y.DA[miss.locator]
    }else
    {}
    
    
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    k <- j/100
      if(ceiling(k)==floor(k))
      {
      #### Update the proposal sds
      proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
        if(!fix.rho.S)
        {
        proposal.sd.rho <- common.accceptrates2(accept[3:4], proposal.sd.rho, 40, 50, 0.5)
        }
      accept.all <- accept.all + accept
      accept <- c(0,0,0,0)
      }else
      {}
    
    
  ################################       
  ## print progress to the console
  ################################
    if(j %in% percentage.points & verbose)
    {
    setTxtProgressBar(progressBar, j/n.sample)
    }
  }
  
  
  ##### end timer
  if(verbose)
  {
  cat("\nSummarising results.")
  close(progressBar)
  }else
  {}
  

  
###################################
#### Summarise and save the results 
###################################
#### Compute the acceptance rates
accept.beta <- 100 
accept.phi <- 100 * accept.all[1] / accept.all[2]
  if(!fix.rho.S)
  {
  accept.rho <- 100 * accept.all[3] / accept.all[4]
  }else
  {
  accept.rho <- NA    
  }
  accept.Sigma <- 100
  if(!fix.rho.T)
  {
  accept.alpha <- 100
  }else
  {
  accept.alpha <- NA
  }
accept.final <- c(accept.beta, accept.phi, accept.rho, accept.Sigma, accept.alpha)
names(accept.final) <- c("beta", "phi", "rho.S", "Sigma", "rho.T")
  
  
#### Compute the fitted deviance
mean.beta <- matrix(apply(samples.beta, 2, mean), nrow=p, ncol=J, byrow=F)
mean.phi <- matrix(apply(samples.phi, 2, mean), nrow=NK, ncol=J, byrow=T)
fitted.mean <- X.standardised %*% mean.beta + mean.phi + offset
nu2.mean <- apply(samples.nu2,2,mean)
deviance.fitted <- -2 * sum(dnorm(as.numeric(t(Y)), mean = as.numeric(t(fitted.mean)), sd=rep(sqrt(nu2.mean), K*N), log = TRUE), na.rm=TRUE)

  
#### Model fit criteria
modelfit <- common.modelfit(samples.loglike, deviance.fitted)
  
  
#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- samples.beta
  for(r in 1:J)
  {
  samples.beta.orig[ ,((r-1)*p+1):(r*p)] <- common.betatransform(samples.beta[ ,((r-1)*p+1):(r*p) ], X.indicator, X.mean, X.sd, p, FALSE)
  }
  
  
#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
col.name <- rep(NA, p*(J-1))
  
  if(is.null(colnames(Y)))
  {
    for(r in 1:J)
    {
    col.name[((r-1)*p+1):(r*p)] <- paste("Variable ", r,  " - ", colnames(X), sep="")   
    }
  }else
  {
    for(r in 1:J)
    {
    col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[r],  " - ", colnames(X), sep="")   
    }
  }
rownames(summary.beta) <- col.name
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
  
summary.hyper <- array(NA, c((2*J+3) ,7))
summary.hyper[1:J, 1:3] <-t(apply(samples.nu2, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.hyper[1:J, 4] <- rep(n.keep, J)
summary.hyper[1:J, 5] <- rep(100, J)
summary.hyper[1:J, 6] <- apply(samples.nu2, 2, effectiveSize)
summary.hyper[1:J, 7] <- geweke.diag(samples.nu2)$z
  
summary.hyper[(J+1):(2*J), 1] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.5)))
summary.hyper[(J+1):(2*J), 2] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.025)))
summary.hyper[(J+1):(2*J), 3] <- diag(apply(samples.Sigma, c(2,3), quantile, c(0.975)))
summary.hyper[(J+1):(2*J), 4] <- rep(n.keep, J)
summary.hyper[(J+1):(2*J), 5] <- rep(100, J)
summary.hyper[(J+1):(2*J), 6] <- diag(apply(samples.Sigma, c(2,3), effectiveSize))
  for(r in 1:J)
  {
  summary.hyper[J+r, 7] <- geweke.diag(samples.Sigma[ ,r,r])$z    
  }
  
  if(!fix.rho.S)
  {
  summary.hyper[(2*J+1), 1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
  summary.hyper[(2*J+1), 4:5] <- c(n.keep, accept.rho)
  summary.hyper[(2*J+1), 6:7] <- c(effectiveSize(samples.rho), geweke.diag(samples.rho)$z)
  }else
  {
  summary.hyper[(2*J+1), 1:3] <- c(rho, rho, rho)
  summary.hyper[(2*J+1), 4:5] <- rep(NA, 2)
  summary.hyper[(2*J+1), 6:7] <- rep(NA, 2)
  }
  
  if(!fix.rho.T)
  {
  summary.hyper[(2*J+2), 1:3] <- quantile(samples.alpha[ ,1], c(0.5, 0.025, 0.975))
  summary.hyper[(2*J+2), 4:5] <- c(n.keep, accept.alpha)
  summary.hyper[(2*J+2), 6:7] <- c(effectiveSize(samples.alpha[ ,1]), geweke.diag(samples.alpha[ ,1])$z)
  summary.hyper[(2*J+3), 1:3] <- quantile(samples.alpha[ ,2], c(0.5, 0.025, 0.975))
  summary.hyper[(2*J+3), 4:5] <- c(n.keep, accept.alpha)
  summary.hyper[(2*J+3), 6:7] <- c(effectiveSize(samples.alpha[ ,2]), geweke.diag(samples.alpha[ ,2])$z)
  }else
  {
  summary.hyper[(2*J+2), 1:3] <- c(alpha[1], alpha[1], alpha[1])
  summary.hyper[(2*J+2), 4:5] <- rep(NA, 2)
  summary.hyper[(2*J+2), 6:7] <- rep(NA, 2)
  summary.hyper[(2*J+3), 1:3] <- c(alpha[2], alpha[2], alpha[2])
  summary.hyper[(2*J+3), 4:5] <- rep(NA, 2)
  summary.hyper[(2*J+3), 6:7] <- rep(NA, 2)
  }
  
summary.results <- rbind(summary.beta, summary.hyper)
rownames(summary.results)[((J*p)+1): nrow(summary.results)] <- c(paste(rep("nu2",J), 1:J, sep=""), paste(rep("Sigma",J), 1:J, 1:J, sep=""), "rho.S", "rho1.T", "rho2.T")
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
  
  
#### Create the fitted values and residuals
fitted.values <- matrix(apply(samples.fitted, 2, mean), nrow=NK, ncol=J, byrow=T)
response.residuals <- Y - fitted.values
nu.mat <- matrix(rep(sqrt(nu2.mean), N*K), nrow=N*K, byrow=T)
pearson.residuals <- response.residuals / nu.mat
residuals <- list(response=response.residuals, pearson=pearson.residuals)
  
  
#### Compile and return the results
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nRandom effects model - Multivariate Autoregressive order 2 CAR model\n")


#### Harmonise samples in case of them not being generated
  if(fix.rho.S & fix.rho.T)
  {
  samples.rhoext <- NA
  }else if(fix.rho.S & !fix.rho.T)
  {
  samples.rhoext <- samples.alpha
  colnames(samples.rhoext) <- c("rho1.T", "rho2.T")
  }else if(!fix.rho.S & fix.rho.T)
  {
  samples.rhoext <- samples.rho  
  names(samples.rhoext) <- "rho.S"
  }else
  {
  samples.rhoext <- cbind(samples.rho, samples.alpha)
  colnames(samples.rhoext) <- c("rho.S", "rho1.T", "rho2.T")
  }

  if(n.miss==0) samples.Y = NA
samples <- list(beta=samples.beta.orig, phi=mcmc(samples.phi), Sigma=samples.Sigma, nu2=mcmc(samples.nu2), rho=mcmc(samples.rhoext), fitted=mcmc(samples.fitted), Y=mcmc(samples.Y))
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL,  formula=formula, model=model.string, X=X)
class(results) <- "CARBayesST"
  
  
#### Finish by stating the time taken    
  if(verbose)
  {
  b<-proc.time()
  cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
  }else
  {}
return(results)
}
