gaussian.MVCARar2MCMC <- function(Y, offset, X.standardised, W,  rho, alpha, fix.rho.S, fix.rho.T, K, N, NK, J, N.all, p, miss.locator, n.miss, burnin, n.sample, thin, prior.mean.beta, prior.var.beta, prior.nu2, prior.Sigma.df, prior.Sigma.scale, verbose, chain)
{
#Rcpp::sourceCpp("src/CARBayesST.cpp")   
#source("R/common.functions.R")
#library(spdep)
#library(truncnorm)
#library(MCMCpack)
#     
#     
############################################
#### Set up the key elements before sampling
############################################
#### Generate the initial parameter values
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
Sigma.a <- rep(1, J)
regression <- X.standardised %*% beta
fitted <- regression + phi + offset
Y.DA <- Y  


#### Matrices to store samples    
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, J*p))
samples.nu2 <- array(NA, c(n.keep, J))
samples.phi <- array(NA, c(n.keep, N.all))
samples.Sigma <- array(NA, c(n.keep, J, J))
samples.Sigma.a <- array(NA, c(n.keep, J))
  if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
  if(!fix.rho.T) samples.alpha <- array(NA, c(n.keep, 2))
samples.loglike <- array(NA, c(n.keep, N.all))
samples.fitted <- array(NA, c(n.keep, N.all))
  if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
  
  
#### Metropolis quantities
accept <- rep(0,4)
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.02
Sigma.post.df <- prior.Sigma.df + J - 1 + K * N  
Sigma.a.post.shape <- (prior.Sigma.df + J) / 2
nu2.posterior.shape <- prior.nu2[1] + 0.5 * K * N


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
W.list<- mat2listw(W, style = "B")
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
    if(rho==1 & alpha[1]==2 & alpha[2]==-1) 
    {
    Sigma.post.df <- prior.Sigma.df + ((N-2) * (K-n.islands)) + J - 1   
    }else if(rho==1)
    {
    Sigma.post.df <- prior.Sigma.df + (N * (K-n.islands)) + J - 1        
    }else if(alpha[1]==2 & alpha[2]==-1)
    {
    Sigma.post.df <- prior.Sigma.df + ((N-2) * K) + J - 1             
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
  

#### Start timer
    if(verbose)
    {
    cat("\nMarkov chain", chain,  "- generating", n.keep, "post burnin and thinned samples.\n", sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage.points<-round((1:100/100)*n.sample)
    }else
    {
    percentage.points<-round((1:100/100)*n.sample)     
    }



##############################
#### Generate the MCMC samples
##############################
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
    Sigma.post.scale <- 2 * prior.Sigma.df * diag(1 / Sigma.a) + t(phi[1:K, ]) %*% Q %*% phi[1:K, ] + t(phi[(K+1):(2*K), ]) %*% Q %*% phi[(K+1):(2*K), ]
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
    ## Sample from Sigma.a
    ######################
    Sigma.a.posterior.scale <- prior.Sigma.df * diag(Sigma.inv) + 1 / prior.Sigma.scale^2
    Sigma.a <- 1 / rgamma(J, Sigma.a.post.shape, scale=(1/Sigma.a.posterior.scale))  

    
    
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
    samples.Sigma.a[ele, ] <- Sigma.a
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
      if(ceiling(j/100)==floor(j/100) & j < burnin)
      {
      #### Update the proposal sds
      proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
        if(!fix.rho.S)
        {
        proposal.sd.rho <- common.accceptrates2(accept[3:4], proposal.sd.rho, 40, 50, 0.5)
        }
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



############################################
#### Return the results to the main function
############################################
#### Compile the results
    if(n.miss==0) samples.Y <- NA
    if(fix.rho.S) samples.rho <- NA
    if(fix.rho.T) samples.alpha <- NA
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.nu2=samples.nu2, samples.Sigma=samples.Sigma, samples.Sigma.a=samples.Sigma.a,  samples.rho=samples.rho, samples.alpha=samples.alpha, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept)

  
#### Return the results
return(chain.results)
}