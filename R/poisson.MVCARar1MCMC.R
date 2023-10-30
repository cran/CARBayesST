poisson.MVCARar1MCMC <- function(Y, offset, X.standardised, W,  rho, alpha, fix.rho.S, fix.rho.T, K, N, NK, J, N.all, p, miss.locator, n.miss, burnin, n.sample, thin, MALA, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.Sigma.df, prior.Sigma.scale, verbose, chain)
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
  for(i in 1:J)
  {
  mod.glm <- glm(Y[ ,i]~X.standardised-1, offset=offset[ ,i], family="quasipoisson")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  beta[ ,i] <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
  }
  
log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
res.temp <- log.Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi.vec <- rnorm(n=N.all, mean=0, sd=res.sd)
phi <- matrix(phi.vec, ncol=J, byrow=TRUE)
Sigma <- cov(phi)
Sigma.inv <- solve(Sigma)
Sigma.a <- rep(1, J)
regression <- X.standardised %*% beta
fitted <- exp(regression + phi + offset)
Y.DA <- Y  


#### Matrices to store samples    
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, J*p))
samples.phi <- array(NA, c(n.keep, N.all))
samples.Sigma <- array(NA, c(n.keep, J, J))
samples.Sigma.a <- array(NA, c(n.keep, J))
    if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
    if(!fix.rho.T) samples.alpha <- array(NA, c(n.keep, 1))
samples.loglike <- array(NA, c(n.keep, N.all))
samples.fitted <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
  
  
#### Metropolis quantities
accept <- rep(0,4)
accept.beta <- rep(0,2*J)
proposal.sd.beta <- rep(0.01, J)
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.02
Sigma.post.df <- prior.Sigma.df + J - 1 + K * N  
Sigma.a.post.shape <- (prior.Sigma.df + J) / 2


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
  if(rho==1 & alpha==1) 
  {
  Sigma.post.df <- prior.Sigma.df + ((N-1) * (K-n.islands)) + J - 1
  }else if(rho==1)
  {
  Sigma.post.df <- prior.Sigma.df + (N * (K-n.islands)) + J - 1        
  }else if(alpha==1)
  {
  Sigma.post.df <- prior.Sigma.df + ((N-1) * K) + J - 1          
  }else
  {}


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
  #### Create the MCMC samples      
  for(j in 1:n.sample)
  {
    ####################################
    ## Sample from Y - data augmentation
    ####################################
    if(n.miss>0)
    {
      Y.DA[miss.locator] <- rpois(n=n.miss, lambda=fitted[miss.locator])    
    }else
    {}
    
    
    
    ###################
    ## Sample from beta
    ###################
    offset.temp <- phi + offset
    for(r in 1:J)
    {
      if(MALA)
      {
        temp <- poissonbetaupdateMALA(X.standardised, NK, p, beta[ ,r], offset.temp[ ,r], Y.DA[ ,r], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block)
      }else
      {
        temp <- poissonbetaupdateRW(X.standardised, NK, p, beta[ ,r], offset.temp[ ,r], Y.DA[ ,r], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block)
      }
      beta[ ,r] <- temp[[1]]
      accept.beta[r] <- accept.beta[r] + temp[[2]]
      accept.beta[(r+J)] <- accept.beta[(r+J)] + n.beta.block  
    }
    regression <- X.standardised %*% beta          
    
    
    
    ##################
    ## Sample from phi
    ##################
    #### Create the offset elements
    den.offset <- rho * W.triplet.sum + 1 - rho
    phi.offset <- regression + offset
    
    #### Create the random draws to create the proposal distribution
    Chol.Sigma <- t(chol(proposal.sd.phi*Sigma))
    z.mat <- matrix(rnorm(n=N.all, mean=0, sd=1), nrow=J, ncol=NK)
    innovations <- t(Chol.Sigma %*% z.mat)

    #### Update the elements of phi
    temp1 <- poissonmvar1carupdateRW(W.triplet, W.begfin, W.triplet.sum, K, N, J, phi, alpha, rho, Sigma.inv, Y.DA, innovations,  phi.offset, den.offset)      
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
    Sigma.post.scale <- 2 * prior.Sigma.df * diag(1 / Sigma.a) + t(phi[1:K, ]) %*% Q %*% phi[1:K, ]
      for(t in 2:N)
      {
      phit <- phi[((t-1)*K+1):(t*K), ]
      phitminus1 <- phi[((t-2)*K+1):((t-1)*K), ]
      temp1 <- phit - alpha * phitminus1
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
      temp  <- MVSTrhoTAR1compute(W.triplet, W.triplet.sum, n.triplet, den.offset, K, N, J, phi, rho, Sigma.inv)
      num <- temp[[1]]
      denom <- temp[[2]]
      alpha <- rnorm(n=1, mean = (num / denom), sd=sqrt(1 / denom))
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
      temp1.QF  <- MVSTrhoSAR1compute(W.triplet, W.triplet.sum, n.triplet, den.offset, K, N, J, phi, rho, alpha, Sigma.inv)
      temp2.QF  <- MVSTrhoSAR1compute(W.triplet, W.triplet.sum, n.triplet, proposal.den.offset, K, N, J, phi, proposal.rho, alpha, Sigma.inv)        
 
 
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
    fitted <- exp(regression + phi + offset)
    loglike <- dpois(x=as.numeric(t(Y)), lambda=as.numeric(t(fitted)), log=TRUE)
    
    
    
    ###################
    ## Save the results
    ###################
    if(j > burnin & (j-burnin)%%thin==0)
    {
      ele <- (j - burnin) / thin
      samples.beta[ele, ] <- as.numeric(beta)
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
      for(r in 1:J)
      {
        if(p>2)
        {
          proposal.sd.beta[r] <- common.accceptrates1(accept.beta[c(r, (r+J))], proposal.sd.beta[r], 40, 50)
        }else
        {
          proposal.sd.beta[r] <- common.accceptrates1(accept.beta[c(r, (r+J))], proposal.sd.beta[r], 30, 40)    
        }
      }
      
      proposal.sd.phi <- common.accceptrates1(accept[1:2], proposal.sd.phi, 40, 50)
        if(!fix.rho.S)
        {
        proposal.sd.rho <- common.accceptrates2(accept[3:4], proposal.sd.rho, 40, 50, 0.5)
        }
      accept <- c(0,0,0,0)
      accept.beta <- rep(0,2*J)
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
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.Sigma=samples.Sigma, samples.Sigma.a=samples.Sigma.a,  samples.rho=samples.rho, samples.alpha=samples.alpha, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept, accept.beta=accept.beta)

  
#### Return the results
return(chain.results)
}