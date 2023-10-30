binomial.MVCARar2MCMC <- function(Y, failures, trials, offset, X.standardised, W,  rho, alpha, fix.rho.S, fix.rho.T, K, N, NK, J, N.all, p, miss.locator, n.miss, burnin, n.sample, thin, MALA, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.Sigma.df, prior.Sigma.scale, verbose, chain)
{
#cpp::sourceCpp("src/CARBayesST.cpp")   
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
  mod.glm <- glm(cbind(Y[ ,i], failures[ ,i])~X.standardised-1, offset=offset[ ,i], family="quasibinomial")
  beta.mean <- mod.glm$coefficients
  beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
  beta[ ,i] <- rnorm(n=p, mean=beta.mean, sd=beta.sd)
  }

theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - X.standardised %*% beta - offset
phi <- res.temp
phi[is.na(phi)] <- rnorm(n=sum(is.na(phi)), mean=0, sd=sd(res.temp, na.rm=T))
Sigma <- cov(phi)
Sigma.inv <- solve(Sigma)
Sigma.a <- rep(1, J)
regression <- X.standardised %*% beta
lp <- regression + phi + offset
prob <- exp(lp)  / (1 + exp(lp))
Y.DA <- Y  
failures.DA <- failures


#### Matrices to store samples    
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, J*p))
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
    Y.DA[miss.locator] <- rbinom(n=n.miss, size=trials[miss.locator], prob=prob[miss.locator])    
    failures.DA <- trials - Y.DA
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
        temp <- binomialbetaupdateMALA(X.standardised, NK, p, beta[ ,r], offset.temp[ ,r], Y.DA[ ,r], failures.DA[ ,r], trials[ ,r], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block)
      }else
      {
        temp <- binomialbetaupdateRW(X.standardised, NK, p, beta[ ,r], offset.temp[ ,r], Y.DA[ ,r], failures.DA[ ,r], prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta[r], list.block)
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
  temp1 <- binomialmvar2carupdateRW(W.triplet, W.begfin, W.triplet.sum, K, N, J, phi, alpha[1], alpha[2], rho, Sigma.inv, Y.DA, failures.DA, innovations, phi.offset, den.offset)      
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
    lp <- regression + phi + offset
    prob <- exp(lp)  / (1 + exp(lp))
    fitted <- trials * prob
    loglike <- dbinom(x=as.numeric(t(Y)), size=as.numeric(t(trials)), prob=as.numeric(t(prob)), log=TRUE)

    
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