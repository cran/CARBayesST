poisson.CARar2MCMC <- function(Y, offset, X.standardised, W,  rho, alpha, fix.rho.S, fix.rho.T, K, N, N.all, p, which.miss, n.miss, burnin, n.sample, thin, MALA, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.tau2, verbose, chain)
{
#Rcpp::sourceCpp("src/CARBayesST.cpp")   
#source("R/common.functions.R")
#library(spdep)
#library(truncnorm)  
#library(MASS)
#     
#     
############################################
#### Set up the key elements before sampling
############################################
#### Generate the initial parameter values
mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
log.Y <- log(Y)
log.Y[Y==0] <- -0.1  
res.temp <- log.Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=N.all, mean=0, sd = res.sd)
tau2 <- var(phi)/10
  
  
#### Specify matrix quantities
Y.DA <- Y  
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)   
fitted <- exp(as.numeric(offset.mat + regression.mat + phi.mat))


#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, N.all))
samples.tau2 <- array(NA, c(n.keep, 1))
  if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
  if(!fix.rho.T) samples.alpha <- array(NA, c(n.keep, 2))
samples.fitted <- array(NA, c(n.keep, N.all))
samples.loglike <- array(NA, c(n.keep, N.all))
  if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))
  
  
#### Specify the Metropolis quantities
accept <- rep(0,6)
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.05
proposal.sd.beta <- 0.01
tau2.shape <- prior.tau2[1] + N.all/2


#### CAR quantities
W.quants <- common.Wcheckformat.leroux(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
W.n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin


#### Create the determinant     
    if(!fix.rho.S) 
    {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
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
  tau2.shape <- prior.tau2[1] + prior.tau2[1] + ((N-2) * (K-n.islands))/2
  }else if(rho==1)
  {
  tau2.shape <- prior.tau2[1] + prior.tau2[1] + (N * (K-n.islands))/2        
  }else if(alpha[1]==2 & alpha[2]==-1)
  {
  tau2.shape <- prior.tau2[1] + prior.tau2[1] + ((N-2) * K)/2          
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
    Y.DA[which.miss==0] <- rpois(n=n.miss, lambda=fitted[which.miss==0])    
    }else
    {}
    Y.DA.mat <- matrix(Y.DA, nrow=K, ncol=N, byrow=FALSE)
    
    
    
  ####################
  ## Sample from beta
  ####################
  offset.temp <- as.numeric(offset.mat + phi.mat)     
    if(MALA)
    {
    temp <- poissonbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
    }else
    {
    temp <- poissonbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
    }
  beta <- temp[[1]]
  accept[1] <- accept[1] + temp[[2]]
  accept[2] <- accept[2] + n.beta.block  
  regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)  
    
    
    
  ####################
  ## Sample from phi
  ####################
  phi.offset <- offset.mat + regression.mat
  den.offset <- rho * W.triplet.sum + 1 - rho
  temp1 <- poissonar2carupdateRW(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, alpha[1], alpha[2], rho, Y.DA.mat, proposal.sd.phi, phi.offset, den.offset)      
  phi.temp <- temp1[[1]]
  phi <- as.numeric(phi.temp)  - mean(as.numeric(phi.temp))
  phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
  accept[3] <- accept[3] + temp1[[2]]
  accept[4] <- accept[4] + K*N
    
    
    ####################
    ## Sample from alpha
    ####################
    if(!fix.rho.T)
    {
    #### Construct the quadratic forms
    temp2 <- alphaquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho, tau2)

    #### Construct the precision matrix
    alpha.prec <- array(c(temp2[[1]], temp2[[3]], temp2[[3]], temp2[[2]]), c(2,2))
    alpha.var <- solve(alpha.prec)
    
    #### Construct the mean vector
    U2 <- (temp2[[1]] * temp2[[5]] - temp2[[3]] * temp2[[4]]) / (temp2[[2]] * temp2[[1]] - temp2[[3]]^2) 
    U1 <- (1 / temp2[[3]]) * (temp2[[5]] - temp2[[2]] * U2)
    alpha.mean <- c(U1, U2)
    alpha <- mvrnorm(n=1, mu=alpha.mean, Sigma=alpha.var)
    }else
    {}
  
  
  
    ####################
    ## Samples from tau2
    ####################
    temp3 <- tauquadformcomputear2(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho, alpha[1], alpha[2])
    tau2.scale <- temp3 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.shape, scale=(1/tau2.scale))  
    
    

    
    ##################
    ## Sample from rho
    ##################
    if(!fix.rho.S)
    {
      proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)
      temp4 <- tauquadformcomputear2(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, proposal.rho, alpha[1], alpha[2])
      det.Q.W.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))
      logprob.current <- N * det.Q.W - temp3 / tau2
      logprob.proposal <- N * det.Q.W.proposal - temp4 / tau2
      hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
      prob <- exp(logprob.proposal - logprob.current + hastings)
      
      if(prob > runif(1))
      {
        rho <- proposal.rho
        det.Q.W <- det.Q.W.proposal
        accept[5] <- accept[5] + 1           
      }else
      {
      }              
      accept[6] <- accept[6] + 1       
    }else
    {}
    
    
    
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(offset.mat + regression.mat + phi.mat)
    fitted <- exp(lp)
    loglike <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)
    
    
    ###################
    ## Save the results
    ###################
    if(j > burnin & (j-burnin)%%thin==0)
    {
      ele <- (j - burnin) / thin
      samples.beta[ele, ] <- beta
      samples.phi[ele, ] <- as.numeric(phi)
      samples.tau2[ele, ] <- tau2
      if(!fix.rho.S) samples.rho[ele, ] <- rho
      if(!fix.rho.T) samples.alpha[ele, ] <- alpha
      samples.fitted[ele, ] <- fitted
      samples.loglike[ele, ] <- loglike
      if(n.miss>0) samples.Y[ele, ] <- Y.DA[which.miss==0]
    }else
    {}
    
    
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
    if(ceiling(j/100)==floor(j/100) & j < burnin)
    {
      #### Update the proposal sds
      if(p>2)
      {
        proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 40, 50)
      }else
      {
        proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 30, 40)    
      }
      proposal.sd.phi <- common.accceptrates1(accept[3:4], proposal.sd.phi, 40, 50)
      if(!fix.rho.S) proposal.sd.rho <- common.accceptrates2(accept[5:6], proposal.sd.rho, 40, 50, 0.5) 
      accept <- rep(0,6)  
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
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.tau2=samples.tau2, samples.rho=samples.rho, samples.alpha=samples.alpha, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept)


#### Return the results
return(chain.results)
}