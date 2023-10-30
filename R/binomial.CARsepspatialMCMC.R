binomial.CARsepspatialMCMC <- function(Y, failures, trials, offset, X.standardised, W, rho, lambda, fix.rho.S, fix.rho.T, K, N, N.all, p, burnin, n.sample, thin, MALA, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.tau2, verbose, chain)
{
#Rcpp::sourceCpp("src/CARBayesST.cpp")   
#source("R/common.functions.R")
#library(spdep)
#library(truncnorm)    
#     
#     
############################################
#### Set up the key elements before sampling
############################################
#### Generate the initial parameter values
dat <- cbind(Y, failures)
mod.glm <- glm(dat~X.standardised-1, offset=offset, family="quasibinomial")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
  
theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=N.all, mean=0, sd = res.sd)
phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)  
delta <- rnorm(n=N, mean=0, sd = res.sd)
tau2 <- apply(phi.mat, 2, var) / 10
sig2 <- var(delta)/10

  
#### Matrix versions
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE)
trials.mat <- matrix(trials, nrow=K, ncol=N, byrow=FALSE)
failures.mat <- matrix(failures, nrow=K, ncol=N, byrow=FALSE)
delta.mat <- matrix(delta, nrow=K, ncol=N, byrow=TRUE)


#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, N.all))
samples.tau2 <- array(NA, c(n.keep, N))
samples.sig2 <- array(NA, c(n.keep, 1))
    if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
    if(!fix.rho.T) samples.lambda <- array(NA, c(n.keep, 1))
samples.delta <- array(NA, c(n.keep, N))     
samples.fitted <- array(NA, c(n.keep, N.all))
samples.loglike <- array(NA, c(n.keep, N.all))


#### Specify the Metropolis quantities
accept <- rep(0,10)
proposal.sd.phi <- 0.1
proposal.sd.rho <- 0.05
proposal.sd.beta <- 0.01
proposal.sd.delta <- 0.05
proposal.sd.lambda <- 0.02
tau2.shape <- prior.tau2[1] + K/2
sig2.shape <- prior.tau2[1] + N/2


#### CAR quantities
W.quants <- common.Wcheckformat.leroux(W)
K <- W.quants$n
N <- N.all / K
W <- W.quants$W
W.triplet <- W.quants$W.triplet
W.n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin


#### Spatial determinant
    if(!fix.rho.S) 
    {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))     
    }else
    {}

  
#### .T quantities
D <-array(0, c(N,N))
    for(i in 1:N)
    {
        for(j in 1:N)
        {
        if(abs((i-j))==1)  D[i,j] <- 1 
        }    
    }

D.triplet <- c(NA, NA, NA)
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            if(D[i,j]>0)
            {
            D.triplet <- rbind(D.triplet, c(i,j, D[i,j]))     
            }else{}
        }
    }
D.triplet <- D.triplet[-1, ]     
D.n.triplet <- nrow(D.triplet) 
D.triplet.sum <- tapply(D.triplet[ ,3], D.triplet[ ,1], sum)
D.neighbours <- tapply(D.triplet[ ,3], D.triplet[ ,1], length)
 
D.begfin <- array(NA, c(N, 2))     
temp <- 1
    for(i in 1:N)
    {
    D.begfin[i, ] <- c(temp, (temp + D.neighbours[i]-1))
    temp <- temp + D.neighbours[i]
    }

    if(!fix.rho.T) 
    {
    Dstar <- diag(apply(D,1,sum)) - D
    Dstar.eigen <- eigen(Dstar)
    Dstar.val <- Dstar.eigen$values
    det.Q.D <-  0.5 * sum(log((lambda * Dstar.val + (1-lambda))))    
    }else
    {} 


#### Check for islands
W.list<- mat2listw(W, style = "B")
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
n.island1 <- length(which(islands==1))
    if(rho==1) tau2.shape <- prior.tau2[1] + 0.5 * (K-n.islands)   
    if(lambda==1) sig2.shape <- prior.tau2[1] + 0.5 * (N-1)     
  

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
    ###################
    ## Sample from beta
    ###################
    offset.temp <- as.numeric(offset.mat + phi.mat + delta.mat)    
    if(MALA)
    {
      temp <- binomialbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y, failures, trials, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
    }else
    {
      temp <- binomialbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y, failures, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
    }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)  

    ####################
    ## Sample from phi
    ####################
    phi.offset <- offset.mat + regression.mat + delta.mat
    den.offset <- rho * W.triplet.sum + 1 - rho
    temp1 <- binomialsrecarupdateRW(W.triplet, W.begfin, W.triplet.sum, K, N, phi.mat, rho, Y.mat, failures.mat, proposal.sd.phi, phi.offset, den.offset, tau2)
    phi.temp <- temp1[[1]]
    phi.mean <- apply(phi.temp,2,mean)
    if(rho<1)
    {
      phi <- as.numeric(phi.temp) - kronecker(phi.mean, rep(1,K))
    }else
    {
      phi.temp[which(islands==1), ] <- phi.temp[which(islands==1), ] - matrix(kronecker(phi.mean, rep(1,n.island1)), ncol=N, byrow=F) 
      phi <- as.numeric(phi.temp)
    }
    phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
    
    accept[3] <- accept[3] + temp1[[2]]
    accept[4] <- accept[4] + N.all
    
    #####################
    ## Samples from delta
    #####################
    delta.offset <- t(offset.mat + phi.mat + regression.mat)
    temp2 <- binomialcarupdateRW(D.triplet, D.begfin, D.triplet.sum, N, delta, sig2, t(Y.mat), t(failures.mat), proposal.sd.delta, lambda, delta.offset, K, rep(1,K))
    delta <- temp2[[1]]
    delta <- delta - mean(delta)
    delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)
    accept[7] <- accept[7] + temp2[[2]]
    accept[8] <- accept[8] + N      

    ####################
    ## Samples from tau2
    ####################
    tau2.temp <- tauquadformcompute2(W.triplet, W.triplet.sum, W.n.triplet, K, N, phi.mat, rho)
    tau2 <- tau2compute(tau2, tau2.temp, tau2.shape, prior.tau2[2], N)
    
    ####################
    ## Samples from sig2
    ####################
    temp2.delta <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, lambda)
    sig2.scale <- temp2.delta + prior.tau2[2] 
    sig2 <- 1 / rgamma(1, sig2.shape, scale=(1/sig2.scale))
    
    ##################
    ## Sample from rho
    ##################
    if(!fix.rho.S)
    {
      temp3 <- rhoquadformcompute(W.triplet, W.triplet.sum, W.n.triplet, K, N, phi.mat, rho, tau2)
      proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)
      temp4 <- rhoquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, proposal.rho, tau2)
      det.Q.W.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))
      logprob.current <- N * det.Q.W - temp3
      logprob.proposal <- N * det.Q.W.proposal - temp4
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
    
    #####################
    ## Sample from lambda
    #####################
    if(!fix.rho.T)
    {
      proposal.lambda <- rtruncnorm(n=1, a=0, b=1, mean=lambda, sd=proposal.sd.lambda)   
      temp3 <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, proposal.lambda)
      det.Q.proposal <- 0.5 * sum(log((proposal.lambda * Dstar.val + (1-proposal.lambda))))              
      logprob.current <- det.Q.D - temp2.delta / sig2
      logprob.proposal <- det.Q.proposal - temp3 / sig2
      hastings <- log(dtruncnorm(x=lambda, a=0, b=1, mean=proposal.lambda, sd=proposal.sd.lambda)) - log(dtruncnorm(x=proposal.lambda, a=0, b=1, mean=lambda, sd=proposal.sd.lambda)) 
      prob <- exp(logprob.proposal - logprob.current + hastings)
      
      #### Accept or reject the proposal
      if(prob > runif(1))
      {
        lambda <- proposal.lambda
        det.Q.D <- det.Q.proposal
        accept[9] <- accept[9] + 1           
      }else
      {
      }              
      accept[10] <- accept[10] + 1           
    }else
    {}
    
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(offset.mat + regression.mat + phi.mat + delta.mat)
    prob <- exp(lp) / (1+exp(lp))
    fitted <- trials * prob
    loglike <- dbinom(x=Y, size=trials, prob=prob, log=TRUE)


    
    ###################
    ## Save the results
    ###################
    if(j > burnin & (j-burnin)%%thin==0)
    {
      ele <- (j - burnin) / thin
      samples.beta[ele, ] <- beta
      samples.phi[ele, ] <- as.numeric(phi)
      if(!fix.rho.S) samples.rho[ele, ] <- rho
      if(!fix.rho.T) samples.lambda[ele, ] <- lambda
      samples.tau2[ele, ] <- tau2
      samples.sig2[ele, ] <- sig2
      samples.delta[ele, ] <- delta
      samples.fitted[ele, ] <- fitted
      samples.loglike[ele, ] <- loglike
    }else
    {
    }
    
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
      proposal.sd.delta <- common.accceptrates1(accept[7:8], proposal.sd.delta, 40, 50)
      if(!fix.rho.S) proposal.sd.rho <- common.accceptrates2(accept[5:6], proposal.sd.rho, 40, 50, 0.5)
      if(!fix.rho.T) proposal.sd.lambda <- common.accceptrates2(accept[9:10], proposal.sd.lambda, 40, 50, 0.5)
      accept <- rep(0,10)
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
    if(fix.rho.S) samples.rho <- NA
    if(fix.rho.T) samples.lambda <- NA
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.delta=samples.delta, samples.lambda=samples.lambda, samples.tau2=samples.tau2, samples.rho=samples.rho, samples.sig2=samples.sig2, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    accept=accept)

         
#### Return the results
return(chain.results)
}