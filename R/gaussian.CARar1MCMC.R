gaussian.CARar1MCMC <- function(Y, offset, X.standardised, W,  rho, gamma, fix.rho.S, fix.rho.T, K, N, N.all, p, which.miss, n.miss, burnin, n.sample, thin, prior.mean.beta, prior.var.beta, prior.tau2, prior.nu2, verbose, chain)
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
mod.glm <- glm(Y~X.standardised-1, offset=offset)
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)

res.temp <- Y - X.standardised %*% beta - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=N.all, mean=0, sd = res.sd)
tau2 <- var(phi)/10
nu2 <- runif(1, 0, res.sd)


#### Matrix versions of quantites
Y.DA <- Y  
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)   
fitted <- as.numeric(offset.mat + regression.mat + phi.mat)


#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, N.all))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.nu2 <- array(NA, c(n.keep, 1))
    if(!fix.rho.S) samples.rho <- array(NA, c(n.keep, 1))
    if(!fix.rho.T) samples.gamma <- array(NA, c(n.keep, 1))
samples.fitted <- array(NA, c(n.keep, N.all))
samples.loglike <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

    
#### Specify the Metropolis quantities
accept <- rep(0,2)
proposal.sd.rho <- 0.05
tau2.shape <- prior.tau2[1] + N.all/2
nu2.shape <- prior.nu2[1] + N.all/2  


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
  if(rho==1 & gamma==1) 
  {
  tau2.shape <- prior.tau2[1] + prior.tau2[1] + ((N-1) * (K-n.islands))/2
  }else if(rho==1)
  {
  tau2.shape <- prior.tau2[1] + prior.tau2[1] + (N * (K-n.islands))/2        
  }else if(gamma==1)
  {
  tau2.shape <- prior.tau2[1] + prior.tau2[1] + ((N-1) * K)/2          
  }else
  {}


#### Beta update quantities
data.precision.beta <- t(X.standardised) %*% X.standardised
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
        Y.DA[which.miss==0] <- rnorm(n=n.miss, mean=fitted[which.miss==0], sd=sqrt(nu2))    
        }else
        {}
    Y.DA.mat <- matrix(Y.DA, nrow=K, ncol=N, byrow=FALSE)
        
        
        
    ##################
    ## Sample from nu2
    ##################
    nu2.offset <- as.numeric(Y.DA.mat - offset.mat - regression.mat - phi.mat)
    nu2.scale <- prior.nu2[2]  + sum(nu2.offset^2)/2
    nu2 <- 1 / rgamma(1, nu2.shape, scale=(1/nu2.scale)) 

    
        
    ####################
    ## Sample from beta
    ####################
    fc.precision <- prior.precision.beta + data.precision.beta / nu2
    fc.var <- solve(fc.precision)
    beta.offset <- as.numeric(Y.DA.mat - offset.mat - phi.mat)
    beta.offset2 <- t(X.standardised) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
    fc.mean <- fc.var %*% beta.offset2
    chol.var <- t(chol(fc.var))
    beta <- fc.mean + chol.var %*% rnorm(p)        
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)  

        
        
    ####################
    ## Sample from phi
    ####################
    phi.offset <- Y.DA.mat - offset.mat - regression.mat
    den.offset <- rho * W.triplet.sum + 1 - rho
    phi.temp <- gaussianar1carupdate(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, nu2, gamma, rho, phi.offset, den.offset)      
    phi <- as.numeric(phi.temp)  - mean(as.numeric(phi.temp))
    phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)

        
        
    ####################
    ## Sample from gamma
    ####################
        if(!fix.rho.T)
        {
        temp2 <- gammaquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho)
        mean.gamma <- temp2[[1]] / temp2[[2]]
        sd.gamma <- sqrt(tau2 / temp2[[2]])
        gamma <- rtruncnorm(n=1, a=0, b=1, mean=mean.gamma, sd=sd.gamma)
        }else
        {}
        
        
        
    ####################
    ## Samples from tau2
    ####################
    temp3 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, rho, gamma)
    tau2.scale <- temp3 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.shape, scale=(1/tau2.scale)) 
        
        
        
    ##################
    ## Sample from rho
    ##################
        if(!fix.rho.S)
        {
        proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)
        temp4 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, proposal.rho, gamma)
        det.Q.W.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))
        logprob.current <- N * det.Q.W - temp3 / tau2
        logprob.proposal <- N * det.Q.W.proposal - temp4 / tau2
        hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
        prob <- exp(logprob.proposal - logprob.current + hastings)
            if(prob > runif(1))
            {
            rho <- proposal.rho
            det.Q.W <- det.Q.W.proposal
            accept[1] <- accept[1] + 1           
            }else
            {}              
        accept[2] <- accept[2] + 1       
        }else
        {}
        
    
        
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(offset.mat + regression.mat + phi.mat)
    loglike <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),N.all), log=TRUE)

        
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- as.numeric(phi)
            if(!fix.rho.S) samples.rho[ele, ] <- rho
            if(!fix.rho.T) samples.gamma[ele, ] <- gamma
        samples.tau2[ele, ] <- tau2
        samples.nu2[ele, ] <- nu2
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
            if(!fix.rho.S) proposal.sd.rho <- common.accceptrates2(accept[1:2], proposal.sd.rho, 40, 50, 0.5)
        accept <- rep(0,2)
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
    if(fix.rho.T) samples.gamma <- NA
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.tau2=samples.tau2, samples.nu2=samples.nu2, samples.rho=samples.rho, samples.gamma=samples.gamma, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept)


#### Return the results
return(chain.results)
}