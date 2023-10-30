gaussian.CARlinearMCMC <- function(Y, offset, X.standardised, W,  rho, lambda, fix.rho.int, fix.rho.slo, K, N, N.all, p, which.miss, n.miss, burnin, n.sample, thin, prior.mean.beta, prior.var.beta, prior.mean.alpha, prior.var.alpha, prior.tau2, prior.nu2, verbose, chain)
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
time <-(1:N - mean(1:N))/N
time.all <- kronecker(time, rep(1,K))
mod.glm <- glm(Y~X.standardised-1 + time.all, offset=offset)
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
temp <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
beta <- temp[1:p]
alpha <- temp[(p+1)]

res.temp <- Y - as.numeric(X.standardised %*% beta) - time.all * alpha - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=0, sd = res.sd)
delta <- rnorm(n=K, mean=0, sd = res.sd)
tau2.phi <- var(phi)/10
tau2.delta <- var(delta)/10
nu2 <- runif(1, 0, res.sd)


#### Specify matrix quantities
Y.DA <- Y  
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)
time.mat <- matrix(rep(time, K), byrow=TRUE, nrow=K)    
delta.time.mat <- apply(time.mat, 2, "*", delta)
alpha.offset1 <- sum(time.mat^2)
fitted <- as.numeric(offset.mat + regression.mat + phi.mat + delta.time.mat + alpha * time.mat)


#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.alpha <- array(NA, c(n.keep, 1))
samples.phi <- array(NA, c(n.keep, K))
samples.delta <- array(NA, c(n.keep, K))
    if(!fix.rho.int) samples.rho <- array(NA, c(n.keep, 1))
    if(!fix.rho.slo) samples.lambda <- array(NA, c(n.keep, 1))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 2))
colnames(samples.tau2) <- c("tau2.int", "tau2.slo")
samples.fitted <- array(NA, c(n.keep, N.all))
samples.loglike <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

    
#### Specify the Metropolis quantities
accept <- rep(0,4)
proposal.sd.rho <- 0.02
proposal.sd.lambda <- 0.02
nu2.shape <- prior.nu2[1] + N*K/2    
tau2.phi.shape <- prior.tau2[1] + K/2
tau2.delta.shape <- prior.tau2[1] + K/2


#### CAR quantities
W.quants <- common.Wcheckformat.leroux(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
W.n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin


#### Create the determinant     
    if(!fix.rho.int | !fix.rho.slo) 
    {
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    }else
    {}
    if(!fix.rho.int) det.Q.rho <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))    
    if(!fix.rho.slo) det.Q.lambda <-  0.5 * sum(log((lambda * Wstar.val + (1-lambda))))     


#### Check for islands
W.list<- mat2listw(W, style = "B")
W.nb <- W.list$neighbours
W.islands <- n.comp.nb(W.nb)
islands <- W.islands$comp.id
n.islands <- max(W.islands$nc)
    if(rho==1) tau2.phi.shape <- prior.tau2[1] + 0.5 * (K-n.islands)   
    if(lambda==1) tau2.delta.shape <- prior.tau2[1] + 0.5 * (K-n.islands)     


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
    nu2.offset <- as.numeric(Y.DA.mat - offset.mat - regression.mat - phi.mat - delta.time.mat - alpha * time.mat)
    nu2.scale <- prior.nu2[2]  + sum(nu2.offset^2)/2
    nu2 <- 1 / rgamma(1, nu2.shape, scale=(1/nu2.scale)) 

    
    
    ####################
    ## Sample from beta
    ####################
    fc.precision <- prior.precision.beta + data.precision.beta / nu2
    fc.var <- solve(fc.precision)
    beta.offset <- as.numeric(Y.DA.mat - offset.mat - phi.mat - delta.time.mat - alpha * time.mat)
    beta.offset2 <- t(X.standardised) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
    fc.mean <- fc.var %*% beta.offset2
    chol.var <- t(chol(fc.var))
    beta <- fc.mean + chol.var %*% rnorm(p)        
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)  
    

        
    ####################
    ## Sample from alpha
    ####################
    fc.var <- 1 / (1 / prior.var.alpha + alpha.offset1 / nu2)
    alpha.offset <- (Y.DA.mat - offset.mat - regression.mat - phi.mat - delta.time.mat) * time.mat
    alpha.offset2 <- sum(alpha.offset, na.rm=TRUE) / nu2
    fc.mean <- fc.var * (alpha.offset2 +  prior.mean.alpha / prior.var.alpha)
    alpha <- rnorm(n=1, mean=fc.mean, sd=sqrt(fc.var))
         
      
    
    ####################
    ## Sample from phi
    ####################
    phi.offset <- Y.DA.mat - offset.mat - regression.mat - delta.time.mat - alpha * time.mat
    phi.offset2 <- apply(phi.offset,1, sum, na.rm=TRUE)
    temp1 <- gaussiancarupdate(W.triplet, W.begfin, W.triplet.sum, K, phi, tau2.phi, nu2, phi.offset2, rho, N)
    phi <- temp1
        if(rho<1)
        {
        phi <- phi - mean(phi)
        }else
        {
        phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
    phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)    

        
        
    ####################
    ## Sample from delta
    ####################
    delta.offset <- (Y.DA.mat - offset.mat - regression.mat - phi.mat - alpha * time.mat) * time.mat
    delta.offset2 <- apply(delta.offset,1, sum, na.rm=TRUE)
    temp2 <- gaussiancarupdate(W.triplet, W.begfin, W.triplet.sum, K, delta, tau2.delta, nu2, delta.offset2, lambda, sum(time^2))
    delta <- temp2
    if(lambda <1)
    {
    delta <- delta - mean(delta)
    }else
    {
    delta[which(islands==1)] <- delta[which(islands==1)] - mean(delta[which(islands==1)])   
    }
    delta.time.mat <- apply(time.mat, 2, "*", delta)
        
        
        
    #######################
    ## Sample from tau2.phi
    #######################
    temp2.phi <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, rho)
    tau2.phi.scale <- temp2.phi + prior.tau2[2] 
    tau2.phi <- 1 / rgamma(1, tau2.phi.shape, scale=(1/tau2.phi.scale))
        
        
    #########################
    ## Sample from tau2.delta
    #########################
    temp2.delta <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, delta, delta, lambda)
    tau2.delta.scale <- temp2.delta + prior.tau2[2] 
    tau2.delta <- 1 / rgamma(1, tau2.delta.shape, scale=(1/tau2.delta.scale))
        
        
        
    ##################
    ## Sample from rho
    ##################
        if(!fix.rho.int)
        {
        proposal.rho <- rtruncnorm(n=1, a=0, b=1, mean=rho, sd=proposal.sd.rho)   
        temp3 <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, proposal.rho)
        det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
        logprob.current <- det.Q.rho - temp2.phi / tau2.phi
        logprob.proposal <- det.Q.proposal - temp3 / tau2.phi
        hastings <- log(dtruncnorm(x=rho, a=0, b=1, mean=proposal.rho, sd=proposal.sd.rho)) - log(dtruncnorm(x=proposal.rho, a=0, b=1, mean=rho, sd=proposal.sd.rho)) 
        prob <- exp(logprob.proposal - logprob.current + hastings)
        
        #### Accept or reject the proposal
            if(prob > runif(1))
            {
            rho <- proposal.rho
            det.Q.rho <- det.Q.proposal
            accept[1] <- accept[1] + 1           
            }else
            {}              
        accept[2] <- accept[2] + 1           
        }else
        {}
        
     
       
    #####################
    ## Sample from lambda
    #####################
        if(!fix.rho.slo)
        {
        proposal.lambda <- rtruncnorm(n=1, a=0, b=1, mean=lambda, sd=proposal.sd.lambda)   
        temp3 <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, delta, delta, proposal.lambda)
        det.Q.proposal <- 0.5 * sum(log((proposal.lambda * Wstar.val + (1-proposal.lambda))))              
        logprob.current <- det.Q.lambda - temp2.delta / tau2.delta
        logprob.proposal <- det.Q.proposal - temp3 / tau2.delta
        hastings <- log(dtruncnorm(x=lambda, a=0, b=1, mean=proposal.lambda, sd=proposal.sd.lambda)) - log(dtruncnorm(x=proposal.lambda, a=0, b=1, mean=lambda, sd=proposal.sd.lambda)) 
        prob <- exp(logprob.proposal - logprob.current + hastings)
        
        #### Accept or reject the proposal
            if(prob > runif(1))
            {
            lambda <- proposal.lambda
            det.Q.lambda <- det.Q.proposal
            accept[3] <- accept[3] + 1           
            }else
            {}              
        accept[4] <- accept[4] + 1           
        }else
        {}
        
     
       
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(offset.mat + regression.mat + phi.mat + delta.time.mat + alpha * time.mat)
    loglike <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),N.all), log=TRUE)
        
    
    
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.beta[ele, ] <- beta
        samples.phi[ele, ] <- phi
        samples.delta[ele, ] <- delta
        samples.alpha[ele, ] <- alpha
            if(!fix.rho.int) samples.rho[ele, ] <- rho
            if(!fix.rho.slo) samples.lambda[ele, ] <- lambda
        samples.nu2[ele, ] <- nu2
        samples.tau2[ele, ] <- c(tau2.phi, tau2.delta)
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
            if(!fix.rho.int) proposal.sd.rho <- common.accceptrates2(accept[1:2], proposal.sd.rho, 40, 50, 0.5)
            if(!fix.rho.slo) proposal.sd.lambda <- common.accceptrates2(accept[3:4], proposal.sd.lambda, 40, 50, 0.5)
        accept <- rep(0,4)            
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
    if(fix.rho.int) samples.rho <- NA
    if(fix.rho.slo) samples.lambda <- NA
chain.results <- list(samples.beta=samples.beta, samples.alpha=samples.alpha, samples.phi=samples.phi, samples.delta=samples.delta, samples.tau2=samples.tau2, samples.nu2=samples.nu2, samples.rho=samples.rho, samples.lambda=samples.lambda, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept)


#### Return the results
return(chain.results)
}