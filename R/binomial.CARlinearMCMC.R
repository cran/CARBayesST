binomial.CARlinearMCMC <- function(Y, failures, trials, offset, X.standardised, W,  rho, lambda, fix.rho.int, fix.rho.slo, K, N, N.all, p, which.miss, n.miss, burnin, n.sample, thin, MALA, n.beta.block, list.block, prior.mean.beta, prior.var.beta, prior.mean.alpha, prior.var.alpha, prior.tau2, verbose, chain)
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
dat <- cbind(Y, failures)
mod.glm <- glm(dat~X.standardised-1 + time.all, offset=offset, family="quasibinomial")
beta.mean <- mod.glm$coefficients
beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
temp <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
beta <- temp[1:p]
alpha <- temp[(p+1)]

theta.hat <- Y / trials
theta.hat[theta.hat==0] <- 0.01
theta.hat[theta.hat==1] <- 0.99
res.temp <- log(theta.hat / (1 - theta.hat)) - as.numeric(X.standardised %*% beta) - time.all * alpha - offset
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi <- rnorm(n=K, mean=0, sd = res.sd)
delta <- rnorm(n=K, mean=0, sd = res.sd)
tau2.phi <- var(phi)/10
tau2.delta <- var(delta)/10

#### Matrix versions
Y.DA <- Y  
failures.DA <- failures
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
trials.mat <- matrix(trials, nrow=K, ncol=N, byrow=FALSE)
phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)
time.mat <- matrix(rep(time, K), byrow=TRUE, nrow=K)    
delta.time.mat <- apply(time.mat, 2, "*", delta)
lp <- as.numeric(offset.mat + regression.mat + phi.mat + delta.time.mat + alpha * time.mat)
prob <- exp(lp)  / (1 + exp(lp))


#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.alpha <- array(NA, c(n.keep, 1))
samples.phi <- array(NA, c(n.keep, K))
samples.delta <- array(NA, c(n.keep, K))
    if(!fix.rho.int) samples.rho <- array(NA, c(n.keep, 1))
    if(!fix.rho.slo) samples.lambda <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 2))
colnames(samples.tau2) <- c("tau2.int", "tau2.slo")
samples.fitted <- array(NA, c(n.keep, N.all))
samples.loglike <- array(NA, c(n.keep, N.all))
    if(n.miss>0) samples.Y <- array(NA, c(n.keep, n.miss))

    
#### Specify the Metropolis quantities
accept <- rep(0,12)
proposal.sd.beta <- 0.01
proposal.sd.phi <- 0.1
proposal.sd.delta <- 0.1
proposal.sd.alpha <- 0.1
proposal.sd.rho <- 0.02
proposal.sd.lambda <- 0.02
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
        Y.DA[which.miss==0] <- rbinom(n=n.miss, size=trials[which.miss==0], prob=prob[which.miss==0])
        failures.DA <- trials - Y.DA
        }else
        {}
    Y.DA.mat <- matrix(Y.DA, nrow=K, ncol=N, byrow=FALSE)
    failures.DA.mat <- matrix(failures.DA, nrow=K, ncol=N, byrow=FALSE)        
        
        
        
    ####################
    ## Sample from beta
    ####################
    offset.temp <- offset + as.numeric(phi.mat) + as.numeric(delta.time.mat) + as.numeric(alpha * time.mat)      
        if(MALA)
        {
        temp <- binomialbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y.DA, failures.DA, trials, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }else
        {
        temp <- binomialbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y.DA, failures.DA, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
        }
    beta <- temp[[1]]
    accept[1] <- accept[1] + temp[[2]]
    accept[2] <- accept[2] + n.beta.block  
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)           

        
    
    ####################
    ## Sample from alpha
    ####################
    proposal.alpha <- rnorm(n=1, mean=alpha, sd=proposal.sd.alpha)
    prob1 <- 0.5 * (alpha - prior.mean.alpha)^2 / prior.var.alpha - 0.5 * (proposal.alpha - prior.mean.alpha)^2 / prior.var.alpha
    lp.current <- offset + as.numeric(regression.mat) + as.numeric(phi.mat) + as.numeric(delta.time.mat) + as.numeric(alpha * time.mat)     
    lp.proposal <- offset + as.numeric(regression.mat) + as.numeric(phi.mat) + as.numeric(delta.time.mat) + as.numeric(proposal.alpha * time.mat)            
    p.current <- exp(lp.current) / (1 + exp(lp.current))
    p.proposal <- exp(lp.proposal) / (1 + exp(lp.proposal))
    like.current <- Y.DA * log(p.current) + failures.DA * log(1-p.current)
    like.proposal <- Y.DA * log(p.proposal) + failures.DA * log(1-p.proposal)
    prob2 <- sum(like.proposal - like.current, na.rm=TRUE)
    prob <- exp(prob1 + prob2)
        if(prob > runif(1))
        {
        alpha <- proposal.alpha
        accept[3] <- accept[3] + 1           
        }else
        {}              
    accept[4] <- accept[4] + 1           
        
        
        
    ####################
    ## Sample from phi
    ####################
    phi.offset <- offset.mat + regression.mat + delta.time.mat + alpha * time.mat
    temp1 <- binomialcarupdateRW(W.triplet, W.begfin, W.triplet.sum, K, phi, tau2.phi,Y.DA.mat, failures.DA.mat, proposal.sd.phi, rho, phi.offset, N, rep(1,N))
    phi <- temp1[[1]]
        if(rho<1)
        {
        phi <- phi - mean(phi)
        }else
        {
        phi[which(islands==1)] <- phi[which(islands==1)] - mean(phi[which(islands==1)])   
        }
    phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)    
    accept[5] <- accept[5] + temp1[[2]]
    accept[6] <- accept[6] + K  
        
        
        
    ####################
    ## Sample from delta
    ####################
    delta.offset <- offset.mat + regression.mat + phi.mat +  alpha * time.mat
    temp2 <- binomialcarupdateRW(W.triplet, W.begfin, W.triplet.sum, K, delta, tau2.delta,Y.DA.mat, failures.DA.mat, proposal.sd.delta, lambda, delta.offset, N, time)
    delta <- temp2[[1]]
        if(lambda <1)
        {
        delta <- delta - mean(delta)
        }else
        {
        delta[which(islands==1)] <- delta[which(islands==1)] - mean(delta[which(islands==1)])   
        }
    delta.time.mat <- apply(time.mat, 2, "*", delta)
    accept[7] <- accept[7] + temp2[[2]]
    accept[8] <- accept[8] + K      
        

        
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
        accept[9] <- accept[9] + 1           
        }else
        {
        }              
        accept[10] <- accept[10] + 1           
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
        accept[11] <- accept[11] + 1           
        }else
        {
        }              
        accept[12] <- accept[12] + 1           
        }else
        {}
        
        
        
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(offset.mat + regression.mat + phi.mat + delta.time.mat + alpha * time.mat)
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
        samples.phi[ele, ] <- phi
        samples.delta[ele, ] <- delta
        samples.alpha[ele, ] <- alpha
            if(!fix.rho.int) samples.rho[ele, ] <- rho
            if(!fix.rho.slo) samples.lambda[ele, ] <- lambda
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
        #### Update the proposal sds
            if(p>2)
            {
            proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 40, 50)
            }else
            {
            proposal.sd.beta <- common.accceptrates1(accept[1:2], proposal.sd.beta, 30, 40)    
            }
        proposal.sd.alpha <- common.accceptrates1(accept[3:4], proposal.sd.alpha, 30, 40) 
        proposal.sd.phi <- common.accceptrates1(accept[5:6], proposal.sd.phi, 40, 50)
        proposal.sd.delta <- common.accceptrates1(accept[7:8], proposal.sd.delta, 40, 50)
            if(!fix.rho.int) proposal.sd.rho <- common.accceptrates2(accept[9:10], proposal.sd.rho, 40, 50, 0.5)
            if(!fix.rho.slo) proposal.sd.lambda <- common.accceptrates2(accept[11:12], proposal.sd.lambda, 40, 50, 0.5)
        accept <- rep(0,12)
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
chain.results <- list(samples.beta=samples.beta, samples.alpha=samples.alpha, samples.phi=samples.phi, samples.delta=samples.delta, samples.tau2=samples.tau2, samples.rho=samples.rho, samples.lambda=samples.lambda, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    samples.Y=samples.Y, accept=accept)


#### Return the results
return(chain.results)
}