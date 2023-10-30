poisson.CARlocalisedMCMC <- function(Y, offset, X.standardised, W, G, Gstar, K, N, N.all, p, burnin, n.sample, thin, MALA, prior.mean.beta, prior.var.beta, prior.delta, prior.tau2, verbose, chain)
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
#### Compute the blocking structure for beta     
    if(!is.null(X.standardised))
    {
    ## Compute the blocking structure for beta     
    block.temp <- common.betablock(p)
    beta.beg  <- block.temp[[1]]
    beta.fin <- block.temp[[2]]
    n.beta.block <- block.temp[[3]]
    list.block <- as.list(rep(NA, n.beta.block*2))
        for(r in 1:n.beta.block)
        {
        list.block[[r]] <- beta.beg[r]:beta.fin[r]-1
        list.block[[r+n.beta.block]] <- length(list.block[[r]])
        }
    }else
    {}
    

#### Compute a starting value for beta
    if(!is.null(X.standardised))
    {
    mod.glm <- glm(Y~X.standardised-1, offset=offset, family="quasipoisson")
    beta.mean <- mod.glm$coefficients
    beta.sd <- sqrt(diag(summary(mod.glm)$cov.scaled))
    beta <- rnorm(n=length(beta.mean), mean=beta.mean, sd=beta.sd)
    regression.vec <- X.standardised %*% beta
    }else
    {
    regression.vec <- rep(0, N.all)    
    }

    
#### Generate the initial parameter values
log.Y <- log(Y)
log.Y[Y==0] <- -0.1   
res.temp <- log.Y - regression.vec - offset
clust <- kmeans(res.temp,G)
lambda <- clust$centers[order(clust$centers)]
lambda.mat <- matrix(rep(lambda, N), nrow=N, byrow=TRUE)
Z <- rep(1, N.all)
    for(j in 2:G)
    {
    Z[clust$cluster==order(clust$centers)[j]] <- j    
    }
Z.mat <- matrix(Z, nrow=K, ncol=N, byrow=FALSE)
mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)
res.sd <- sd(res.temp, na.rm=TRUE)/5
phi.mat <- matrix(rnorm(n=N.all, mean=0, sd = res.sd), nrow=K, byrow=FALSE)
phi <- as.numeric(phi.mat)
tau2 <- var(phi)/10
gamma <- runif(1)
delta <- runif(1,1, min(2, prior.delta))


#### Specify matrix quantities
Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE) 
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE)    
regression.mat <- matrix(regression.vec, nrow=K, ncol=N, byrow=FALSE)


#### Matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.Z <- array(NA, c(n.keep, N.all))
samples.lambda <- array(NA, c(n.keep, G))
samples.delta <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 1))
samples.gamma <- array(NA, c(n.keep, 1))
samples.phi <- array(NA, c(n.keep, N.all))
samples.fitted <- array(NA, c(n.keep, N.all))
samples.loglike <- array(NA, c(n.keep, N.all))


#### Specify the Metropolis quantities
    if(!is.null(X.standardised))
    {
    samples.beta <- array(NA, c(n.keep, p))
    accept <- rep(0,8)
    proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
    chol.proposal.corr.beta <- chol(proposal.corr.beta) 
    proposal.sd.beta <- 0.01
    }else
    {
    accept <- rep(0,6)    
    }
    
proposal.sd.lambda <- 0.1
proposal.sd.delta <- 0.1
proposal.sd.phi <- 0.1
Y.extend <- matrix(rep(Y, G), byrow=F, ncol=G)
delta.update <- matrix(rep(1:G, N.all-K), ncol=G, byrow=T)
tau2.posterior.shape <- prior.tau2[1] + N * (K-1) /2


#### CAR quantities
W.quants <- common.Wcheckformat.leroux(W)
W <- W.quants$W
W.triplet <- W.quants$W.triplet
W.n.triplet <- W.quants$n.triplet
W.triplet.sum <- W.quants$W.triplet.sum
n.neighbours <- W.quants$n.neighbours 
W.begfin <- W.quants$W.begfin


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
    ####################
    ## Sample from beta
    ####################
        if(!is.null(X.standardised))
        {
        offset.temp <- offset + as.numeric(mu) + as.numeric(phi.mat)   
            if(MALA)
            {
            temp <- poissonbetaupdateMALA(X.standardised, N.all, p, beta, offset.temp, Y, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
            }else
            {
            temp <- poissonbetaupdateRW(X.standardised, N.all, p, beta, offset.temp, Y, prior.mean.beta, prior.var.beta, n.beta.block, proposal.sd.beta, list.block)
            }
        beta <- temp[[1]]
        accept[7] <- accept[7] + temp[[2]]
        accept[8] <- accept[8] + n.beta.block  
        regression.vec <- X.standardised %*% beta
        regression.mat <- matrix(regression.vec, nrow=K, ncol=N, byrow=FALSE)  
        }else
        {}
        
        
        
    #######################     
    #### Sample from lambda
    #######################
    #### Propose a new value
    proposal.extend <- c(-100, lambda, 100) 
        for(r in 1:G)
        {
        proposal.extend[(r+1)] <- rtruncnorm(n=1, a=proposal.extend[r], b=proposal.extend[(r+2)], mean=proposal.extend[(r+1)], sd=proposal.sd.lambda)
        }
    proposal <- proposal.extend[-c(1, (G+2))]
        
    #### Compute the data likelihood
    lp.current <- lambda[Z] + offset + as.numeric(regression.mat) + as.numeric(phi.mat)       
    lp.proposal <- proposal[Z] + offset + as.numeric(regression.mat) + as.numeric(phi.mat)     
    like.current <- Y * lp.current - exp(lp.current)
    like.proposal <- Y * lp.proposal - exp(lp.proposal)
    prob <- exp(sum(like.proposal - like.current))
        if(prob > runif(1))
        {
        lambda <- proposal
        lambda.mat <- matrix(rep(lambda, N), nrow=N, byrow=TRUE)
        mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)
        accept[1] <- accept[1] + 1  
        }else
        {}
    accept[2] <- accept[2] + 1           
        
        
        
    ##################     
    #### Sample from Z
    ##################
    prior.offset <- rep(NA, G)
        for(r in 1:G)
        {
        prior.offset[r] <-  log(sum(exp(-delta * ((1:G - r)^2 + (1:G - Gstar)^2)))) 
        }
    mu.offset <- exp(offset.mat + regression.mat + phi.mat)
    test <- Zupdatesqpoi(Z=Z.mat, Offset=mu.offset, Y=Y.mat, delta=delta, lambda=lambda, nsites=K, ntime=N, G=G, SS=1:G, prioroffset=prior.offset, Gstar=Gstar)          
    Z.mat <- test
    Z <- as.numeric(Z.mat)
    mu <- matrix(lambda[Z], nrow=K, ncol=N, byrow=FALSE)
         
        
        
    ######################
    #### Sample from delta
    ######################
    proposal.delta <-  rtruncnorm(n=1, a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)
    sum.delta1 <- sum((Z - Gstar)^2)
    sum.delta2 <- sum((Z.mat[ ,-1] - Z.mat[ ,-N])^2)
    current.fc1 <- -delta * (sum.delta1 + sum.delta2) - K *  log(sum(exp(-delta * (1:G - Gstar)^2))) 
    proposal.fc1 <- -proposal.delta * (sum.delta1 + sum.delta2) - K *  log(sum(exp(-proposal.delta * (1:G - Gstar)^2)))                 
    Z.temp <- matrix(rep(as.numeric(Z.mat[ ,-N]),G), ncol=G, byrow=FALSE)
    Z.temp2 <- (delta.update - Z.temp)^2 + (delta.update - Gstar)^2
    current.fc <- current.fc1 - sum(log(apply(exp(-delta * Z.temp2),1,sum)))
    proposal.fc <- proposal.fc1 - sum(log(apply(exp(-proposal.delta * Z.temp2),1,sum)))    
    hastings <- log(dtruncnorm(x=delta, a=1, b=prior.delta, mean=proposal.delta, sd=proposal.sd.delta)) - log(dtruncnorm(x=proposal.delta, a=1, b=prior.delta, mean=delta, sd=proposal.sd.delta)) 
    prob <- exp(proposal.fc - current.fc + hastings)       
        if(prob > runif(1))
        {
        delta <- proposal.delta
        accept[3] <- accept[3] + 1  
        }else
        {}
    accept[4] <- accept[4] + 1   
        
        
        
    ####################
    #### Sample from phi
    ####################
    phi.offset <- mu + offset.mat + regression.mat
    temp1 <- poissonar1carupdateRW(W.triplet, W.begfin, W.triplet.sum,  K, N, phi.mat, tau2, gamma, 1, Y.mat, proposal.sd.phi, phi.offset, W.triplet.sum)      
    phi.temp <- temp1[[1]]
    phi <- as.numeric(phi.temp)
        for(i in 1:G)
        {
        phi[which(Z==i)] <- phi[which(Z==i)] - mean(phi[which(Z==i)])
        }
    phi.mat <- matrix(phi, nrow=K, ncol=N, byrow=FALSE)
    accept[5] <- accept[5] + temp1[[2]]
    accept[6] <- accept[6] + K*N    
        
        
        
    ####################
    ## Sample from gamma
    ####################
    temp2 <- gammaquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, 1)
    mean.gamma <- temp2[[1]] / temp2[[2]]
    sd.gamma <- sqrt(tau2 / temp2[[2]]) 
    gamma <- rtruncnorm(n=1, a=0, b=1, mean=mean.gamma, sd=sd.gamma)   
        
        
        
    ####################
    ## Samples from tau2
    ####################
    temp3 <- tauquadformcompute(W.triplet, W.triplet.sum, W.n.triplet,  K, N, phi.mat, 1, gamma)
    tau2.posterior.scale <- temp3 + prior.tau2[2] 
    tau2 <- 1 / rgamma(1, tau2.posterior.shape, scale=(1/tau2.posterior.scale))          
        
        
        
    #########################
    ## Calculate the deviance
    #########################
    lp <- as.numeric(mu + offset.mat + regression.mat + phi.mat)  
    fitted <- exp(lp)
    loglike <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)
 
        
    ###################
    ## Save the results
    ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
        ele <- (j - burnin) / thin
        samples.delta[ele, ] <- delta
        samples.lambda[ele, ] <- lambda
        samples.Z[ele, ] <- Z
        samples.phi[ele, ] <- as.numeric(phi.mat)
        samples.tau2[ele, ] <- tau2
        samples.gamma[ele, ] <- gamma
        samples.fitted[ele, ] <- fitted
        samples.loglike[ele, ] <- loglike
            if(!is.null(X.standardised)) samples.beta[ele, ] <- beta        
        }else
        {}
        
        
    
    ########################################
    ## Self tune the acceptance probabilties
    ########################################
        if(ceiling(j/100)==floor(j/100) & j < burnin)
        {
            if(!is.null(X.standardised))
            {
                if(p>2)
                {
                proposal.sd.beta <- common.accceptrates1(accept[7:8], proposal.sd.beta, 40, 50)
                }else
                {
                proposal.sd.beta <- common.accceptrates1(accept[7:8], proposal.sd.beta, 30, 40)    
                }
            proposal.sd.phi <- common.accceptrates1(accept[5:6], proposal.sd.phi, 40, 50)
            proposal.sd.lambda <- common.accceptrates2(accept[1:2], proposal.sd.lambda, 20, 40, 10)
            proposal.sd.delta <- common.accceptrates2(accept[3:4], proposal.sd.delta, 40, 50, prior.delta/6)
            accept <- rep(0,8)  
            }else
            {
            proposal.sd.phi <- common.accceptrates1(accept[5:6], proposal.sd.phi, 40, 50)
            proposal.sd.lambda <- common.accceptrates2(accept[1:2], proposal.sd.lambda, 20, 40, 10)
            proposal.sd.delta <- common.accceptrates2(accept[3:4], proposal.sd.delta, 40, 50, prior.delta/6)
            accept <- rep(0,6)     
            }
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
    if(is.null(X.standardised)) samples.beta <- NA
chain.results <- list(samples.beta=samples.beta, samples.phi=samples.phi, samples.Z=samples.Z, samples.lambda=samples.lambda, samples.tau2=samples.tau2, samples.delta=samples.delta, samples.gamma=samples.gamma, samples.loglike=samples.loglike, samples.fitted=samples.fitted,
                    accept=accept)

#### Return the results
return(chain.results)
}