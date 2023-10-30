poisson.CARsepspatial <- function(formula, data=NULL, W, burnin, n.sample, thin=1,  n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, rho.S=NULL, rho.T=NULL, MALA=TRUE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)  
  

#### Frame object
frame.results <- common.frame(formula, data, "poisson")
N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
n.miss <- frame.results$n.miss  
    if(n.miss>0) stop("the response has missing 'NA' values.", call.=FALSE)


#### Determine the number of spatial and temporal units
W.quants <- common.Wcheckformat.leroux(W)
K <- W.quants$n
N <- N.all / K
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 


#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE) 
  

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
    lambda <- runif(1)
    fix.rho.T <- FALSE   
    }else
    {
    lambda <- rho.T
    fix.rho.T <- TRUE
    }
    if(!is.numeric(lambda)) stop("rho.T is fixed but is not numeric.", call.=FALSE)  
    if(lambda<0 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  
    if(lambda>1 ) stop("rho.T is outside the range [0, 1].", call.=FALSE)  


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
prior.beta.check(prior.mean.beta, prior.var.beta, p)
prior.var.check(prior.tau2)

  
#### Compute the blocking structure for beta     
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
  

#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)      
  
  
  
########################
#### Run the MCMC chains
########################
   if(n.chains==1)
   {
   #### Only 1 chain
   results <- poisson.CARsepspatialMCMC(Y=Y, offset=offset, X.standardised=X.standardised, W=W, rho=rho, lambda=lambda, fix.rho.S=fix.rho.S, fix.rho.T=fix.rho.T, K=K, N=N, N.all=N.all, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- poisson.CARsepspatialMCMC(Y=Y, offset=offset, X.standardised=X.standardised, W=W, rho=rho, lambda=lambda, fix.rho.S=fix.rho.S, fix.rho.T=fix.rho.T, K=K, N=N, N.all=N.all, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, verbose=verbose, chain=i)
         }
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores>1 & ceiling(n.cores)==floor(n.cores))   
   {   
   #### Multiple chains in parallel
   results <- as.list(rep(NA, n.chains))
      if(verbose)
      {
      compclust <- makeCluster(n.cores, outfile="CARBayesSTprogress.txt")
      cat("The current progress of the model fitting algorithm has been output to CARBayesSTprogress.txt in the working directory")
      }else
      {
      compclust <- makeCluster(n.cores)
      }
   results <- clusterCall(compclust, fun=poisson.CARsepspatialMCMC, Y=Y, offset=offset, X.standardised=X.standardised, W=W, rho=rho, lambda=lambda, fix.rho.S=fix.rho.S, fix.rho.T=fix.rho.T, K=K, N=N, N.all=N.all, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, verbose=verbose, chain="all")
   stopCluster(compclust)
   }else
   {
   stop("n.chains or n.cores are not positive integers.", call.=FALSE)  
   }


#### end timer
    if(verbose)
    {
    cat("\nSummarising results.\n")
    }else
    {}

  
  
###################################
#### Summarise and save the results 
###################################
    if(n.chains==1)
    {
    #### If n.chains==1
    ## Compute the acceptance rates
    accept.final <- rep(NA, 5)
    names(accept.final) <- c("beta", "phi", "delta", "rho.S", "rho.T")
    accept.final[1] <- 100 * results$accept[1] / results$accept[2]
    accept.final[2] <- 100 * results$accept[3] / results$accept[4]
    accept.final[3] <- 100 * results$accept[7] / results$accept[8]
        if(!fix.rho.S) accept.final[4] <- 100 * results$accept[5] / results$accept[6]
        if(!fix.rho.T) accept.final[5] <- 100 * results$accept[9] / results$accept[10]

    ## Compute the fitted deviance
    mean.beta <- apply(results$samples.beta,2,mean)
    regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)   
    mean.phi <- matrix(apply(results$samples.phi, 2, mean), nrow=K, ncol=N)
    mean.delta <- apply(results$samples.delta,2,mean)
    delta.mat <- matrix(mean.delta, nrow=K, ncol=N, byrow=TRUE)
    fitted.mean <- as.numeric(exp(offset.mat + mean.phi + regression.mat + delta.mat))
    deviance.fitted <- -2 * sum(dpois(x=as.numeric(Y), lambda=fitted.mean, log=TRUE), na.rm=TRUE)
    modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

    ## Create the fitted values and residuals
    fitted.values <- apply(results$samples.fitted, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values)
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

    ## Transform the parameters back to the origianl covariate scale.
    samples.beta.orig <- common.betatransform(results$samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    
    ## Create the samples object
        if(fix.rho.S & fix.rho.T)
        {
        samples.rhoext <- NA
        }else if(fix.rho.S & !fix.rho.T)
        {
        samples.rhoext <- results$samples.lambda
        names(samples.rhoext) <- "rho.T"
        }else if(!fix.rho.S & fix.rho.T)
        {
        samples.rhoext <- results$samples.rho  
        names(samples.rhoext) <- "rho.S"
        }else
        {
        samples.rhoext <- cbind(results$samples.rho, results$samples.lambda)
        colnames(samples.rhoext) <- c("rho.S", "rho.T")
        }
     samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(results$samples.phi), delta=mcmc(results$samples.delta), tau2=mcmc(results$samples.tau2), tau2.T=mcmc(results$samples.sig2), rho=mcmc(samples.rhoext),  fitted=mcmc(results$samples.fitted))

    ## Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975))))     
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.final[names(accept.final)=="beta"],p), effectiveSize(samples$beta), geweke.diag(samples$beta)$z)
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")

    summary.hyper <- array(NA, c(3 + N, 7))    
        for (tt in 1:N)
        {
        summary.hyper[tt,1:3] <- c(mean(results$samples.tau2[, tt]), quantile(results$samples.tau2[, tt], c(0.025, 0.975)))
        summary.hyper[tt, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(results$samples.tau2[, tt])), geweke.diag(mcmc(results$samples.tau2[, tt]))$z) 
        }
    summary.hyper[N+1,1:3] <- c(mean(results$samples.sig2), quantile(results$samples.sig2, c(0.025, 0.975)))
    summary.hyper[N+1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(results$samples.sig2)), geweke.diag(mcmc(results$samples.sig2))$z)  
  
        if(!fix.rho.S)
        {
        summary.hyper[N+2, 1:3] <- c(mean(results$samples.rho), quantile(results$samples.rho, c(0.025, 0.975)))
        summary.hyper[N+2, 4:7] <- c(n.keep, accept.final[names(accept.final)=="rho.S"], effectiveSize(results$samples.rho), geweke.diag(results$samples.rho)$z)
        }else
        {
        summary.hyper[N+2, 1:3] <- c(rho, rho, rho)
        summary.hyper[N+2, 4:7] <- rep(NA, 4)
        }
        if(!fix.rho.T)
        {
        summary.hyper[N+3, 1:3] <- c(mean(results$samples.lambda), quantile(results$samples.lambda, c(0.025, 0.975)))
        summary.hyper[N+3, 4:7] <- c(n.keep, accept.final[names(accept.final)=="rho.T"], effectiveSize(results$samples.lambda), geweke.diag(results$samples.lambda)$z)
        }else
        {
        summary.hyper[N+3, 1:3] <- c(lambda, lambda, lambda)
        summary.hyper[N+3, 4:7] <- rep(NA, 4)
        }    
    rownames(summary.hyper) <- c(paste("tau2.", c(1:N), sep = ""), "tau2.T", "rho.S","rho.T")  
    summary.results <- rbind(summary.beta, summary.hyper)
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    }else
    {
    #### If n.chains > 1
    ## Compute the acceptance rates
    accept.final <- rep(NA, 5)
    names(accept.final) <- c("beta", "phi", "delta", "rho.S", "rho.T")
    accept.temp <- lapply(results, function(l) l[["accept"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.final[1] <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
    accept.final[2] <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])
    accept.final[3] <- 100 * sum(accept.temp2[ ,7]) / sum(accept.temp2[ ,8])
        if(!fix.rho.S) accept.final[4] <- 100 * sum(accept.temp2[ ,5]) / sum(accept.temp2[ ,6])
        if(!fix.rho.T) accept.final[5] <- 100 * sum(accept.temp2[ ,9]) / sum(accept.temp2[ ,10])

    ## Extract the samples into separate matrix and list objects
    samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)   
    samples.phi.list <- lapply(results, function(l) l[["samples.phi"]])
    samples.phi.matrix <- do.call(what=rbind, args=samples.phi.list)
    samples.delta.list <- lapply(results, function(l) l[["samples.delta"]])
    samples.delta.matrix <- do.call(what=rbind, args=samples.delta.list)
        if(!fix.rho.S)
        {
        samples.rho.list <- lapply(results, function(l) l[["samples.rho"]])
        samples.rho.matrix <- do.call(what=rbind, args=samples.rho.list)
        }
        if(!fix.rho.T)
        {
        samples.lambda.list <- lapply(results, function(l) l[["samples.lambda"]])
        samples.lambda.matrix <- do.call(what=rbind, args=samples.lambda.list)
        }
    samples.tau2.list <- lapply(results, function(l) l[["samples.tau2"]])
    samples.tau2.matrix <- do.call(what=rbind, args=samples.tau2.list)
    samples.sig2.list <- lapply(results, function(l) l[["samples.sig2"]])
    samples.sig2.matrix <- do.call(what=rbind, args=samples.sig2.list)   
    samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
    samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
    samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
    samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)

    ## Compute the fitted deviance
    mean.beta <- apply(samples.beta.matrix,2,mean)
    regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)   
    mean.phi <- matrix(apply(samples.phi.matrix, 2, mean), nrow=K, ncol=N)
    mean.delta <- apply(samples.delta.matrix,2,mean)
    delta.mat <- matrix(mean.delta, nrow=K, ncol=N, byrow=TRUE)
    fitted.mean <- as.numeric(exp(offset.mat + mean.phi + regression.mat + delta.mat))
    deviance.fitted <- -2 * sum(dpois(x=as.numeric(Y), lambda=fitted.mean, log=TRUE), na.rm=TRUE)
    modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)

    ## Create the fitted values and residuals
    fitted.values <- apply(samples.fitted.matrix, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values)
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

    ## Transform the parameters back to the original covariate scale.
    samples.beta.list <- samples.beta.list
        for(j in 1:n.chains)
        {
        samples.beta.list[[j]] <- common.betatransform(samples.beta.list[[j]], X.indicator, X.mean, X.sd, p, FALSE)  
        }
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)

    ## Create MCMC objects
    beta.mcmc <- mcmc.list(lapply(samples.beta.list, mcmc))
    phi.mcmc <- mcmc.list(lapply(samples.phi.list, mcmc))
    delta.mcmc <- mcmc.list(lapply(samples.delta.list, mcmc))
    fitted.mcmc <- mcmc.list(lapply(samples.fitted.list, mcmc))
    tau2.mcmc <- mcmc.list(lapply(samples.tau2.list, mcmc))    
    sig2.mcmc <- mcmc.list(lapply(samples.sig2.list, mcmc))    
        if(fix.rho.S & fix.rho.T)
        {
        rhoext.mcmc <- NA
        }else if(fix.rho.S & !fix.rho.T)
        {
            for(j in 1:n.chains)
            {
            colnames(samples.lambda.list[[j]]) <- c("rho.T")  
            } 
        rhoext.mcmc <- mcmc.list(lapply(samples.lambda.list, mcmc))
        }else if(!fix.rho.S & fix.rho.T)
        {
            for(j in 1:n.chains)
            {
            colnames(samples.rho.list[[j]]) <- c("rho.S")  
            } 
        rhoext.mcmc <- mcmc.list(lapply(samples.rho.list, mcmc))
        }else
        {
        rho.temp <- as.list(rep(NA, n.chains))    
            for(j in 1:n.chains)
            {
            rho.temp[[j]] <- cbind(samples.rho.list[[j]], samples.lambda.list[[j]])
            colnames(rho.temp[[j]]) <- c("rho.S", "rho.T")
            }
        rhoext.mcmc <- mcmc.list(lapply(rho.temp, mcmc))
        }
    samples <- list(beta=beta.mcmc, phi=phi.mcmc, delta=delta.mcmc, rho=rhoext.mcmc, tau2=tau2.mcmc, tau2.T=sig2.mcmc, fitted=fitted.mcmc)
    
    ## create a summary object
    n.keep <- floor((n.sample - burnin)/thin) * n.chains
    summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.final[names(accept.final)=="beta"],p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")

    summary.hyper <- array(NA, c(3 + N, 7))    
        for (tt in 1:N)
        {
        summary.hyper[tt,1:3] <- c(mean(samples.tau2.matrix[, tt]), quantile(samples.tau2.matrix[, tt], c(0.025, 0.975)))
        summary.hyper[tt, 4:7] <- c(n.keep, 100, effectiveSize(tau2.mcmc[ ,tt]), gelman.diag(tau2.mcmc[ ,tt])$psrf[ ,2]) 
        }
    summary.hyper[N+1,1:3] <- c(mean(samples.sig2.matrix), quantile(samples.sig2.matrix, c(0.025, 0.975)))
    summary.hyper[N+1, 4:7] <- c(n.keep, 100, effectiveSize(sig2.mcmc), gelman.diag(sig2.mcmc)$psrf[ ,2])  
  
        if(!fix.rho.S)
        {
        temp <- mcmc.list(lapply(samples.rho.list, mcmc))
        summary.hyper[N+2, 1:3] <- c(mean(samples.rho.matrix), quantile(samples.rho.matrix, c(0.025, 0.975)))
        summary.hyper[N+2, 4:7] <- c(n.keep, accept.final[names(accept.final)=="rho.S"], effectiveSize(temp), gelman.diag(temp)$psrf[ ,2])
        }else
        {
        summary.hyper[N+2, 1:3] <- c(rho, rho, rho)
        summary.hyper[N+2, 4:7] <- rep(NA, 4)
        }
        if(!fix.rho.T)
        {
        temp <- mcmc.list(lapply(samples.lambda.list, mcmc))
        summary.hyper[N+3, 1:3] <- c(mean(samples.lambda.matrix), quantile(samples.lambda.matrix, c(0.025, 0.975)))
        summary.hyper[N+3, 4:7] <- c(n.keep, accept.final[names(accept.final)=="rho.T"], effectiveSize(temp), gelman.diag(temp)$psrf[ ,2])
        }else
        {
        summary.hyper[N+3, 1:3] <- c(lambda, lambda, lambda)
        summary.hyper[N+3, 4:7] <- rep(NA, 4)
        }    
    rownames(summary.hyper) <- c(paste("tau2.", c(1:N), sep = ""), "tau2.T", "rho.S","rho.T")  
    summary.results <- rbind(summary.beta, summary.hyper)
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    }

  

###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Poisson (log link function)", "\nLatent structure model - An overall time trend with temporal specific spatial effects\n")
n.total <- floor((n.sample - burnin) / thin) * n.chains
mcmc.info <- c(n.total, n.sample, burnin, thin, n.chains)
names(mcmc.info) <- c("Total samples", "n.sample", "burnin", "thin", "n.chains")
results.final <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL, formula=formula, model=model.string,  mcmc.info=mcmc.info, X=X)
class(results.final) <- "CARBayesST"
     if(verbose)
     {
     b<-proc.time()
     cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
     }else
     {}
 return(results.final)
}