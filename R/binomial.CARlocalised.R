binomial.CARlocalised <- function(formula, data=NULL, G, trials,  W, burnin, n.sample, thin=1,  n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.delta=NULL, prior.tau2=NULL, MALA=TRUE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)  
    
    
#### Frame object
frame.results <- common.frame.localised(formula, data, "binomial")
N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- as.numeric(!is.na(Y))
n.miss <- N.all - sum(which.miss)
    if(n.miss>0) stop("the response has missing 'NA' values.", call.=FALSE)



#### Determine the number of spatial and temporal units
W.quants <- common.Wcheckformat.leroux(W)
K <- W.quants$n
N <- N.all / K
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE)   


#### Check on MALA argument
    if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
    if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  


#### Check the trials argument
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- N.all-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)
failures <- trials - Y


#### Format and check the number of clusters G     
    if(length(G)!=1) stop("G is the wrong length.", call.=FALSE)    
    if(!is.numeric(G)) stop("G is not numeric.", call.=FALSE)    
    if(G<=1) stop("G is less than 2.", call.=FALSE)    
    if(G!=round(G)) stop("G is not an integer.", call.=FALSE) 
    if(floor(G/2)==ceiling(G/2))
    {
    Gstar <- G/2    
    }else
    {
    Gstar <- (G+1)/2          
    }


#### Priors
    if(!is.null(X.standardised))
    {
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    prior.beta.check(prior.mean.beta, prior.var.beta, p)
    }else
    {}
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
prior.var.check(prior.tau2)
    if(is.null(prior.delta)) prior.delta <- 10
    if(length(prior.delta)!=1) stop("the prior value for delta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.delta)) stop("the prior value for delta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.delta))!=0) stop("the prior value for delta has missing values.", call.=FALSE)    
    if(prior.delta<=0) stop("the prior value for delta is not positive.", call.=FALSE)    


#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)  



########################
#### Run the MCMC chains
########################
   if(n.chains==1)
   {
   #### Only 1 chain
   results <- binomial.CARlocalisedMCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, G=G, Gstar=Gstar, K=K, N=N, N.all=N.all, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.delta=prior.delta, prior.tau2=prior.tau2, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- binomial.CARlocalisedMCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, G=G, Gstar=Gstar, K=K, N=N, N.all=N.all, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.delta=prior.delta, prior.tau2=prior.tau2, verbose=verbose, chain=i)
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
   results <- clusterCall(compclust, fun=binomial.CARlocalisedMCMC, Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, G=G, Gstar=Gstar, K=K, N=N, N.all=N.all, p=p, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.delta=prior.delta, prior.tau2=prior.tau2, verbose=verbose, chain="all")
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
    accept.lambda <- 100 * results$accept[1] / results$accept[2]
    accept.delta <- 100 * results$accept[3] / results$accept[4]
    accept.phi <- 100 * results$accept[5] / results$accept[6]
    accept.gamma <- 100
        if(!is.null(X.standardised))
        {
        accept.beta <- 100 * results$accept[7] / results$accept[8]   
        accept.final <- c(accept.beta, accept.lambda, accept.delta, accept.phi, accept.gamma)
        names(accept.final) <- c("beta", "lambda", "delta", "phi", "rho.T")  
        }else
        {
        accept.final <- c(accept.lambda,  accept.delta, accept.phi, accept.gamma)
        names(accept.final) <- c("lambda", "delta", "phi", "rho.T")   
        }

    ## Compute the fitted deviance
    mean.Z <- round(apply(results$samples.Z,2,mean), 0)       
    mean.lambda <- apply(results$samples.lambda, 2, mean)
    mean.mu <- matrix(mean.lambda[mean.Z], nrow=K, ncol=N, byrow=FALSE)
        if(!is.null(X.standardised))
        {
        mean.beta <- apply(results$samples.beta,2,mean)
        regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)     
        }else
        {
        regression.mat <- matrix(0, nrow=K, ncol=N, byrow=FALSE)    
        }
    mean.phi <- matrix(apply(results$samples.phi, 2, mean), nrow=K, byrow=FALSE)
    lp.mean <- as.numeric(mean.mu + offset.mat + mean.phi + regression.mat)   
    mean.prob <- exp(lp.mean)  / (1 + exp(lp.mean))
    fitted.mean <- trials * mean.prob
    deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mean.prob, log=TRUE))
    modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

    ## Create the fitted values and residuals
    fitted.values <- apply(results$samples.fitted, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

    ## Transform the parameters back to the origianl covariate scale.
        if(!is.null(X.standardised))
        {    
        samples.beta.orig <- common.betatransform(results$samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
        }else
        {
        samples.beta.orig = NA    
        }
    
    ## Create the samples object
    samples <- list(beta=mcmc(samples.beta.orig), lambda=mcmc(results$samples.lambda),  Z=mcmc(results$samples.Z), delta=mcmc(results$samples.delta), phi = mcmc(results$samples.phi), tau2=mcmc(results$samples.tau2), rho.T=mcmc(results$samples.gamma), fitted=mcmc(results$samples.fitted))

    ## Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.hyper <- array(NA, c(3, 7))     
    summary.hyper[1,1:3] <- c(mean(results$samples.delta), quantile(results$samples.delta, c(0.025, 0.975)))
    summary.hyper[2,1:3] <- c(mean(results$samples.tau2), quantile(results$samples.tau2, c(0.025, 0.975)))
    summary.hyper[3,1:3] <- c(mean(results$samples.gamma), quantile(results$samples.gamma, c(0.025, 0.975)))
    rownames(summary.hyper) <- c("delta", "tau2", "rho.T")      
    summary.hyper[1, 4:7] <- c(n.keep, accept.delta, effectiveSize(mcmc(results$samples.delta)), geweke.diag(mcmc(results$samples.delta))$z)   
    summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(results$samples.tau2)), geweke.diag(mcmc(results$samples.tau2))$z)   
    summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(results$samples.gamma)), geweke.diag(mcmc(results$samples.gamma))$z)   
    
    summary.lambda <- t(rbind(apply(results$samples.lambda, 2, mean), apply(results$samples.lambda, 2, quantile, c(0.025, 0.975)))) 
    summary.lambda <- cbind(summary.lambda, rep(n.keep, G), rep(accept.lambda, G), effectiveSize(mcmc(results$samples.lambda)), geweke.diag(mcmc(results$samples.lambda))$z)
    summary.lambda <- matrix(summary.lambda, ncol=7)
    rownames(summary.lambda) <- paste("lambda", 1:G, sep="")

        if(!is.null(X.standardised))
        {
        samples.beta.orig <- mcmc(samples.beta.orig)
        summary.beta <- t(rbind(apply(samples.beta.orig, 2, mean), apply(samples.beta.orig, 2, quantile, c(0.025, 0.975)))) 
        summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
        rownames(summary.beta) <- colnames(X)
        colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")
        summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)    
        }else
        {
        summary.results <- rbind(summary.lambda, summary.hyper)    
        }
    
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    colnames(summary.results) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")    
    }else
    {
    #### If n.chains > 1
    ## Compute the acceptance rates
    accept.temp <- lapply(results, function(l) l[["accept"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.lambda <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
    accept.delta <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])
    accept.phi <- 100 * sum(accept.temp2[ ,5]) / sum(accept.temp2[ ,6])
    accept.gamma <- 100
        if(!is.null(X.standardised))
        {
        accept.beta <- 100 * sum(accept.temp2[ ,7]) / sum(accept.temp2[ ,8])  
        accept.final <- c(accept.beta, accept.lambda, accept.delta, accept.phi, accept.gamma)
        names(accept.final) <- c("beta", "lambda", "delta", "phi", "rho.T")  
        }else
        {
        accept.final <- c(accept.lambda,  accept.delta, accept.phi, accept.gamma)
        names(accept.final) <- c("lambda", "delta", "phi", "rho.T")   
        }

    ## Extract the samples into separate matrix and list objects
        if(!is.null(X.standardised))
        {
        samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
        samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)   
        }else
        {}
    samples.phi.list <- lapply(results, function(l) l[["samples.phi"]])
    samples.phi.matrix <- do.call(what=rbind, args=samples.phi.list)
    samples.Z.list <- lapply(results, function(l) l[["samples.Z"]])
    samples.Z.matrix <- do.call(what=rbind, args=samples.Z.list)
    samples.gamma.list <- lapply(results, function(l) l[["samples.gamma"]])
    samples.gamma.matrix <- do.call(what=rbind, args=samples.gamma.list)
    samples.delta.list <- lapply(results, function(l) l[["samples.delta"]])
    samples.delta.matrix <- do.call(what=rbind, args=samples.delta.list)
    samples.lambda.list <- lapply(results, function(l) l[["samples.lambda"]])
    samples.lambda.matrix <- do.call(what=rbind, args=samples.lambda.list)
    samples.tau2.list <- lapply(results, function(l) l[["samples.tau2"]])
    samples.tau2.matrix <- do.call(what=rbind, args=samples.tau2.list)   
    samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
    samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
    samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
    samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)

    ## Compute the fitted deviance
    mean.Z <- round(apply(samples.Z.matrix,2,mean), 0)       
    mean.lambda <- apply(samples.lambda.matrix, 2, mean)
    mean.mu <- matrix(mean.lambda[mean.Z], nrow=K, ncol=N, byrow=FALSE)
        if(!is.null(X.standardised))
        {
        mean.beta <- apply(samples.beta.matrix,2,mean)
        regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)     
        }else
        {
        regression.mat <- matrix(0, nrow=K, ncol=N, byrow=FALSE)    
        }
    mean.phi <- matrix(apply(samples.phi.matrix, 2, mean), nrow=K, byrow=FALSE)
    lp.mean <- as.numeric(mean.mu + offset.mat + mean.phi + regression.mat)   
    mean.prob <- exp(lp.mean)  / (1 + exp(lp.mean))
    fitted.mean <- trials * mean.prob
    deviance.fitted <- -2 * sum(dbinom(x=Y, size=trials, prob=mean.prob, log=TRUE))
    modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)

    ## Create the fitted values and residuals
    fitted.values <- apply(samples.fitted.matrix, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

    ## Transform the parameters back to the original covariate scale.
        if(!is.null(X.standardised))
        {
        samples.beta.list <- samples.beta.list
            for(j in 1:n.chains)
            {
            samples.beta.list[[j]] <- common.betatransform(samples.beta.list[[j]], X.indicator, X.mean, X.sd, p, FALSE)  
            }
        samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)
        beta.mcmc <- mcmc.list(lapply(samples.beta.list, mcmc))
        }else
        {
        beta.mcmc = mcmc(NA)    
        }
 
    ## Create MCMC objects
    phi.mcmc <- mcmc.list(lapply(samples.phi.list, mcmc))
    fitted.mcmc <- mcmc.list(lapply(samples.fitted.list, mcmc))
    tau2.mcmc <- mcmc.list(lapply(samples.tau2.list, mcmc))
    Z.mcmc <- mcmc.list(lapply(samples.Z.list, mcmc))
    gamma.mcmc <- mcmc.list(lapply(samples.gamma.list, mcmc))
    delta.mcmc <- mcmc.list(lapply(samples.delta.list, mcmc))    
    lambda.mcmc <- mcmc.list(lapply(samples.lambda.list, mcmc))    
    samples <- list(beta=beta.mcmc, phi=phi.mcmc, Z=Z.mcmc, rho.T=gamma.mcmc, lambda=lambda.mcmc, delta=delta.mcmc, tau2=tau2.mcmc, fitted=fitted.mcmc)
 
    ## create a summary object
    n.keep <- floor((n.sample - burnin)/thin) * n.chains
    summary.hyper <- array(NA, c(3, 7))     
    summary.hyper[1,1:3] <- c(mean(samples.delta.matrix), quantile(samples.delta.matrix, c(0.025, 0.975)))
    summary.hyper[2,1:3] <- c(mean(samples.tau2.matrix), quantile(samples.tau2.matrix, c(0.025, 0.975)))
    summary.hyper[3,1:3] <- c(mean(samples.gamma.matrix), quantile(samples.gamma.matrix, c(0.025, 0.975)))
    rownames(summary.hyper) <- c("delta", "tau2", "rho.T")      
    summary.hyper[1, 4:7] <- c(n.keep, accept.delta, effectiveSize(delta.mcmc), gelman.diag(delta.mcmc)$psrf[ ,2])   
    summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(tau2.mcmc), gelman.diag(tau2.mcmc)$psrf[ ,2])   
    summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(gamma.mcmc), gelman.diag(gamma.mcmc)$psrf[ ,2])   
    
    summary.lambda <- t(rbind(apply(samples.lambda.matrix, 2, mean), apply(samples.lambda.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.lambda <- cbind(summary.lambda, rep(n.keep, G), rep(accept.lambda, G), effectiveSize(lambda.mcmc), gelman.diag(lambda.mcmc)$psrf[ ,2])
    summary.lambda <- matrix(summary.lambda, ncol=7)
    rownames(summary.lambda) <- paste("lambda", 1:G, sep="")

        if(!is.null(X.standardised))
        {
        summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
        summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.final[names(accept.final)=="beta"],p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
        rownames(summary.beta) <- colnames(X)
        colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")
        summary.results <- rbind(summary.beta, summary.lambda, summary.hyper)    
        }else
        {
        summary.results <- rbind(summary.lambda, summary.hyper)    
        }
    
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    colnames(summary.results) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")    
    }



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Binomial (logit link function)", "\nLatent structure model - Localised autoregressive order 1 CAR model\n")
n.total <- floor((n.sample - burnin) / thin) * n.chains
mcmc.info <- c(n.total, n.sample, burnin, thin, n.chains)
names(mcmc.info) <- c("Total samples", "n.sample", "burnin", "thin", "n.chains")
results.final <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=mean.Z, formula=formula, model=model.string,  mcmc.info=mcmc.info, X=X)
class(results.final) <- "CARBayesST"
     if(verbose)
     {
     b<-proc.time()
     cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
     }else
     {}
 return(results.final)
}
