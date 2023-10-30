gaussian.CARlinear <- function(formula, data=NULL, W, burnin, n.sample, thin=1,  n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.mean.alpha=NULL, prior.var.alpha=NULL, prior.nu2=NULL, prior.tau2=NULL, rho.slo=NULL, rho.int=NULL,  verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)  
    
    
#### Frame object
frame.results <- common.frame(formula, data, "gaussian")
N.all <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
which.miss <- frame.results$which.miss
n.miss <- frame.results$n.miss  


#### Determine the number of spatial and temporal units
W.quants <- common.Wcheckformat.leroux(W)
K <- W.quants$n
N <- N.all / K
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
time <-(1:N - mean(1:N))/N
time.mat <- matrix(rep(time, K), byrow=TRUE, nrow=K)   

    
#### Check on the rho arguments
    if(is.null(rho.int))
    {
    rho <- runif(1)
    fix.rho.int <- FALSE   
    }else
    {
    rho <- rho.int
    fix.rho.int <- TRUE
    }
    if(!is.numeric(rho)) stop("rho.int is fixed but is not numeric.", call.=FALSE)  
    if(rho<0 ) stop("rho.int is outside the range [0, 1].", call.=FALSE)  
    if(rho>1 ) stop("rho.int is outside the range [0, 1].", call.=FALSE)    

    if(is.null(rho.slo))
    {
    lambda <- runif(1)
    fix.rho.slo <- FALSE   
    }else
    {
    lambda <- rho.slo
    fix.rho.slo <- TRUE
    }
    if(!is.numeric(lambda)) stop("rho.slo is fixed but is not numeric.", call.=FALSE)  
    if(lambda<0 ) stop("rho.slo is outside the range [0, 1].", call.=FALSE)  
    if(lambda>1 ) stop("rho.slo is outside the range [0, 1].", call.=FALSE)  


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(1, 0.01)
    if(is.null(prior.nu2)) prior.nu2 <- c(1, 0.01)
    if(is.null(prior.mean.alpha)) prior.mean.alpha <- rep(0, 1)
    if(is.null(prior.var.alpha)) prior.var.alpha <- rep(100000, 1)
prior.beta.check(prior.mean.beta, prior.var.beta, p)
prior.var.check(prior.tau2)
prior.var.check(prior.nu2)
    if(length(prior.mean.alpha)!=1) stop("the prior mean for alpha is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.alpha)) stop("the  prior mean for alpha is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.alpha))!=0) stop("the prior mean for alpha has missing values.", call.=FALSE)       
    if(length(prior.var.alpha)!=1) stop("the prior variance for alpha is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.alpha)) stop("the  prior variance for alpha is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.alpha))!=0) stop("the  prior variance for alpha has missing values.", call.=FALSE)    
    if(min(prior.var.alpha) <=0) stop("the prior variance for alpha has elements less than zero", call.=FALSE)

    
#### MCMC quantities - burnin, n.sample, thin
common.burnin.nsample.thin.check(burnin, n.sample, thin)



########################
#### Run the MCMC chains
########################
   if(n.chains==1)
   {
   #### Only 1 chain
   results <- gaussian.CARlinearMCMC(Y=Y, offset=offset, X.standardised=X.standardised, W=W, rho=rho, lambda=lambda, fix.rho.int=fix.rho.int, fix.rho.slo=fix.rho.slo, K=K, N=N, N.all=N.all, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.alpha=prior.mean.alpha, prior.var.alpha=prior.var.alpha, prior.tau2=prior.tau2, prior.nu2=prior.nu2, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- gaussian.CARlinearMCMC(Y=Y, offset=offset, X.standardised=X.standardised, W=W, rho=rho, lambda=lambda, fix.rho.int=fix.rho.int, fix.rho.slo=fix.rho.slo, K=K, N=N, N.all=N.all, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.alpha=prior.mean.alpha, prior.var.alpha=prior.var.alpha, prior.tau2=prior.tau2, prior.nu2=prior.nu2, verbose=verbose, chain=i)
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
   results <- clusterCall(compclust, fun=gaussian.CARlinearMCMC, Y=Y, offset=offset, X.standardised=X.standardised, W=W, rho=rho, lambda=lambda, fix.rho.int=fix.rho.int, fix.rho.slo=fix.rho.slo, K=K, N=N, N.all=N.all, p=p, which.miss=which.miss, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.alpha=prior.mean.alpha, prior.var.alpha=prior.var.alpha, prior.tau2=prior.tau2, prior.nu2=prior.nu2, verbose=verbose, chain="all")
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
    accept.final <- rep(NA, 6)
    names(accept.final) <- c("beta", "alpha", "phi", "delta", "rho.int", "rho.slo")
    accept.final[1:4] <- 100
        if(!fix.rho.int) accept.final[5] <- 100 * results$accept[1] / results$accept[2]
        if(!fix.rho.slo) accept.final[6] <- 100 * results$accept[3] / results$accept[4]

    ## Compute the fitted deviance
    mean.phi <- apply(results$samples.phi, 2, mean)
    mean.delta <- apply(results$samples.delta, 2, mean)  
    mean.phi.mat <- matrix(rep(mean.phi, N), byrow=F, nrow=K)
    delta.time.mat <- apply(time.mat, 2, "*", mean.delta)
    mean.alpha <- mean(results$samples.alpha)
    mean.beta <- apply(results$samples.beta,2,mean)
    regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)
    nu2.mean <- mean(results$samples.nu2)
    fitted.mean <- offset.mat + regression.mat + mean.phi.mat + delta.time.mat + mean.alpha * time.mat
    deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.mean, sd = rep(sqrt(nu2.mean),N.all), log = TRUE), na.rm=TRUE)    
    modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

    ## Create the fitted values and residuals
    fitted.values <- apply(results$samples.fitted, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    pearson.residuals <- response.residuals /sqrt(nu2.mean)
    residuals <- data.frame(response=response.residuals, pearson=pearson.residuals)

    ## Transform the parameters back to the origianl covariate scale.
    samples.beta.orig <- common.betatransform(results$samples.beta, X.indicator, X.mean, X.sd, p, FALSE)
    
    ## Create the samples object
        if(fix.rho.int & fix.rho.slo)
        {
        samples.rhoext <- NA
        }else if(fix.rho.int & !fix.rho.slo)
        {
        samples.rhoext <- results$samples.lambda
        names(samples.rhoext) <- "rho.slo"
        }else if(!fix.rho.int & fix.rho.slo)
        {
        samples.rhoext <- results$samples.rho  
        names(samples.rhoext) <- "rho.int"
        }else
        {
        samples.rhoext <- cbind(results$samples.rho, results$samples.lambda)
        colnames(samples.rhoext) <- c("rho.int", "rho.slo")
        }
    colnames(results$samples.tau2) <- c("tau2.int", "tau2.slo")
    samples <- list(beta=mcmc(samples.beta.orig), alpha=mcmc(results$samples.alpha), phi=mcmc(results$samples.phi), delta=mcmc(results$samples.delta), tau2=mcmc(results$samples.tau2), nu2=mcmc(results$samples.nu2), rho=mcmc(samples.rhoext),  fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))

    ## Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples$beta, 2, mean), apply(samples$beta, 2, quantile, c(0.025, 0.975))))     
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.final[names(accept.final)=="beta"],p), effectiveSize(samples$beta), geweke.diag(samples$beta)$z)
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")

    summary.tau2 <- cbind(apply(results$samples.tau2, 2, mean), t(apply(results$samples.tau2, 2, quantile, c(0.025, 0.975))), rep(n.keep, 2), rep(100, 2),
                        effectiveSize(samples$tau2), geweke.diag(samples$tau2)$z)
    summary.alpha <- c(mean(results$samples.alpha), quantile(results$samples.alpha, c(0.025, 0.975)), n.keep, accept.final[names(accept.final)=="alpha"],
                        effectiveSize(samples$alpha), geweke.diag(samples$alpha)$z)
    summary.nu2 <- c(mean(results$samples.nu2), quantile(results$samples.nu2, c(0.025, 0.975)), n.keep, 100,
                        effectiveSize(samples$nu2), geweke.diag(samples$nu2)$z)
    summary.combine <- rbind(summary.alpha, summary.nu2, summary.tau2)
    rownames(summary.combine) <- c("alpha", "nu2", "tau2.int", "tau2.slo")

    summary.rho <- array(NA, c(2,7))
    row.names(summary.rho) <- c("rho.int", "rho.slo")
            if(!fix.rho.int)
            {
            summary.rho[1, 1:3] <- c(mean(results$samples.rho), quantile(results$samples.rho, c(0.025, 0.975)))
            summary.rho[1, 4:7] <- c(n.keep, accept.final[names(accept.final)=="rho.int"], effectiveSize(results$samples.rho), geweke.diag(results$samples.rho)$z)
            }else
            {
            summary.rho[1, 1:3] <- c(rho, rho, rho)
            summary.rho[1, 4:7] <- rep(NA, 4)
            }
            if(!fix.rho.slo)
            {
            summary.rho[2, 1:3] <- c(mean(results$samples.lambda), quantile(results$samples.lambda, c(0.025, 0.975)))
            summary.rho[2, 4:7] <- c(n.keep, accept.final[names(accept.final)=="rho.slo"], effectiveSize(results$samples.lambda), geweke.diag(results$samples.lambda)$z)
            }else
            {
            summary.rho[2, 1:3] <- c(lambda, lambda, lambda)
            summary.rho[2, 4:7] <- rep(NA, 4)
            }
    
    summary.results <- rbind(summary.beta, summary.combine, summary.rho)
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1) 
    }else
    {
    #### If n.chains > 1
    ## Compute the acceptance rates
    accept.final <- rep(NA, 6)
    names(accept.final) <- c("beta", "alpha", "phi", "delta", "rho.int", "rho.slo")
    accept.temp <- lapply(results, function(l) l[["accept"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.final[1:4] <- 100
        if(!fix.rho.int) accept.final[5] <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
        if(!fix.rho.slo) accept.final[6] <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])

    ## Extract the samples into separate matrix and list objects
    samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)   
    samples.phi.list <- lapply(results, function(l) l[["samples.phi"]])
    samples.phi.matrix <- do.call(what=rbind, args=samples.phi.list)
    samples.delta.list <- lapply(results, function(l) l[["samples.delta"]])
    samples.delta.matrix <- do.call(what=rbind, args=samples.delta.list)
    samples.alpha.list <- lapply(results, function(l) l[["samples.alpha"]])
    samples.alpha.matrix <- do.call(what=rbind, args=samples.alpha.list)
        if(!fix.rho.int)
        {
        samples.rho.list <- lapply(results, function(l) l[["samples.rho"]])
        samples.rho.matrix <- do.call(what=rbind, args=samples.rho.list)
        }
        if(!fix.rho.slo)
        {
        samples.lambda.list <- lapply(results, function(l) l[["samples.lambda"]])
        samples.lambda.matrix <- do.call(what=rbind, args=samples.lambda.list)
        }
    samples.tau2.list <- lapply(results, function(l) l[["samples.tau2"]])
    samples.tau2.matrix <- do.call(what=rbind, args=samples.tau2.list)   
    samples.nu2.list <- lapply(results, function(l) l[["samples.nu2"]])
    samples.nu2.matrix <- do.call(what=rbind, args=samples.nu2.list)   
    samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
    samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
    samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
    samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)
        if(n.miss>0) samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])

    ## Compute the fitted deviance
    mean.phi <- apply(samples.phi.matrix, 2, mean)
    mean.delta <- apply(samples.delta.matrix, 2, mean)  
    mean.phi.mat <- matrix(rep(mean.phi, N), byrow=F, nrow=K)
    delta.time.mat <- apply(time.mat, 2, "*", mean.delta)
    mean.alpha <- mean(samples.alpha.matrix)
    mean.beta <- apply(samples.beta.matrix,2,mean)
    regression.mat <- matrix(X.standardised %*% mean.beta, nrow=K, ncol=N, byrow=FALSE)
    nu2.mean <- mean(samples.nu2.matrix)
    fitted.mean <- offset.mat + regression.mat + mean.phi.mat + delta.time.mat + mean.alpha * time.mat
    deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.mean, sd = rep(sqrt(nu2.mean),N.all), log = TRUE), na.rm=TRUE)    
    modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)

    ## Create the fitted values and residuals
    fitted.values <- apply(samples.fitted.matrix, 2, mean)
    response.residuals <- as.numeric(Y) - fitted.values
    pearson.residuals <- response.residuals /sqrt(nu2.mean)
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
    alpha.mcmc <- mcmc.list(lapply(samples.alpha.list, mcmc))
    phi.mcmc <- mcmc.list(lapply(samples.phi.list, mcmc))
    delta.mcmc <- mcmc.list(lapply(samples.delta.list, mcmc))
    nu2.mcmc <- mcmc.list(lapply(samples.nu2.list, mcmc))
    fitted.mcmc <- mcmc.list(lapply(samples.fitted.list, mcmc))
        for(j in 1:n.chains)
        {
        colnames(samples.tau2.list[[j]]) <- c("tau2.int", "tau2.slo")  
        }
    tau2.mcmc <- mcmc.list(lapply(samples.tau2.list, mcmc))
        if(n.miss>0) 
        {    
        Y.mcmc <- mcmc.list(lapply(samples.Y.list, mcmc))
        }else
        {
        Y.mcmc <- NA    
        }
        if(fix.rho.int & fix.rho.slo)
        {
        rhoext.mcmc <- NA
        }else if(fix.rho.int & !fix.rho.slo)
        {
            for(j in 1:n.chains)
            {
            colnames(samples.lambda.list[[j]]) <- c("rho.slo")  
            } 
        rhoext.mcmc <- mcmc.list(lapply(samples.lambda.list, mcmc))
        }else if(!fix.rho.int & fix.rho.slo)
        {
            for(j in 1:n.chains)
            {
            colnames(samples.rho.list[[j]]) <- c("rho.int")  
            } 
        rhoext.mcmc <- mcmc.list(lapply(samples.rho.list, mcmc))
        }else
        {
        rho.temp <- as.list(rep(NA, n.chains))    
            for(j in 1:n.chains)
            {
            rho.temp[[j]] <- cbind(samples.rho.list[[j]], samples.lambda.list[[j]])
            colnames(rho.temp[[j]]) <- c("rho.int", "rho.slo")
            }
        rhoext.mcmc <- mcmc.list(lapply(rho.temp, mcmc))
        }
    samples <- list(beta=beta.mcmc, alpha=alpha.mcmc, phi=phi.mcmc, delta=delta.mcmc, rho=rhoext.mcmc, tau2=tau2.mcmc, nu2=nu2.mcmc, fitted=fitted.mcmc, Y=Y.mcmc)

    ## create a summary object
    n.keep <- floor((n.sample - burnin)/thin) * n.chains
    summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.final[names(accept.final)=="beta"],p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
    rownames(summary.beta) <- colnames(X)
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")

    summary.tau2 <- cbind(apply(samples.tau2.matrix, 2, mean), t(apply(samples.tau2.matrix, 2, quantile, c(0.025, 0.975))), rep(n.keep, 2), rep(100, 2),
                        effectiveSize(tau2.mcmc), gelman.diag(tau2.mcmc)$psrf[ ,2])
    summary.alpha <- c(mean(samples.alpha.matrix), quantile(samples.alpha.matrix, c(0.025, 0.975)), n.keep, accept.final[names(accept.final)=="alpha"],
                        effectiveSize(alpha.mcmc), gelman.diag(alpha.mcmc)$psrf[ ,2])
    summary.nu2 <- c(mean(samples.nu2.matrix), quantile(samples.nu2.matrix, c(0.025, 0.975)), n.keep, 100,
                        effectiveSize(nu2.mcmc), gelman.diag(nu2.mcmc)$psrf[ ,2])
    summary.combine <- rbind(summary.alpha, summary.nu2, summary.tau2)
    rownames(summary.combine) <- c("alpha", "nu2", "tau2.int", "tau2.slo")

    summary.rho <- array(NA, c(2,7))
    row.names(summary.rho) <- c("rho.int", "rho.slo")
            if(!fix.rho.int)
            {
            temp <- mcmc.list(lapply(samples.rho.list, mcmc))
            summary.rho[1, 1:3] <- c(mean(samples.rho.matrix), quantile(samples.rho.matrix, c(0.025, 0.975)))
            summary.rho[1, 4:7] <- c(n.keep, accept.final[names(accept.final)=="rho.int"], effectiveSize(temp), gelman.diag(temp)$psrf[ ,2])
            }else
            {
            summary.rho[1, 1:3] <- c(rho, rho, rho)
            summary.rho[1, 4:7] <- rep(NA, 4)
            }
            if(!fix.rho.slo)
            {
            temp <- mcmc.list(lapply(samples.lambda.list, mcmc))
            summary.rho[2, 1:3] <- c(mean(samples.lambda.matrix), quantile(samples.lambda.matrix, c(0.025, 0.975)))
            summary.rho[2, 4:7] <- c(n.keep, accept.final[names(accept.final)=="rho.slo"], effectiveSize(temp), gelman.diag(temp)$psrf[ ,2])
            }else
            {
            summary.rho[2, 1:3] <- c(lambda, lambda, lambda)
            summary.rho[2, 4:7] <- rep(NA, 4)
            }
    
    summary.results <- rbind(summary.beta, summary.combine, summary.rho)
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    }



###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nLatent structure model - Spatially autocorrelated linear time trends\n")
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


