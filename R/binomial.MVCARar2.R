binomial.MVCARar2 <- function(formula, data=NULL,  trials, W, burnin, n.sample, thin=1, n.chains=1, n.cores=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.Sigma.df=NULL, prior.Sigma.scale=NULL, rho.S=NULL, rho.T=NULL, MALA=TRUE, verbose=TRUE)
{
##############################################
#### Format the arguments and check for errors
##############################################
#### Verbose
a <- common.verbose(verbose)
  
  
#### Check on MALA argument
  if(length(MALA)!=1) stop("MALA is not length 1.", call.=FALSE)
  if(!is.logical(MALA)) stop("MALA is not logical.", call.=FALSE)  
  
  
#### Frame object
frame.results <- common.frame.MVST(formula, data, "binomial")
NK <- frame.results$n
p <- frame.results$p
X <- frame.results$X
X.standardised <- frame.results$X.standardised
X.sd <- frame.results$X.sd
X.mean <- frame.results$X.mean
X.indicator <- frame.results$X.indicator 
offset <- frame.results$offset
Y <- frame.results$Y
N.all <- length(Y)
J <- ncol(Y)
which.miss <- frame.results$which.miss
n.miss <- N.all - sum(which.miss)

  
#### W matrix
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
W.quants <- common.Wcheckformat.leroux(W)
K <- W.quants$n
N <- NK / K
    if(ceiling(N)!= floor(N)) stop("The number of data points in Y divided by the number of rows in W is not a whole number.", call.=FALSE)
 

  
#### Create a missing list
    if(n.miss>0)
    {
    miss.locator <- array(NA, c(n.miss, 2))
    colnames(miss.locator) <- c("row", "column")
    locations <- which(which.miss==0)
    miss.locator[ ,1] <- ceiling(locations/J)
    miss.locator[ ,2] <- locations - (miss.locator[ ,1]-1) * J
    }else
    {
    miss.locator <- NA    
    }


#### Check and format the trials argument
    if(ncol(trials)!=J) stop("trials has the wrong number of columns.", call.=FALSE)
    if(nrow(trials)!=NK) stop("trials has the wrong number of rows.", call.=FALSE)
    if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values.", call.=FALSE)
int.check <- N.all-sum(ceiling(trials)==floor(trials))
    if(int.check > 0) stop("the numbers of trials has non-integer values.", call.=FALSE)
    if(min(trials)<=0) stop("the numbers of trials has zero or negative values.", call.=FALSE)
failures <- trials - Y
    if(sum(Y>trials, na.rm=TRUE)>0) stop("the response variable has larger values that the numbers of trials.", call.=FALSE)

  
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
    alpha <- c(runif(1), runif(1))
    fix.rho.T <- FALSE   
    }else
    {
    alpha <- rho.T
    fix.rho.T <- TRUE
    }
    if(!is.numeric(alpha)) stop("rho.T is fixed but is not numeric.", call.=FALSE)  
    if(length(alpha)!=2) stop("rho.T is fixed but is not of length 2.", call.=FALSE)  


#### Priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(100000, p)
    if(is.null(prior.Sigma.df)) prior.Sigma.df <- 2
    if(is.null(prior.Sigma.scale)) prior.Sigma.scale <- rep(100000, J)
prior.beta.check(prior.mean.beta, prior.var.beta, p)
    if(!is.numeric(prior.Sigma.scale)) stop("prior.Sigma.scale has non-numeric values.", call.=FALSE)    
    if(sum(is.na(prior.Sigma.scale))!=0) stop("prior.Sigma.scale has missing values.", call.=FALSE)   


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
   results <- binomial.MVCARar2MCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, rho=rho, alpha=alpha, fix.rho.S=fix.rho.S, fix.rho.T=fix.rho.T, K=K, N=N, NK=NK, J=J, N.all=N.all, p=p, miss.locator=miss.locator, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale, verbose=verbose, chain=1)
   }else if(n.chains > 1 & ceiling(n.chains)==floor(n.chains) & n.cores==1)
   {
   #### Multiple chains in  series
   results <- as.list(rep(NA, n.chains))
         for(i in 1:n.chains)
         {
         results[[i]] <- binomial.MVCARar2MCMC(Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, rho=rho, alpha=alpha, fix.rho.S=fix.rho.S, fix.rho.T=fix.rho.T, K=K, N=N, NK=NK, J=J, N.all=N.all, p=p, miss.locator=miss.locator, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df,  prior.Sigma.scale=prior.Sigma.scale, verbose=verbose, chain=i)
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
   results <- clusterCall(compclust, fun=binomial.MVCARar2MCMC, Y=Y, failures=failures, trials=trials, offset=offset, X.standardised=X.standardised, W=W, rho=rho, alpha=alpha, fix.rho.S=fix.rho.S, fix.rho.T=fix.rho.T, K=K, N=N, NK=NK, J=J, N.all=N.all, p=p, miss.locator=miss.locator, n.miss=n.miss, burnin=burnin, n.sample=n.sample, thin=thin, MALA=MALA, n.beta.block=n.beta.block, list.block=list.block, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.Sigma.df=prior.Sigma.df, prior.Sigma.scale=prior.Sigma.scale, verbose=verbose, chain="all")
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
    names(accept.final) <- c("beta", "phi", "rho.S", "rho.T", "Sigma")
    accept.final[1] <- 100 * sum(results$accept.beta[1:J]) / sum(results$accept.beta[(J+1):(2*J)])
    accept.final[2] <- 100 * results$accept[1] / results$accept[2]
        if(!fix.rho.S) accept.final[3] <- 100 * results$accept[3] / results$accept[4]
        if(!fix.rho.T) accept.final[4] <- 100
    accept.final[5] <- 100

    ## Compute the fitted deviance
    mean.beta <- matrix(apply(results$samples.beta, 2, mean), nrow=p, ncol=J, byrow=F)
    mean.phi <- matrix(apply(results$samples.phi, 2, mean), nrow=NK, ncol=J, byrow=T)
    mean.logit <- X.standardised %*% mean.beta + mean.phi + offset
    mean.prob <- exp(mean.logit)  / (1 + exp(mean.logit))
    fitted.mean <- trials * mean.prob
    deviance.fitted <- -2 * sum(dbinom(x=as.numeric(t(Y)), size=as.numeric(t(trials)), prob=as.numeric(t(mean.prob)), log=TRUE), na.rm=TRUE)
    modelfit <- common.modelfit(results$samples.loglike, deviance.fitted)

    ## Create the fitted values and residuals
    fitted.values <- matrix(apply(results$samples.fitted, 2, mean), nrow=NK, ncol=J, byrow=T)
    response.residuals <- Y - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
    residuals <- list(response=response.residuals, pearson=pearson.residuals)
 
    ## Transform the parameters back to the original covariate scale
    samples.beta.orig <- results$samples.beta
        for(r in 1:J)
        {
        samples.beta.orig[ ,((r-1)*p+1):(r*p)] <- common.betatransform(results$samples.beta[ ,((r-1)*p+1):(r*p) ], X.indicator, X.mean, X.sd, p, FALSE)
        }

    ## Create the samples object
        if(fix.rho.S & fix.rho.T)
        {
        samples.rhoext <- NA
        }else if(fix.rho.S & !fix.rho.T)
        {
        samples.rhoext <- results$samples.alpha
        colnames(samples.rhoext) <- c("rho1.T", "rho2.T")
        }else if(!fix.rho.S & fix.rho.T)
        {
        samples.rhoext <- results$samples.rho  
        names(samples.rhoext) <- "rho.S"
        }else
        {
        samples.rhoext <- cbind(results$samples.rho, results$samples.alpha)
        colnames(samples.rhoext) <- c("rho.S", "rho1.T", "rho2.T")
        }
    samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(results$samples.phi), Sigma=results$samples.Sigma, rho=mcmc(samples.rhoext), fitted=mcmc(results$samples.fitted), Y=mcmc(results$samples.Y))
    
    ## Create a summary object
    n.keep <- floor((n.sample - burnin)/thin)
    summary.beta <- t(rbind(apply(samples.beta.orig, 2, mean), apply(samples.beta.orig, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.final[names(accept.final)=="beta"],p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
    col.name <- rep(NA, p*(J-1))
        if(is.null(colnames(Y)))
        {
            for(r in 1:J)
            {
            col.name[((r-1)*p+1):(r*p)] <- paste("Variable ", r,  " - ", colnames(X), sep="")   
            }
        }else
        {
            for(r in 1:J)
            {
            col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[r],  " - ", colnames(X), sep="")   
            }
        }
    rownames(summary.beta) <- col.name
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")

    summary.hyper <- array(NA, c((J+3) ,7))
        for(r in 1:J)
        {
        summary.hyper[r, 1] <- mean(results$samples.Sigma[ ,r,r])
        summary.hyper[r, 2:3] <- quantile(results$samples.Sigma[ ,r,r], c(0.025, 0.975)) 
        summary.hyper[r, 4] <- n.keep
        summary.hyper[r, 5] <- 100
        summary.hyper[r, 6] <- effectiveSize(results$samples.Sigma[ ,r,r])
        summary.hyper[r, 7] <- geweke.diag(results$samples.Sigma[ ,r,r])$z    
        }
 
        if(!fix.rho.S)
        {
        summary.hyper[(J+1), 1:3] <- c(mean(results$samples.rho), quantile(results$samples.rho, c(0.025, 0.975)))
        summary.hyper[(J+1), 4:5] <- c(n.keep, accept.final[names(accept.final)=="rho.S"])
        summary.hyper[(J+1), 6:7] <- c(effectiveSize(results$samples.rho), geweke.diag(results$samples.rho)$z)
        }else
        {
        summary.hyper[(J+1), 1:3] <- c(rho, rho, rho)
        summary.hyper[(J+1), 4:5] <- rep(NA, 2)
        summary.hyper[(J+1), 6:7] <- rep(NA, 2)
        }
  
        if(!fix.rho.T)
        {
        summary.hyper[(J+2), 1:3] <- c(mean(results$samples.alpha[ ,1]), quantile(results$samples.alpha[ ,1], c(0.025, 0.975)))
        summary.hyper[(J+2), 4:5] <- c(n.keep, accept.final[names(accept.final)=="rho.T"])
        summary.hyper[(J+2), 6:7] <- c(effectiveSize(results$samples.alpha[ ,1]), geweke.diag(results$samples.alpha[ ,1])$z)
        summary.hyper[(J+3), 1:3] <- c(mean(results$samples.alpha[ ,2]), quantile(results$samples.alpha[ ,2], c(0.025, 0.975)))
        summary.hyper[(J+3), 4:5] <- c(n.keep, accept.final[names(accept.final)=="rho.T"])
        summary.hyper[(J+3), 6:7] <- c(effectiveSize(results$samples.alpha[ ,2]), geweke.diag(results$samples.alpha[ ,2])$z)
        }else
        {
        summary.hyper[(J+2), 1:3] <- c(alpha[1], alpha[1], alpha[1])
        summary.hyper[(J+2), 4:5] <- rep(NA, 2)
        summary.hyper[(J+2), 6:7] <- rep(NA, 2)
        summary.hyper[(J+3), 1:3] <- c(alpha[2], alpha[2], alpha[2])
        summary.hyper[(J+3), 4:5] <- rep(NA, 2)
        summary.hyper[(J+3), 6:7] <- rep(NA, 2)
        }

    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[((J*p)+1): nrow(summary.results)] <- c(paste(rep("Sigma",J), 1:J, 1:J, sep=""), "rho.S", "rho1.T", "rho2.T")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    }else
    {
    #### If n.chains > 1
    ## Compute the acceptance rates
    accept.final <- rep(NA, 5)
    names(accept.final) <- c("beta", "phi", "rho.S", "rho.T", "Sigma")
    accept.final[5] <- 100
    accept.temp <- lapply(results, function(l) l[["accept.beta"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.final[1] <- 100 * sum(accept.temp2[ ,1:J]) / sum(accept.temp2[ ,(J+1):(2*J)])
    accept.temp <- lapply(results, function(l) l[["accept"]])
    accept.temp2 <- do.call(what=rbind, args=accept.temp)
    accept.final[2] <- 100 * sum(accept.temp2[ ,1]) / sum(accept.temp2[ ,2])
        if(!fix.rho.S) accept.final[3] <- 100 * sum(accept.temp2[ ,3]) / sum(accept.temp2[ ,4])
        if(!fix.rho.T) accept.final[4] <- 100

    ## Extract the samples into separate matrix and list objects
    samples.beta.list <- lapply(results, function(l) l[["samples.beta"]])
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)   
    samples.phi.list <- lapply(results, function(l) l[["samples.phi"]])
    samples.phi.matrix <- do.call(what=rbind, args=samples.phi.list)
        if(!fix.rho.S)
        {
        samples.rho.list <- lapply(results, function(l) l[["samples.rho"]])
        samples.rho.matrix <- do.call(what=rbind, args=samples.rho.list)
        }
        if(!fix.rho.T)
        {
        samples.alpha.list <- lapply(results, function(l) l[["samples.alpha"]])
        samples.alpha.matrix <- do.call(what=rbind, args=samples.alpha.list)
        }
    samples.Sigma.list <- lapply(results, function(l) l[["samples.Sigma"]])
    samples.loglike.list <- lapply(results, function(l) l[["samples.loglike"]])
    samples.loglike.matrix <- do.call(what=rbind, args=samples.loglike.list)
    samples.fitted.list <- lapply(results, function(l) l[["samples.fitted"]])
    samples.fitted.matrix <- do.call(what=rbind, args=samples.fitted.list)
        if(n.miss>0) samples.Y.list <- lapply(results, function(l) l[["samples.Y"]])

    ## Compute the fitted deviance
    mean.beta <- matrix(apply(samples.beta.matrix, 2, mean), nrow=p, ncol=J, byrow=F)
    mean.phi <- matrix(apply(samples.phi.matrix, 2, mean), nrow=NK, ncol=J, byrow=T)
    mean.logit <- X.standardised %*% mean.beta + mean.phi + offset
    mean.prob <- exp(mean.logit)  / (1 + exp(mean.logit))
    fitted.mean <- trials * mean.prob
    deviance.fitted <- -2 * sum(dbinom(x=as.numeric(t(Y)), size=as.numeric(t(trials)), prob=as.numeric(t(mean.prob)), log=TRUE), na.rm=TRUE)
    modelfit <- common.modelfit(samples.loglike.matrix, deviance.fitted)

    ## Create the fitted values and residuals
    fitted.values <- matrix(apply(samples.fitted.matrix, 2, mean), nrow=NK, ncol=J, byrow=T)
    response.residuals <- Y - fitted.values
    pearson.residuals <- response.residuals /sqrt(fitted.values * (1 - mean.prob))
    residuals <- list(response=response.residuals, pearson=pearson.residuals)

    ## Transform the parameters back to the original covariate scale.
    samples.beta.list <- samples.beta.list
        for(j in 1:n.chains)
        {
            for(r in 1:J)
            {
            samples.beta.list[[j]][ ,((r-1)*p+1):(r*p)] <- common.betatransform(samples.beta.list[[j]][ ,((r-1)*p+1):(r*p) ], X.indicator, X.mean, X.sd, p, FALSE)
            }
        }
    samples.beta.matrix <- do.call(what=rbind, args=samples.beta.list)

    ## Create MCMC objects
    beta.mcmc <- mcmc.list(lapply(samples.beta.list, mcmc))
    phi.mcmc <- mcmc.list(lapply(samples.phi.list, mcmc))
    fitted.mcmc <- mcmc.list(lapply(samples.fitted.list, mcmc))
        if(n.miss>0) 
        {    
        Y.mcmc <- mcmc.list(lapply(samples.Y.list, mcmc))
        }else
        {
        Y.mcmc <- NA    
        }
        if(fix.rho.S & fix.rho.T)
        {
        rhoext.mcmc <- NA
        }else if(fix.rho.S & !fix.rho.T)
        {
            for(j in 1:n.chains)
            {
            colnames(samples.alpha.list[[j]]) <- c("rho1.T", "rho2.T")
            } 
        rhoext.mcmc <- mcmc.list(lapply(samples.alpha.list, mcmc))
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
            rho.temp[[j]] <- cbind(samples.rho.list[[j]], samples.alpha.list[[j]])
            colnames(rho.temp[[j]]) <- c("rho.S", "rho1.T", "rho2.T")
            }
        rhoext.mcmc <- mcmc.list(lapply(rho.temp, mcmc))
        }
    samples <- list(beta=beta.mcmc, phi=phi.mcmc, rho=rhoext.mcmc, Sigma=samples.Sigma.list, fitted=fitted.mcmc, Y=Y.mcmc)

    ## create a summary object
    n.keep <- floor((n.sample - burnin)/thin) * n.chains
    summary.beta <- t(rbind(apply(samples.beta.matrix, 2, mean), apply(samples.beta.matrix, 2, quantile, c(0.025, 0.975)))) 
    summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.final[names(accept.final)=="beta"],p), effectiveSize(beta.mcmc), gelman.diag(beta.mcmc)$psrf[ ,2])
    col.name <- rep(NA, p*(J-1))
        if(is.null(colnames(Y)))
        {
            for(r in 1:J)
            {
            col.name[((r-1)*p+1):(r*p)] <- paste("Variable ", r,  " - ", colnames(X), sep="")   
            }
        }else
        {
            for(r in 1:J)
            {
            col.name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[r],  " - ", colnames(X), sep="")   
            }
        }
    rownames(summary.beta) <- col.name
    colnames(summary.beta) <- c("Mean", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "PSRF (upper 95% CI)")

    summary.hyper <- array(NA, c((J+3) ,7))
        for(r in 1:J)
        {
        temp <- NA
        temp2 <- as.list(rep(NA, n.chains))
            for(v in 1:n.chains)
            {
            temp <- c(temp, samples.Sigma.list[[v]][ ,r,r])
            temp2[[v]] <- mcmc(samples.Sigma.list[[v]][ ,r,r])
            }
        temp <- temp[-1]    
        summary.hyper[r, 1] <- mean(temp)
        summary.hyper[r, 2:3] <- quantile(temp, c(0.025, 0.975)) 
        summary.hyper[r, 4] <- n.keep
        summary.hyper[r, 5] <- 100
        summary.hyper[r, 6] <- effectiveSize(mcmc.list(temp2))
        summary.hyper[r, 7] <- gelman.diag(mcmc.list(temp2))$psrf[ ,2]
        }
 
        if(!fix.rho.S)
        {
        temp <- mcmc.list(lapply(samples.rho.list, mcmc))
        summary.hyper[(J+1), 1:3] <- c(mean(samples.rho.matrix), quantile(samples.rho.matrix, c(0.025, 0.975)))
        summary.hyper[(J+1), 4:5] <- c(n.keep, accept.final[names(accept.final)=="rho.S"])
        summary.hyper[(J+1), 6:7] <- c(effectiveSize(temp), gelman.diag(temp)$psrf[ ,2])
        }else
        {
        summary.hyper[(J+1), 1:3] <- c(rho, rho, rho)
        summary.hyper[(J+1), 4:5] <- rep(NA, 2)
        summary.hyper[(J+1), 6:7] <- rep(NA, 2)
        }
  
        if(!fix.rho.T)
        {
        temp <- mcmc.list(lapply(samples.alpha.list, mcmc))
        summary.hyper[(J+2), 1:3] <- c(mean(samples.alpha.matrix[ ,1]), quantile(samples.alpha.matrix[ ,1], c(0.025, 0.975)))
        summary.hyper[(J+2), 4:5] <- c(n.keep, accept.final[names(accept.final)=="rho.T"])
        summary.hyper[(J+2), 6:7] <- c(effectiveSize(temp)[1], gelman.diag(temp)$psrf[ ,2][1])
        summary.hyper[(J+3), 1:3] <- c(mean(samples.alpha.matrix[ ,2]), quantile(samples.alpha.matrix[ ,2], c(0.025, 0.975)))
        summary.hyper[(J+3), 4:5] <- c(n.keep, accept.final[names(accept.final)=="rho.T"])
        summary.hyper[(J+3), 6:7] <- c(effectiveSize(temp)[2], gelman.diag(temp)$psrf[ ,2][2])
        }else
        {
        summary.hyper[(J+2), 1:3] <- c(alpha[1], alpha[1], alpha[1])
        summary.hyper[(J+2), 4:5] <- rep(NA, 2)
        summary.hyper[(J+2), 6:7] <- rep(NA, 2)
        summary.hyper[(J+3), 1:3] <- c(alpha[2], alpha[2], alpha[2])
        summary.hyper[(J+3), 4:5] <- rep(NA, 2)
        summary.hyper[(J+3), 6:7] <- rep(NA, 2)
        }
    summary.results <- rbind(summary.beta, summary.hyper)
    rownames(summary.results)[((J*p)+1): nrow(summary.results)] <- c(paste(rep("Sigma",J), 1:J, 1:J, sep=""), "rho.S", "rho1.T", "rho2.T")
    summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
    summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    }
    


###################################
#### Compile and return the results
###################################
model.string <- c("Likelihood model - Binomial (logit link function)", "\nRandom effects model - Multivariate Autoregressive order 2 CAR model\n")
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
