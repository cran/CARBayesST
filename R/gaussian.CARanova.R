gaussian.CARanova <- function(formula, data=NULL, W, burnin, n.sample, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, verbose=TRUE)
{
#### Check on the verbose option
    if(is.null(verbose)) verbose=TRUE     
    if(!is.logical(verbose)) stop("the verbose option is not logical.", call.=FALSE)
    
    if(verbose)
    {
        cat("Setting up the model\n")
        a<-proc.time()
    }else{}
    
    
    
##############################################
#### Format the arguments and check for errors
##############################################
#### Overall formula object
frame <- try(suppressWarnings(model.frame(formula, data=data, na.action=na.pass)), silent=TRUE)
    if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths or the data object has not been specified.", call.=FALSE)
    
    
#### Design matrix
## Create the matrix
X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
    if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
    if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
    
N.all <- nrow(X)
p <- ncol(X)
    
    ## Check for linearly related columns
    cor.X <- suppressWarnings(cor(X))
    diag(cor.X) <- 0
    
    if(max(cor.X, na.rm=TRUE)==1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    if(min(cor.X, na.rm=TRUE)==-1) stop("the covariate matrix has two exactly linearly related columns.", call.=FALSE)
    
    if(p>1)
    {
        if(sort(apply(X, 2, sd))[2]==0) stop("the covariate matrix has two intercept terms.", call.=FALSE)
    }else
    {
    }
    
    ## Standardise the matrix
    X.standardised <- X
    X.sd <- apply(X, 2, sd)
    X.mean <- apply(X, 2, mean)
    X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back
    
    for(j in 1:p)
    {
        if(length(table(X[ ,j]))>2)
        {
            X.indicator[j] <- 1
            X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
        }else if(length(table(X[ ,j]))==1)
        {
            X.indicator[j] <- 2
        }else
        {
            X.indicator[j] <- 0
        }
    }
    
    
#### Response variable
Y <- model.response(frame)
    if(sum(is.na(Y))>0) stop("the response has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(Y)) stop("the response variable has non-numeric values.", call.=FALSE)
    
    
    
#### Offset variable
offset <- try(model.offset(frame), silent=TRUE)
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,N.all)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
    
    
    
#### Format and check the neighbourhood matrix W
    if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
K <- nrow(W)
    if(ncol(W)!= K) stop("W has the wrong number of columns.", call.=FALSE)
    if(nrow(W)!= K) stop("W has the wrong number of rows.", call.=FALSE)
    if(floor(N.all/K)!=ceiling(N.all/K)) stop("The number of spatial areas is not a multiple of the number of data points.", call.=FALSE)
N <- N.all / K
    if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
    if(min(W)<0) stop("W has negative elements.", call.=FALSE)
    if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
    if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    

        
    
#### Specify the initial parameter values
beta <- lm(Y~X.standardised-1, offset=offset)$coefficients
regression.vec <- X.standardised %*% beta
res.temp <- Y - regression.vec - offset
phi <- rnorm(n=K, mean=0, sd = 2*sd(res.temp))
delta <- rnorm(n=N, mean=0, sd = 2*sd(res.temp))
tau2.phi <- runif(1, min=var(res.temp)/2, max=var(res.temp)*2)
tau2.delta <- runif(1, min=var(res.temp)/2, max=var(res.temp)*2)
nu2 <- runif(1, min=var(res.temp)/2, max=var(res.temp)*2)
rho <- runif(1)
lambda <- runif(1)  
    
    
#### Check and specify the priors
## Put in default priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
    if(is.null(prior.nu2)) prior.nu2 <- c(0.001, 0.001)

    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)       
    
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)
    
    if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
    
    if(length(prior.nu2)!=2) stop("the prior value for nu2 is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.nu2)) stop("the prior value for nu2 is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.nu2))!=0) stop("the prior value for nu2 has missing values.", call.=FALSE)       



#### MCMC quantities
## Checks
    if(is.null(burnin)) stop("the burnin argument is missing", call.=FALSE)
    if(is.null(n.sample)) stop("the n.sample argument is missing", call.=FALSE)
    if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
    if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE) 
    if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
    if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
    if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
    if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)
    if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
    if(n.sample <= thin)  stop("thin is greater than n.sample.", call.=FALSE)
    if(burnin!=round(burnin)) stop("burnin is not an integer.", call.=FALSE) 
    if(n.sample!=round(n.sample)) stop("n.sample is not an integer.", call.=FALSE) 
    if(thin!=round(thin)) stop("thin is not an integer.", call.=FALSE) 
    
    
    
#### Set up matrices to store samples
n.keep <- floor((n.sample - burnin)/thin)
samples.beta <- array(NA, c(n.keep, p))
samples.phi <- array(NA, c(n.keep, K))
samples.delta <- array(NA, c(n.keep, N))
samples.rho <- array(NA, c(n.keep, 1))
samples.lambda <- array(NA, c(n.keep, 1))
samples.nu2 <- array(NA, c(n.keep, 1))
samples.fitted <- array(NA, c(n.keep, N.all))
samples.deviance <- array(NA, c(n.keep, 1))
samples.tau2 <- array(NA, c(n.keep, 2))
colnames(samples.tau2) <- c("tau2.phi", "tau2.delta")    
    
    
    
#### Specify the Metropolis quantities
accept.all <- rep(0,4)
accept <- accept.all
proposal.sd.rho <- 0.02
proposal.sd.lambda <- 0.02
tau2.phi.shape <- prior.tau2[1] + K/2
tau2.delta.shape <- prior.tau2[1] + N/2
nu2.shape <- prior.nu2[1] + N*K/2    
    
    
    
#### Spatial quantities
## Create the triplet object
W.triplet <- c(NA, NA, NA)
    for(i in 1:K)
    {
        for(j in 1:K)
        {
            if(W[i,j]>0)
            {
                W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
            }else{}
        }
    }
    W.triplet <- W.triplet[-1, ]     
    W.n.triplet <- nrow(W.triplet) 
    W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
    W.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
    
    
    ## Create the start and finish points for W updating
    W.begfin <- array(NA, c(K, 2))     
    temp <- 1
    for(i in 1:K)
    {
        W.begfin[i, ] <- c(temp, (temp + W.neighbours[i]-1))
        temp <- temp + W.neighbours[i]
    }
    
    
    ## Create the determinant     
    Wstar <- diag(apply(W,1,sum)) - W
    Wstar.eigen <- eigen(Wstar)
    Wstar.val <- Wstar.eigen$values
    det.Q.W <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))    
    
    
    
    
    #### Temporal quantities
    ## Temporal neighbourhood matrix
    D <-array(0, c(N,N))
    for(i in 1:N)
    {
        for(j in 1:N)
        {
            if(abs((i-j))==1)  D[i,j] <- 1 
        }    
    }
    
    
    ## Create the triplet object
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
    
    
    ## Create the start and finish points for W updating
    D.begfin <- array(NA, c(N, 2))     
    temp <- 1
    for(i in 1:N)
    {
        D.begfin[i, ] <- c(temp, (temp + D.neighbours[i]-1))
        temp <- temp + D.neighbours[i]
    }
    
    
    ## Create the determinant     
    Dstar <- diag(apply(D,1,sum)) - D
    Dstar.eigen <- eigen(Dstar)
    Dstar.val <- Dstar.eigen$values
    det.Q.D <-  0.5 * sum(log((lambda * Dstar.val + (1-lambda))))    
    
    
#### Specify quantities that do not change
offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE)
Y.mat.trans <- t(Y.mat)
phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)
delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)
    


#### Beta update quantities
data.precision.beta <- t(X.standardised) %*% X.standardised
    if(length(prior.var.beta)==1)
    {
    prior.precision.beta <- 1 / prior.var.beta
    }else
    {
    prior.precision.beta <- solve(diag(prior.var.beta))
    }




    ###########################
    #### Run the Bayesian model
    ###########################
    ## Start timer
    if(verbose)
    {
        cat("Collecting", n.sample, "samples\n", sep = " ")
        progressBar <- txtProgressBar(style = 3)
        percentage.points<-round((1:100/100)*n.sample)
    }else
    {
        percentage.points<-round((1:100/100)*n.sample)     
    }
    
    
    for(j in 1:n.sample)
    {
        ##################
        ## Sample from nu2
        ##################
        nu2.offset <- as.numeric(Y.mat - offset.mat - regression.mat - phi.mat - delta.mat)
        nu2.scale <- prior.nu2[2]  + sum(nu2.offset^2)/2
        nu2 <- 1 / rgamma(1, nu2.shape, scale=(1/nu2.scale)) 
        
        
        
        ####################
        ## Sample from beta
        ####################
        fc.precision <- prior.precision.beta + data.precision.beta / nu2
        fc.var <- solve(fc.precision)
        beta.offset <- as.numeric(Y.mat - offset.mat - phi.mat - delta.mat)
        beta.offset2 <- t(X.standardised) %*% beta.offset / nu2 + prior.precision.beta %*% prior.mean.beta
        fc.mean <- fc.var %*% beta.offset2
        chol.var <- t(chol(fc.var))
        beta <- fc.mean + chol.var %*% rnorm(p)        
        regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)  
        
        
        ####################
        ## Sample from phi
        ####################
        phi.offset <- Y.mat - offset.mat - regression.mat -  delta.mat
        phi.offset2 <- apply(phi.offset,1, sum)
        temp1 <- gaussiancarupdate(W.triplet, W.begfin, W.triplet.sum, K, phi, tau2.phi, nu2, phi.offset2, rho, N)
        phi <- temp1
        phi <- phi - mean(phi)
        phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)    
        
        
        
        ####################
        ## Sample from delta
        ####################
        delta.offset <- Y.mat - offset.mat - regression.mat -  phi.mat
        delta.offset2 <- apply(delta.offset,2, sum)
        temp2 <- gaussiancarupdate(D.triplet, D.begfin, D.triplet.sum, N, delta, tau2.delta, nu2, delta.offset2, lambda, K)
        delta <- temp2
        delta <- delta - mean(delta)
        delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)

                
        
        #######################
        ## Sample from tau2.phi
        #######################
        temp2.phi <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, rho)
        tau2.phi.scale <- temp2.phi + prior.tau2[2] 
        tau2.phi <- 1 / rgamma(1, tau2.phi.shape, scale=(1/tau2.phi.scale))
        
        
        #########################
        ## Sample from tau2.delta
        #########################
        temp2.delta <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, lambda)
        tau2.delta.scale <- temp2.delta + prior.tau2[2] 
        tau2.delta <- 1 / rgamma(1, tau2.delta.shape, scale=(1/tau2.delta.scale))
        
        
        
        ##################
        ## Sample from rho
        ##################
        proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=proposal.sd.rho)   
        temp3 <- quadform(W.triplet, W.triplet.sum, W.n.triplet, K, phi, phi, proposal.rho)
        det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))              
        logprob.current <- det.Q.W - temp2.phi / tau2.phi
        logprob.proposal <- det.Q.proposal - temp3 / tau2.phi
        prob <- exp(logprob.proposal - logprob.current)
        
        #### Accept or reject the proposal
        if(prob > runif(1))
        {
            rho <- proposal.rho
            det.Q.W <- det.Q.proposal
            accept[1] <- accept[1] + 1           
        }else
        {
        }              
        accept[2] <- accept[2] + 1           
        
        
        
        #####################
        ## Sample from lambda
        #####################
        proposal.lambda <- rtrunc(n=1, spec="norm", a=0, b=1, mean=lambda, sd=proposal.sd.lambda)   
        temp3 <- quadform(D.triplet, D.triplet.sum, D.n.triplet, N, delta, delta, proposal.lambda)
        det.Q.proposal <- 0.5 * sum(log((proposal.lambda * Dstar.val + (1-proposal.lambda))))              
        logprob.current <- det.Q.D - temp2.delta / tau2.delta
        logprob.proposal <- det.Q.proposal - temp3 / tau2.delta
        prob <- exp(logprob.proposal - logprob.current)
        
        #### Accept or reject the proposal
        if(prob > runif(1))
        {
            lambda <- proposal.lambda
            det.Q.D <- det.Q.proposal
            accept[3] <- accept[3] + 1           
        }else
        {
        }              
        accept[4] <- accept[4] + 1           
        
        
        
        #########################
        ## Calculate the deviance
        #########################
        fitted <- as.numeric(offset.mat + regression.mat + phi.mat  + delta.mat)
        deviance.all <- dnorm(Y, mean = fitted, sd = rep(sqrt(nu2),N.all), log=TRUE)
        deviance <- -2 * sum(deviance.all)  
        
        
        ###################
        ## Save the results
        ###################
        if(j > burnin & (j-burnin)%%thin==0)
        {
            ele <- (j - burnin) / thin
            samples.beta[ele, ] <- beta
            samples.phi[ele, ] <- phi
            samples.delta[ele, ] <- delta
            samples.rho[ele, ] <- rho
            samples.nu2[ele, ] <- nu2
            samples.lambda[ele, ] <- lambda
            samples.deviance[ele, ] <- deviance
            samples.fitted[ele, ] <- fitted
            samples.tau2[ele, ] <- c(tau2.phi, tau2.delta)
        }else
        {
        }
        
        
        ########################################
        ## Self tune the acceptance probabilties
        ########################################
        k <- j/100
        if(ceiling(k)==floor(k))
        {
            #### Determine the acceptance probabilities
            accept.rho <- 100 * accept[1] / accept[2]
            accept.lambda <- 100 * accept[3] / accept[4]
            accept.all <- accept.all + accept
            accept <- rep(0,4)
            
            #### rho tuning parameter
            if(accept.rho > 50)
            {
                proposal.sd.rho <- min(2 * proposal.sd.rho, 0.5)
            }else if(accept.rho < 40)              
            {
                proposal.sd.rho <- 0.5 * proposal.sd.rho
            }else
            {
            }
            #### lambda tuning parameter
            if(accept.lambda > 50)
            {
                proposal.sd.lambda <- min(2 * proposal.sd.lambda, 0.5)
            }else if(accept.lambda < 40)              
            {
                proposal.sd.lambda <- 0.5 * proposal.sd.lambda
            }else
            {
            }      
        }else
        {   
        }
        
        
        ################################       
        ## print progress to the console
        ################################
        if(j %in% percentage.points & verbose)
        {
            setTxtProgressBar(progressBar, j/n.sample)
        }
    }
    
    # end timer
    if(verbose)
    {
        cat("\nSummarising results")
        close(progressBar)
    }else
    {}
    
    
###################################
#### Summarise and save the results 
###################################
## Compute the acceptance rates
accept.rho <- 100 * accept.all[1] / accept.all[2]
accept.lambda <- 100 * accept.all[3] / accept.all[4]
accept.final <- c(rep(100,3), accept.rho, accept.lambda)
names(accept.final) <- c("beta", "phi", "delta", "rho", "lambda")        
    
   

## Compute DIC
median.phi <- apply(samples.phi, 2, median)
median.delta <- apply(samples.delta, 2, median)  
median.phi.mat <- matrix(rep(median.phi, N), byrow=F, nrow=K)
median.delta.mat <- matrix(rep(median.delta, K), byrow=T, nrow=K)
median.beta <- apply(samples.beta,2,median)
regression.mat <- matrix(X.standardised %*% median.beta, nrow=K, ncol=N, byrow=FALSE)   
fitted.median <- as.numeric(offset.mat + regression.mat + median.phi.mat + median.delta.mat)    
nu2.median <- median(samples.nu2)
deviance.fitted <- -2 * sum(dnorm(Y, mean = fitted.median, sd = rep(sqrt(nu2.median),N.all), log = TRUE))
p.d <- median(samples.deviance) - deviance.fitted
DIC <- 2 * median(samples.deviance) - deviance.fitted    
 
    
#### Compute the Conditional Predictive Ordinate  
CPO <- rep(NA, N.all)
    for(j in 1:N.all)
    {
    CPO[j] <- 1/median((1 / dnorm(Y[j], mean=samples.fitted[ ,j], sd=sqrt(samples.nu2))))    
    }
LMPL <- sum(log(CPO))  

    
    
## Create the Fitted values
fitted.values <- apply(samples.fitted, 2, median)
residuals <- as.numeric(Y) - fitted.values
    
    
#### transform the parameters back to the origianl covariate scale.
samples.beta.orig <- samples.beta
number.cts <- sum(X.indicator==1)     
    if(number.cts>0)
    {
        for(r in 1:p)
        {
            if(X.indicator[r]==1)
            {
                samples.beta.orig[ ,r] <- samples.beta[ ,r] / X.sd[r]
            }else if(X.indicator[r]==2 & p>1)
            {
                X.transformed <- which(X.indicator==1)
                samples.temp <- as.matrix(samples.beta[ ,X.transformed])
                for(s in 1:length(X.transformed))
                {
                    samples.temp[ ,s] <- samples.temp[ ,s] * X.mean[X.transformed[s]]  / X.sd[X.transformed[s]]
                }
                intercept.adjustment <- apply(samples.temp, 1,sum) 
                samples.beta.orig[ ,r] <- samples.beta[ ,r] - intercept.adjustment
            }else
            {
            }
        }
    }else
    {
    }
    
    
#### Create a summary object
samples.beta.orig <- mcmc(samples.beta.orig)
summary.beta <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(100,p), effectiveSize(samples.beta.orig), geweke.diag(samples.beta.orig)$z)
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept", "n.effective", "Geweke.diag")

summary.hyper <- array(NA, c(5, 7))     
summary.hyper[1,1:3] <- quantile(samples.tau2[ ,1], c(0.5, 0.025, 0.975))
summary.hyper[2,1:3] <- quantile(samples.tau2[ ,2], c(0.5, 0.025, 0.975))
summary.hyper[3,1:3] <- quantile(samples.nu2, c(0.5, 0.025, 0.975))
summary.hyper[4,1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
summary.hyper[5,1:3] <- quantile(samples.lambda, c(0.5, 0.025, 0.975))
rownames(summary.hyper) <- c("tau2.phi", "tau2.delta",  "nu2", "rho.phi", "rho.delta")     
summary.hyper[1, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,1])), geweke.diag(mcmc(samples.tau2[ ,1]))$z)     
summary.hyper[2, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.tau2[ ,2])), geweke.diag(mcmc(samples.tau2[ ,2]))$z)   
summary.hyper[3, 4:7] <- c(n.keep, 100, effectiveSize(mcmc(samples.nu2)), geweke.diag(mcmc(samples.nu2))$z)  
summary.hyper[4, 4:7] <- c(n.keep, accept.rho, effectiveSize(mcmc(samples.rho)), geweke.diag(mcmc(samples.rho))$z)   
summary.hyper[5, 4:7] <- c(n.keep, accept.lambda, effectiveSize(mcmc(samples.lambda)), geweke.diag(mcmc(samples.lambda))$z)     

summary.results <- rbind(summary.beta, summary.hyper)
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:7] <- round(summary.results[ , 4:7], 1)
    
    
## Compile and return the results
modelfit <- c(DIC, p.d, LMPL)
names(modelfit) <- c("DIC", "p.d", "LMPL")
samples.rhoext <- cbind(samples.rho, samples.lambda)
samples <- list(beta=mcmc(samples.beta.orig), phi=mcmc(samples.phi),  delta=mcmc(samples.delta), tau2=mcmc(samples.tau2), nu2=mcmc(samples.nu2), rho=mcmc(samples.rhoext), fitted=mcmc(samples.fitted))        
model.string <- c("Likelihood model - Gaussian (identity link function)", "\nLatent structure model - spatial and temporal main effects\n")
results <- list(summary.results=summary.results, samples=samples, fitted.values=fitted.values, residuals=residuals, modelfit=modelfit, accept=accept.final, localised.structure=NULL, formula=formula, model=model.string,  X=X)
class(results) <- "carbayesST"
    if(verbose)
    {
        b<-proc.time()
        cat(" finished in ", round(b[3]-a[3], 1), "seconds")
    }else
    {}
    return(results)
}
