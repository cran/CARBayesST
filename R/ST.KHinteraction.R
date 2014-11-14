ST.KHinteraction <- function(formula, data=NULL, W, burnin=0, n.sample=1000, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, verbose=TRUE)
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
int.check <- N.all-sum(ceiling(Y)==floor(Y))
    if(int.check > 0) stop("the respons variable has non-integer values.", call.=FALSE)
    if(min(Y)<0) stop("the response variable has negative values.", call.=FALSE)
log.Y <- log(Y)
log.Y[Y==0] <- -0.1     

     
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
  if(floor(N.all/K)!=ceiling(N.all/K)) stop("The number of spatial areas is not a multiple of the number of data points.", call.=FALSE)
N <- N.all / K
  if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
  if(min(W)<0) stop("W has negative elements.", call.=FALSE)
  if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)

  
  
  #### Format and check the MCMC quantities
  if(!is.numeric(burnin)) stop("burn-in is not a number", call.=FALSE)
  if(!is.numeric(n.sample)) stop("n.sample is not a number", call.=FALSE) 
  if(!is.numeric(thin)) stop("thin is not a number", call.=FALSE)
  if(n.sample <= 0) stop("n.sample is less than or equal to zero.", call.=FALSE)
  if(burnin < 0) stop("burn-in is less than zero.", call.=FALSE)
  if(thin <= 0) stop("thin is less than or equal to zero.", call.=FALSE)
  if(n.sample <= burnin)  stop("Burn-in is greater than n.sample.", call.=FALSE)
  if(burnin!=round(burnin)) stop("burnin is not an integer.", call.=FALSE) 
  if(n.sample!=round(n.sample)) stop("n.sample is not an integer.", call.=FALSE) 
  if(thin!=round(thin)) stop("thin is not an integer.", call.=FALSE) 
  
  
  #### Check and specify the priors
    if(is.null(prior.mean.beta)) prior.mean.beta <- rep(0, p)
    if(length(prior.mean.beta)!=p) stop("the vector of prior means for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.mean.beta)) stop("the vector of prior means for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.mean.beta))!=0) stop("the vector of prior means for beta has missing values.", call.=FALSE)       
     
    if(is.null(prior.var.beta)) prior.var.beta <- rep(1000, p)
    if(length(prior.var.beta)!=p) stop("the vector of prior variances for beta is the wrong length.", call.=FALSE)    
    if(!is.numeric(prior.var.beta)) stop("the vector of prior variances for beta is not numeric.", call.=FALSE)    
    if(sum(is.na(prior.var.beta))!=0) stop("the vector of prior variances for beta has missing values.", call.=FALSE)    
    if(min(prior.var.beta) <=0) stop("the vector of prior variances has elements less than zero", call.=FALSE)

  if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
  if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
  

     
     
  #### Specify the initial parameter values
  beta <- glm(Y~X.standardised-1, offset=offset, family=poisson)$coefficients
  phi <- rnorm(K)
  theta <- rnorm(K)
  alpha <- rnorm(N)
  delta <- rnorm(N)
  gamma <- rnorm(n=(N*K))
  tau2.phi <- runif(1)
  tau2.theta <- runif(1)
  tau2.alpha <- runif(1)
  tau2.delta <- runif(1)
  tau2.gamma <- runif(1)  
  
## Compute the blocking structure for beta     
blocksize.beta <- 5
     if(blocksize.beta >= p)
     {
     n.beta.block <- 1
     beta.beg <- 1
     beta.fin <- p
     }else
     {
     n.standard <- 1 + floor((p-blocksize.beta) / blocksize.beta)
     remainder <- p - n.standard * blocksize.beta
     
          if(remainder==0)
          {
          beta.beg <- c(1,seq((blocksize.beta+1), p, blocksize.beta))
          beta.fin <- seq(blocksize.beta, p, blocksize.beta)
          n.beta.block <- length(beta.beg)
          }else
          {
          beta.beg <- c(1, seq((blocksize.beta+1), p, blocksize.beta))
          beta.fin <- c(seq((blocksize.beta), p, blocksize.beta), p)
          n.beta.block <- length(beta.beg)
          }
     }         

     
#### Set up matrices to store samples
  n.keep <- floor((n.sample - burnin)/thin)
  samples.beta <- array(NA, c(n.keep, p))
  samples.space <- array(NA, c(n.keep, K))
  samples.time <- array(NA, c(n.keep, N))
  samples.gamma <- array(NA, c(n.keep, N.all))
  samples.tau2 <- array(NA, c(n.keep, 5))
  colnames(samples.tau2) <- c("tau2.phi", "tau2.theta", "tau2.alpha", "tau2.beta", "tau2.gamma")
  samples.fitted <- array(NA, c(n.keep, N.all))
  samples.deviance <- array(NA, c(n.keep, 1))

  
  
  #### Specify the Metropolis quantities
  accept.all <- rep(0,12)
  accept <- accept.all
  proposal.sd.phi <- 0.1
  proposal.sd.theta <- 0.1
  proposal.sd.alpha <- 0.1
  proposal.sd.beta <- 0.01
  proposal.sd.gamma <- 0.1
  proposal.sd.delta <- 0.1
  proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
  chol.proposal.corr.beta <- chol(proposal.corr.beta)   
  tau2.phi.shape <- prior.tau2[1] + (K-1)/2
  tau2.theta.shape <- prior.tau2[1] + K/2
  tau2.alpha.shape <- prior.tau2[1] + N/2
  tau2.delta.shape <- prior.tau2[1] + (N-1)/2
  tau2.gamma.shape <- prior.tau2[1] + N*K/2
  
  n.neighbours.space <- as.numeric(apply(W, 1, sum))
  space.duplet <- c(NA, NA)
  for(i in 1:K)
  {
    for(j in 1:K)
    {
      if(W[i,j]==1)
      {
        space.duplet <- rbind(space.duplet, c(i,j))     
      }else{}
    }
  }
  space.duplet <- space.duplet[-1, ]     
  n.space.duplet <- nrow(space.duplet) 
  
  
  #### Create the list form for the W matrix
  spacelist <- as.list(rep(NA,K))     
    for(i in 1:K)
    {
    spacelist[[i]] <- which(W[i, ]==1)     
    }

  #### Create the list form for the temporal neighbours
  timelist <- as.list(rep(NA,N))     
    for(i in 1:N)
    {
    timelist[[i]] <- c(i-1, i+1)     
    }  
timelist[[1]]  <- 2
timelist[[N]] <- N-1     
  
  
  #### Specify quantities that do not change
  offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
  regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
  Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE) 
  Y.area <- apply(Y.mat, 1,sum)
  Y.time <- apply(Y.mat, 2,sum)  
  Y.vec <- as.numeric(Y)
  Y.sum <- sum(Y)
  theta.mat <- matrix(rep(theta, N), byrow=F, nrow=K)
  phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)
  alpha.mat <- matrix(rep(alpha, K), byrow=T, nrow=K)
  delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)
  gamma.mat <- matrix(gamma, byrow=F, nrow=K)
  
  
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
    ####################
    ## Sample from phi
    ####################
    phi.offset <- exp(offset.mat + regression.mat + theta.mat + alpha.mat + delta.mat + gamma.mat)
    phi.offset.sum <- apply(phi.offset,1,sum)
    temp1 <- poissoncarupdate(W_list=spacelist, nsites=K, phi=phi, tau2=tau2.phi, y=Y.area, phi_tune=proposal.sd.phi, rho_num=1, rho_den=1, offset=phi.offset.sum)
    phi <- temp1[[1]] - mean(temp1[[1]])
    phi.mat <- matrix(rep(phi, N), byrow=F, nrow=K)
    accept[1] <- accept[1] + temp1[[2]]
    accept[2] <- accept[2] + K    
    
    
    ####################
    ## Sample from theta
    ####################
    theta.offset <- exp(offset.mat + regression.mat + phi.mat + alpha.mat + delta.mat + gamma.mat)
    theta.offset.sum <- apply(theta.offset,1,sum)
    temp2 <- poissonindepupdate(nsites=K, theta=theta, tau2=tau2.theta, y=Y.area, theta_tune=proposal.sd.theta, offset=theta.offset.sum)
    theta <- temp2[[1]] - mean(temp2[[1]])
    theta.mat <- matrix(rep(theta, N), byrow=F, nrow=K)
    accept[3] <- accept[3] + temp2[[2]]
    accept[4] <- accept[4] + K    
    
    
    ####################
    ## Sample from alpha
    ####################
    alpha.offset <- exp(offset.mat + regression.mat + phi.mat + theta.mat + delta.mat + gamma.mat)
    alpha.offset.sum <- apply(alpha.offset,2,sum)
    temp3 <- poissonindepupdate(nsites=N, theta=alpha, tau2=tau2.alpha, y=Y.time, theta_tune=proposal.sd.alpha, offset=alpha.offset.sum)
    alpha <- temp3[[1]] - mean(temp3[[1]])
    alpha.mat <- matrix(rep(alpha, K), byrow=T, nrow=K)
    accept[5] <- accept[5] + temp3[[2]]
    accept[6] <- accept[6] + N    
    
    
    ####################
    ## Sample from delta
    ####################
    delta.offset <- exp(offset.mat + regression.mat + phi.mat + theta.mat + alpha.mat + gamma.mat)
    delta.offset.sum <- apply(delta.offset,2,sum)
    temp4 <- poissoncarupdate(W_list=timelist, nsites=N, phi=delta, tau2=tau2.delta, y=Y.time, phi_tune=proposal.sd.delta, rho_num=1, rho_den=1, offset=delta.offset.sum)
    delta <- temp4[[1]] - mean(temp4[[1]])
    delta.mat <- matrix(rep(delta, K), byrow=T, nrow=K)
    accept[7] <- accept[7] + temp4[[2]]
    accept[8] <- accept[8] + N      
    
    
    ####################
    ## Sample from gamma
    ####################
    gamma.offset <- exp(offset.mat + regression.mat + phi.mat + theta.mat + alpha.mat + delta.mat)
    gamma.offset.vec <- as.numeric(gamma.offset)
    temp5 <- poissonindepupdate(nsites=N.all, theta=gamma, tau2=tau2.gamma, y=Y.vec, theta_tune=proposal.sd.gamma, offset=gamma.offset.vec)
    gamma <- temp5[[1]] - mean(temp5[[1]])
    gamma.mat <- matrix(gamma, byrow=F, nrow=K)
    accept[9] <- accept[9] + temp5[[2]]
    accept[10] <- accept[10] + N * K  
    
    
    ####################
    ## Sample from beta
    ####################
    proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)  
    proposal.beta <- beta
    offset.temp <- offset + as.numeric(phi.mat) + as.numeric(theta.mat) + as.numeric(alpha.mat) + as.numeric(delta.mat) + as.numeric(gamma.mat)       
       
       for(r in 1:n.beta.block)
       {
       proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
       prob <- poissonbetaupdate(X.standardised, N.all, p, beta, proposal.beta, offset.temp, Y, prior.mean.beta, prior.var.beta)
            if(prob > runif(1))
            {
            beta[beta.beg[r]:beta.fin[r]] <- proposal.beta[beta.beg[r]:beta.fin[r]]
            accept[11] <- accept[11] + 1  
            }else
            {
            proposal.beta[beta.beg[r]:beta.fin[r]] <- beta[beta.beg[r]:beta.fin[r]]
            }
        }

    accept[12] <- accept[12] + n.beta.block    
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   


       
       
    #######################
    ## Sample from tau2.phi
    #######################
    tau2.phi.scale <- prior.tau2[2]  + quadform(W_duplet1=space.duplet[ ,1], W_duplet2=space.duplet[ ,2], n_duplet=n.space.duplet,  nsites=K, phi=phi, theta=phi, nneighbours=n.neighbours.space, diagonal=1, offdiagonal=1)      
    tau2.phi <- 1 / rgamma(1, tau2.phi.shape, scale=(1/tau2.phi.scale)) 
    
    
    #########################
    ## Sample from tau2.theta
    #########################
    tau2.theta.scale <- prior.tau2[2]  + sum(theta^2)/2
    tau2.theta <- 1 / rgamma(1, tau2.theta.shape, scale=(1/tau2.theta.scale)) 
    
    
    #########################
    ## Sample from tau2.alpha
    #########################
    tau2.alpha.scale <- prior.tau2[2]  + sum(alpha^2)/2
    tau2.alpha <- 1 / rgamma(1, tau2.alpha.shape, scale=(1/tau2.alpha.scale)) 
    
    
    ########################
    ## Sample from tau2.delta
    ########################
    tau2.delta.scale <- prior.tau2[2]  + (2*sum(delta^2) - delta[1]^2 - delta[N]^2 - 2 * sum(delta[1:(N-1)] * delta[2:N]))/2       
    tau2.delta <- 1 / rgamma(1, tau2.delta.shape, scale=(1/tau2.delta.scale)) 
    
    
    #########################
    ## Sample from tau2.gamma
    #########################
    tau2.gamma.scale <- prior.tau2[2]  + sum(gamma.mat^2)/2
    tau2.gamma <- 1 / rgamma(1, tau2.gamma.shape, scale=(1/tau2.gamma.scale)) 
    
    
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(exp(offset.mat + regression.mat + phi.mat + theta.mat + alpha.mat + delta.mat + gamma.mat))
    deviance.all <- dpois(x=as.numeric(Y), lambda=fitted)
    deviance.all[deviance.all==0] <- min(deviance.all[deviance.all!=0])
    deviance <- -2 * sum(log(deviance.all))    
    
    
    ###################
    ## Save the results
    ###################
    if(j > burnin & (j-burnin)%%thin==0)
    {
      ele <- (j - burnin) / thin
      samples.beta[ele, ] <- beta
      samples.space[ele, ] <- phi + theta
      samples.time[ele, ] <- alpha + delta
      samples.gamma[ele, ] <- gamma
      samples.tau2[ele, ] <- c(tau2.phi, tau2.theta, tau2.alpha, tau2.delta, tau2.gamma)
      samples.deviance[ele, ] <- deviance
      samples.fitted[ele, ] <- fitted
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
      accept.phi <- 100 * accept[1] / accept[2]
      accept.theta <- 100 * accept[3] / accept[4]
      accept.alpha <- 100 * accept[5] / accept[6]
      accept.delta <- 100 * accept[7] / accept[8]
      accept.gamma <- 100 * accept[9] / accept[10]
      accept.beta <- 100 * accept[11] / accept[12]
      accept.all <- accept.all + accept
      accept <- rep(0,12)
      
      #### phi tuning parameter
      if(accept.phi > 60)
      {
        proposal.sd.phi <- 2 * proposal.sd.phi
      }else if(accept.phi < 40)              
      {
        proposal.sd.phi <- 0.5 * proposal.sd.phi
      }else
      {
      }
      
      #### theta tuning parameter
      if(accept.theta > 60)
      {
        proposal.sd.theta <- 2 * proposal.sd.theta
      }else if(accept.theta < 40)              
      {
        proposal.sd.theta <- 0.5 * proposal.sd.theta
      }else
      {
      }
      
      #### alpha tuning parameter
      if(accept.alpha > 60)
      {
        proposal.sd.alpha <- 2 * proposal.sd.alpha
      }else if(accept.alpha < 40)              
      {
        proposal.sd.alpha <- 0.5 * proposal.sd.alpha
      }else
      {
      }
      
      #### beta tuning parameter
      if(accept.beta > 60)
      {
        proposal.sd.beta <- 2 * proposal.sd.beta
      }else if(accept.beta < 40)              
      {
        proposal.sd.beta <- 0.5 * proposal.sd.beta
      }else
      {
      }
      #### gamma tuning parameter
        if(accept.gamma > 60)
        {
        proposal.sd.gamma <- 2 * proposal.sd.gamma
        }else if(accept.gamma < 40)              
        {
        proposal.sd.gamma <- 0.5 * proposal.sd.gamma
        }else
        {
        }
     #### delta tuning parameter
        if(accept.delta > 60)
        {
        proposal.sd.delta <- 2 * proposal.sd.delta
        }else if(accept.delta < 40)              
        {
        proposal.sd.delta <- 0.5 * proposal.sd.delta
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
  accept.phi <- 100 * accept.all[1] / accept.all[2]
  accept.theta <- 100 * accept.all[3] / accept.all[4]
  accept.alpha <- 100 * accept.all[5] / accept.all[6]
  accept.beta <- 100 * accept.all[7] / accept.all[8]
  accept.gamma <- 100 * accept.all[9] / accept.all[10]
  accept.delta <- 100 * accept.all[11] / accept.all[12]
  accept.final <- c(accept.beta, accept.phi, accept.theta, accept.alpha, accept.gamma, accept.delta)
  names(accept.final) <- c("beta", "phi", "theta", "alpha", "gamma", "delta")
  
  
  ## Compute DIC
  median.space <- apply(samples.space, 2, median)
  median.time <- apply(samples.time, 2, median)  
  median.space.mat <- matrix(rep(median.space, N), byrow=F, nrow=K)
  median.time.mat <- matrix(rep(median.time, K), byrow=T, nrow=K)
  median.gamma <- apply(samples.gamma, 2,median)
  median.beta <- apply(samples.beta,2,median)
  regression.mat <- matrix(X.standardised %*% median.beta, nrow=K, ncol=N, byrow=FALSE)   
  fitted.median <- as.numeric(exp(offset.mat + regression.mat + median.space.mat + median.time.mat + median.gamma))
  deviance.fitted <- -2 * sum(dpois(x=as.numeric(Y), lambda=fitted.median, log=TRUE))
  p.d <- median(samples.deviance) - deviance.fitted
  DIC <- 2 * median(samples.deviance) - deviance.fitted  
     
  
  ## Compute the LMPL
  CPO <- rep(NA, N.all)
     for(j in 1:N.all)
     {
     CPO[j] <- 1/median((1 / dpois(x=Y[j], lambda=samples.fitted[ ,j])))    
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
summary.beta <- cbind(summary.beta, rep(n.keep, p), rep(accept.beta,p))
rownames(summary.beta) <- colnames(X)
colnames(summary.beta) <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")

summary.hyper <- array(NA, c(5, 5))     
summary.hyper[1,1:3] <- quantile(samples.tau2[ ,1], c(0.5, 0.025, 0.975))
summary.hyper[2,1:3] <- quantile(samples.tau2[ ,2], c(0.5, 0.025, 0.975))
summary.hyper[3,1:3] <- quantile(samples.tau2[ ,3], c(0.5, 0.025, 0.975))
summary.hyper[4,1:3] <- quantile(samples.tau2[ ,4], c(0.5, 0.025, 0.975))
summary.hyper[5,1:3] <- quantile(samples.tau2[ ,5], c(0.5, 0.025, 0.975))
rownames(summary.hyper) <- c("tau2.phi", "tau2.theta", "tau2.alpha", "tau2.delta", "tau2.gamma")     
summary.hyper[1, 4:5] <- c(n.keep, 100)     
summary.hyper[2, 4:5] <- c(n.keep, 100)   
summary.hyper[3, 4:5] <- c(n.keep, 100)   
summary.hyper[4, 4:5] <- c(n.keep, 100)   
summary.hyper[5, 4:5] <- c(n.keep, 100)   
     
summary.results <- rbind(summary.beta, summary.hyper)
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)
     

## Compile and return the results
  modelfit <- c(DIC, p.d, LMPL)
  names(modelfit) <- c("DIC", "p.d", "LMPL")
  samples <- list(beta=mcmc(samples.beta.orig), space=mcmc(samples.space),  time=mcmc(samples.time), tau2=mcmc(samples.tau2), gamma=mcmc(samples.gamma), fitted=mcmc(samples.fitted))
model.string <- c("Likelihood model - Poisson (log link function)", "\nLatent structure model - Convolution of spatial and temporal main effects and independent interactions\n")
results <- list(formula=formula, samples=samples, fitted.values=fitted.values, residuals=residuals, stepchange=NULL, modelfit=modelfit, summary.results=summary.results, model=model.string,  accept=accept.final)
  class(results) <- "carbayesST"
     if(verbose)
     {
     b<-proc.time()
     cat(" finished in ", round(b[3]-a[3], 1), "seconds")
     }else
     {}
  return(results)
}
