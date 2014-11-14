ST.ARCARcluster <- function(formula, data=NULL, W, G, burnin=0, n.sample=1000, thin=1,  prior.mean.beta=NULL, prior.var.beta=NULL, prior.tau2=NULL, prior.sigma2=NULL, prior.alpha=NULL, prior.delta=NULL, verbose=TRUE)
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
## Create the offset
offset <- try(model.offset(frame), silent=TRUE)
    if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
    if(is.null(offset))  offset <- rep(0,N.all)
    if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
    if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)

     
  #### Format and check the number of clusters G     
  if(length(G)!=1) stop("G is the wrong length.", call.=FALSE)    
  if(!is.numeric(G)) stop("G is not numeric.", call.=FALSE)    
  if(G<=0) stop("G is not positive.", call.=FALSE)    
  if(G!=round(G)) stop("G is not an integer.", call.=FALSE) 
Gstar <- (G+1)/2  
  
     
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
 
   if(is.null(prior.sigma2)) prior.sigma2 <- c(0.001, 0.001)
  if(length(prior.sigma2)!=2) stop("the prior value for sigma2 is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.sigma2)) stop("the prior value for sigma2 is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.sigma2))!=0) stop("the prior value for sigma2 has missing values.", call.=FALSE)    
  
  if(is.null(prior.tau2)) prior.tau2 <- c(0.001, 0.001)
  if(length(prior.tau2)!=2) stop("the prior value for tau2 is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.tau2)) stop("the prior value for tau2 is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.tau2))!=0) stop("the prior value for tau2 has missing values.", call.=FALSE)    
  
  if(is.null(prior.alpha)) prior.alpha <- 10       
  if(length(prior.alpha)!=1) stop("the prior value for alpha is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.alpha)) stop("the prior value for alpha is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.alpha))!=0) stop("the prior value for alpha has missing values.", call.=FALSE)    
  if(prior.alpha<=0) stop("the prior value for alpha is not positive.", call.=FALSE)    
     
  if(is.null(prior.delta)) prior.delta <- 10       
  if(length(prior.delta)!=1) stop("the prior value for delta is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.delta)) stop("the prior value for delta is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.delta))!=0) stop("the prior value for delta has missing values.", call.=FALSE)    
  if(prior.delta<=0) stop("the prior value for delta is not positive.", call.=FALSE)    
 
     
  #### Specify the initial parameter values
  beta <- glm(Y~X.standardised-1, offset=offset, family=poisson)$coefficients
  beta[which(X.indicator==2)] <- 0     
  res.temp <- log.Y - X.standardised %*% beta - offset
  res.temp.mat <- matrix(res.temp, ncol=N, nrow=K, byrow=FALSE)
  lambda <- array(NA, c(N,G))
  Z <- array(NA, c(K,N)) 
     for(t in 1:N)
     {
    kmean <- kmeans(x=as.numeric(res.temp.mat[ ,t]), centers=G, nstart=1000)     
    lambda.temp <- kmean$centers 
    lambda.order <- order(lambda.temp)
    lambda[t, ] <- sort(lambda.temp)
    cluster.temp <- kmean$cluster     
    Z.vec <- cluster.temp    
          for(i in 1:K)
          {
          Z[i,t] <- which(lambda.order==Z.vec[i])
          }     
     }
   
  alpha <- runif(1,0, prior.alpha)
  delta <- runif(1,0, prior.delta)
  mu <- array(NA, c(K,N))
  phi <- array(NA, c(K,N))
     for(t in 1:N)
     {
     lambdatemp <- lambda[t, ]
     mu[ ,t] <-  lambdatemp[Z[ ,t]]
     phi[ ,t] <- res.temp.mat[ ,t] - mu[ ,t]
     } 
  sigma2 <- var(as.numeric(diff(lambda)))  
  tau2 <- var(as.numeric(phi)) / mean(apply(W,2,sum))        
  rho <- runif(1) 
  gamma <- runif(1)
  offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
  regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
  Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE) 

     
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
  samples.phi <- array(NA, c(n.keep, N.all))
  samples.beta <- array(NA, c(n.keep, p))
  samples.Z <- array(NA, c(n.keep, N.all))
  samples.lambda <- array(NA, c(n.keep, (N*G)))
  samples.tau2 <- array(NA, c(n.keep, 1))
  samples.sigma2 <- array(NA, c(n.keep, 1))   
  samples.alpha <- array(NA, c(n.keep, 1))
  samples.delta <- array(NA, c(n.keep, 1))
  samples.gamma <- array(NA, c(n.keep, 1))
  samples.rho <- array(NA, c(n.keep, 1))
  samples.fitted <- array(NA, c(n.keep, N.all))
  samples.deviance <- array(NA, c(n.keep, 1))
  
     
  #### Specify the Metropolis quantities
  accept.all <- rep(0,12)
  accept <- accept.all
  proposal.sd.phi <- 0.1
  proposal.sd.alpha <- 0.1
  proposal.sd.delta <- 0.1
  proposal.sd.lambda <- 0.1
  proposal.sd.rho <- 0.05
  proposal.sd.beta <- 0.01
  proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
  chol.proposal.corr.beta <- chol(proposal.corr.beta) 
  tau2.shape <- prior.tau2[1] + N.all/2
  sigma2.shape <- prior.sigma2[1] + 0.5 * G * (N-1)      
  
  
  #### Create the sparse forms for the W matrix
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

  spacelist <- as.list(rep(NA,K))     
    for(i in 1:K)
    {
    spacelist[[i]] <- which(W[i, ]==1)     
    }
  
  
  ## Create the determinant for Q     
  Wstar <- diag(n.neighbours.space) - W
  Wstar.eigen <- eigen(Wstar)
  Wstar.val <- Wstar.eigen$values
  det.Q <-  0.5 * sum(log((rho * Wstar.val + (1-rho))))    
  
  
  
  ### Create the first order random walk matrix V
  V <- array(0, c(N,N))
  Vneigh <- array(0, c(N,N))
  diag(V) <- c(1, rep(2, N-2), 1)
     for(i in 2:N)
     {
    V[i, (i-1)] <- -1
    V[(i-1),i] <- -1
    Vneigh[i, (i-1)] <- 1
    Vneigh[(i-1),i] <- 1     
     }
  

  #### Specify quantities that do not change
  offset.mat <- matrix(offset, nrow=K, ncol=N, byrow=FALSE) 
  regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)
  Y.mat <- matrix(Y, nrow=K, ncol=N, byrow=FALSE) 
     
     
     
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
    ## Sample from phi
    ##################
    phi.offset <- exp(mu + offset.mat + regression.mat)
    den.offset <- rho * n.neighbours.space + 1 - rho
    temp1 <- poissonarcarupdate(W_list=spacelist, nsites=K, ntime=N, phi=phi, tau2=tau2, gamma=gamma, rho=rho, ymat=Y.mat, phi_tune=proposal.sd.phi, offset=phi.offset, denoffset=den.offset)      
    phi <- temp1[[1]]
    accept[1] <- accept[1] + temp1[[2]]
    accept[2] <- accept[2] + N.all    

     for(r in 1:G)
     {
     temp1 <- which(as.numeric(Z)==r)    
          if(length(temp1)>0)
          {
          mean.temp <- mean(as.numeric(phi)[temp1])
          phi[Z==r] <- phi[Z==r] - mean.temp
          }else{}
     }
          

       
    ###########################
    ## Sample from tau2 and rho
    ###########################
    proposal.rho <- rtrunc(n=1, spec="norm", a=0, b=1, mean=rho, sd=proposal.sd.rho)    
    QF <- quadform(W_duplet1=space.duplet[ ,1], W_duplet2=space.duplet[ ,2], n_duplet=n.space.duplet,  nsites=K, phi=phi[ ,1], theta=phi[ ,1], nneighbours=n.neighbours.space, diagonal=rho, offdiagonal=rho)
    QF.prop <- quadform(W_duplet1=space.duplet[ ,1], W_duplet2=space.duplet[ ,2], n_duplet=n.space.duplet,  nsites=K, phi=phi[ ,1], theta=phi[ ,1], nneighbours=n.neighbours.space, diagonal=proposal.rho, offdiagonal=proposal.rho)
                
       for(t in 2:N)
       {
       ## Compute the quadratic forms
       phi.scaled <- phi[ ,t] - gamma * phi[ ,(t-1)]
       temp2 <- quadform(W_duplet1=space.duplet[ ,1], W_duplet2=space.duplet[ ,2], n_duplet=n.space.duplet,  nsites=K, phi=phi.scaled, theta=phi.scaled, nneighbours=n.neighbours.space, diagonal=rho, offdiagonal=rho)      
       temp3 <- quadform(W_duplet1=space.duplet[ ,1], W_duplet2=space.duplet[ ,2], n_duplet=n.space.duplet,  nsites=K, phi=phi.scaled, theta=phi.scaled, nneighbours=n.neighbours.space, diagonal=proposal.rho, offdiagonal=proposal.rho)      
       QF <- QF + temp2
       QF.prop <- QF.prop + temp3
       }
    
         
     #### Sample from tau2
     tau2.scale <- QF + prior.tau2[2] 
     tau2 <- 1 / rgamma(1, tau2.shape, scale=(1/tau2.scale)) 

         
     ## Update rho  
     det.Q.proposal <- 0.5 * sum(log((proposal.rho * Wstar.val + (1-proposal.rho))))
     logprob.current <- N * det.Q - QF / tau2
     logprob.proposal <- N * det.Q.proposal - QF.prop / tau2
     prob <- exp(logprob.proposal - logprob.current)
          if(prob > runif(1))
          {
          rho <- proposal.rho
          det.Q <- det.Q.proposal
          accept[7] <- accept[7] + 1           
          }else
          {
          }              
    accept[8] <- accept[8] + 1           



    ####################
    ## Sample from gamma
    ####################
    num <- 0
    den <- 0
       for(t in 2:N)
       {
       ## Compute the quadratic forms
       num <- num + quadform(W_duplet1=space.duplet[ ,1], W_duplet2=space.duplet[ ,2], n_duplet=n.space.duplet,  nsites=K, phi=phi[ ,t], theta=phi[ ,(t-1)], nneighbours=n.neighbours.space, diagonal=rho, offdiagonal=rho)
       den <- den + quadform(W_duplet1=space.duplet[ ,1], W_duplet2=space.duplet[ ,2], n_duplet=n.space.duplet,  nsites=K, phi=phi[ ,(t-1)], theta=phi[ ,(t-1)], nneighbours=n.neighbours.space, diagonal=rho, offdiagonal=rho)
       }
    mean.gamma <- num/den
    sd.gamma <- sqrt(tau2/den)
    gamma <- rnorm(n=1, mean=mean.gamma, sd=sd.gamma)         

           
    
    #######################     
    #### Sample from lambda
    #######################
     for(t in 1:N)
     {
      #### Propose a new value
      lambdatemp <- lambda[t, ]
      proposal.extend <- c(-100, lambdatemp, 100) 
      for(r in 1:G)
      {
        proposal.extend[(r+1)] <- rtrunc(n=1, spec="norm", a=proposal.extend[r], b=proposal.extend[(r+2)], mean=proposal.extend[(r+1)], sd=proposal.sd.lambda)
      }
      proposal <- proposal.extend[-c(1, (G+2))]
      
      #### Compute the data likelihood
      lp.current <- lambdatemp[Z[ ,t]] + as.numeric(offset.mat[ ,t]) + as.numeric(phi[ ,t]) + as.numeric(regression.mat[ ,t])   
      lp.proposal <- proposal[Z[ ,t]] + as.numeric(offset.mat[ ,t]) + as.numeric(phi[ ,t]) + as.numeric(regression.mat[ ,t])   
      prob1 <- sum(as.numeric(Y.mat[ ,t]) * lp.current - exp(lp.current))
      prob2 <- sum(as.numeric(Y.mat[ ,t]) * lp.proposal - exp(lp.proposal))
      
      #### Compute the prior     
      var.prior <- rep(sigma2 / V[t,t],G)               
      mean.prior <- Vneigh[t, ] %*% lambda / V[t,t]               
      prob3 <- sum(((lambdatemp - mean.prior)^2 - (proposal - mean.prior)^2) / (2 * var.prior))                            
      
      #### Compute the acceptance probability     
      prob <- exp(prob2 - prob1 + prob3)
      if(prob > runif(1))
      {
        lambda[t, ] <- proposal
        mu[ ,t] <- proposal[Z[ ,t]] 
        accept[5] <- accept[5] + 1  
      }else
      {
        proposal <-  lambda
      }
      accept[6] <- accept[6] + 1           
    }
    
    
    ##################     
    #### Sample from Z
    ##################
    mu.offset <- exp(offset.mat + regression.mat + phi)
    test <- Zupdate(Z=Z, Offset=mu.offset, Y=Y.mat, alpha=alpha, lambda=lambda, nsites=K, ntime=N, G=G, SS=1:G, Gstar, delta=delta)          
    Z <- test
    for(t in 1:N)
    {
      lambdatemp <- lambda[t, ]
      mu[ ,t] <-  lambdatemp[Z[ ,t]]    
    } 
    
    
    ######################
    #### Sample from alpha
    ######################
    proposal.alpha <-  rtrunc(n=1, spec="norm", a=0, b=prior.alpha, mean=alpha, sd=proposal.sd.alpha)    
    logprob1 <- sum(t(apply(Z, 1, diff))^2) * (alpha - proposal.alpha)
    Z.vec <- as.numeric(Z[ ,-N])       
    norm1 <- norm(Z=Z.vec, G=G, Gstar=Gstar, alpha=alpha, delta=delta, SS=1:G, K=K, Nall=length(Z.vec))
    norm2 <- norm(Z=Z.vec, G=G, Gstar=Gstar, alpha=proposal.alpha, delta=delta, SS=1:G, K=K, Nall=length(Z.vec))       
    prob <- exp(logprob1 + norm1[2] - norm2[2])       
    if(prob > runif(1))
    {
      alpha <- proposal.alpha
      accept[3] <- accept[3] + 1  
    }else
    {
    }
    accept[4] <- accept[4] + 1  
       
    
    ######################
    #### Sample from delta
    ######################
    proposal.delta <-  rtrunc(n=1, spec="norm", a=0, b=prior.delta, mean=delta, sd=proposal.sd.delta)    
    logprob1 <- sum((Z-Gstar)^2) * (delta - proposal.delta)     
    norm1 <- norm(Z=Z.vec, G=G, Gstar=Gstar, alpha=alpha, delta=delta, SS=1:G, K=K, Nall=length(Z.vec))
    norm2 <- norm(Z=Z.vec, G=G, Gstar=Gstar, alpha=alpha, delta=proposal.delta, SS=1:G, K=K, Nall=length(Z.vec))       
    prob <- exp(logprob1 + sum(norm1) - sum(norm2))       
    if(prob > runif(1))
    {
      delta <- proposal.delta
      accept[11] <- accept[11] + 1  
    }else
    {
    }
    accept[12] <- accept[12] + 1           
  
       
       
    #######################
    #### Sample from sigma2
    #######################
    SS <- 0
    for(r in 1:G)
    {
      SS <- SS + t(lambda[ ,r]) %*% V %*% lambda[ ,r]     
    }
    sigma2.scale <- 0.5 * SS + prior.sigma2[2] 
    sigma2 <- 1 / rgamma(1, sigma2.shape, scale=(1/sigma2.scale))
    
    
    
    ####################
    ## Sample from beta
    ####################
    proposal <- beta + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)
    proposal[which(X.indicator==2)] <- 0     
    proposal.beta <- beta
    offset.temp <- as.numeric(phi) + offset + as.numeric(mu)
       
       
       for(r in 1:n.beta.block)
       {
       proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
       prob <- poissonbetaupdate(X.standardised, N.all, p, beta, proposal.beta, offset.temp, Y, prior.mean.beta, prior.var.beta)
            if(prob > runif(1))
            {
            beta[beta.beg[r]:beta.fin[r]] <- proposal.beta[beta.beg[r]:beta.fin[r]]
            accept[9] <- accept[9] + 1  
            }else
            {
            proposal.beta[beta.beg[r]:beta.fin[r]] <- beta[beta.beg[r]:beta.fin[r]]
            }
        }

    accept[10] <- accept[10] + n.beta.block    
    regression.mat <- matrix(X.standardised %*% beta, nrow=K, ncol=N, byrow=FALSE)   
  
       
    #########################
    ## Calculate the deviance
    #########################
    fitted <- as.numeric(exp(mu + offset.mat + phi + regression.mat))             
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
      samples.phi[ele, ] <- as.numeric(phi)
      samples.tau2[ele, ] <- tau2
      samples.sigma2[ele, ] <- sigma2
      samples.rho[ele, ] <- rho
      samples.alpha[ele, ] <- alpha
      samples.gamma[ele, ] <- gamma
      samples.delta[ele, ] <- delta
      samples.lambda[ele, ] <- as.numeric(lambda)
      samples.Z[ele, ] <- as.numeric(Z)
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
      accept.lambda <- 100 * accept[5] / accept[6]
      accept.alpha <- 100 * accept[3] / accept[4]
      accept.rho <- 100 * accept[7] / accept[8]
      accept.beta <- 100 * accept[9] / accept[10]
      accept.delta <- 100 * accept[11] / accept[12]
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
      
      #### lambda tuning parameter
      if(accept.lambda > 70)
      {
        proposal.sd.lambda <- 2 * proposal.sd.lambda
      }else if(accept.lambda < 40)              
      {
        proposal.sd.lambda <- 0.5 * proposal.sd.lambda
      }else
      {
      }
      
      #### alpha tuning parameter               
      if(accept.alpha > 70)
      {
        proposal.sd.alpha <- 2 * proposal.sd.alpha
      }else if(accept.alpha < 40)              
      {
        proposal.sd.alpha <- 0.5 * proposal.sd.alpha
      }else
      {
      }
      
      #### rho tuning parameter
      if(accept.rho > 70)
      {
        proposal.sd.rho <- min(proposal.sd.rho, 0.5)
      }else if(accept.rho < 50)              
      {
        proposal.sd.rho <- 0.5 * proposal.sd.rho
      }else
      {
      }
         
     #### beta tuning parameter
     if(accept.beta > 70)
     {
     proposal.sd.beta <- 2 * proposal.sd.beta
     }else if(accept.beta < 50)              
     {
     proposal.sd.beta <- 0.5 * proposal.sd.beta
     }else
     {
     }
         
      #### delta tuning parameter               
      if(accept.delta > 70)
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
  accept.lambda <- 100 * accept.all[5] / accept.all[6]
  accept.alpha <- 100 * accept.all[3] / accept.all[4]
  accept.rho <- 100 * accept.all[7] / accept.all[8]
  accept.beta <- 100 * accept.all[9] / accept.all[10]
  accept.delta <- 100 * accept.all[11] / accept.all[12]
  accept.final <- c(accept.beta, accept.lambda, accept.alpha, accept.phi, accept.rho, accept.delta)
  names(accept.final) <- c("beta", "lambda", "alpha", "phi", "rho", "delta")
  
  
  ## Summarise the Z results
  median.Z <- matrix(apply(samples.Z,2,median), nrow=K, ncol=N)     
  
  
  ## Compute information criterion (DIC, DIC3, WAIC)
  median.lambda <- matrix(apply(samples.lambda, 2, median), nrow=N, ncol=G)
  median.mu <- array(NA, c(K,N))
  for(t in 1:N)
  {
    lambdatemp <- median.lambda[t, ]
    median.mu[ ,t] <-  lambdatemp[median.Z[ ,t]]    
  }
  
  median.beta <- apply(samples.beta,2,median)
  regression.mat <- matrix(X.standardised %*% median.beta, nrow=K, ncol=N, byrow=FALSE)   
  median.phi <- matrix(apply(samples.phi, 2, median), nrow=K, ncol=N)
  fitted.median <- as.numeric(exp(median.mu + offset.mat + median.phi + regression.mat))
  samples.deviance[samples.deviance==0] <- min(samples.deviance[samples.deviance!=0])
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

summary.hyper <- array(NA, c(6, 5))     
summary.hyper[1,1:3] <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
summary.hyper[2,1:3] <- quantile(samples.rho, c(0.5, 0.025, 0.975))
summary.hyper[3,1:3] <- quantile(samples.sigma2, c(0.5, 0.025, 0.975))
summary.hyper[4,1:3] <- quantile(samples.alpha, c(0.5, 0.025, 0.975))
summary.hyper[5,1:3] <- quantile(samples.delta, c(0.5, 0.025, 0.975))
summary.hyper[6,1:3] <- quantile(samples.gamma, c(0.5, 0.025, 0.975))
rownames(summary.hyper) <- c("tau2", "rho", "sigma2", "alpha", "delta", "gamma")     
summary.hyper[1, 4:5] <- c(n.keep, 100)     
summary.hyper[2, 4:5] <- c(n.keep, accept.rho)   
summary.hyper[3, 4:5] <- c(n.keep, 100)   
summary.hyper[4, 4:5] <- c(n.keep, accept.alpha)   
summary.hyper[5, 4:5] <- c(n.keep, accept.delta)   
summary.hyper[6, 4:5] <- c(n.keep, 100)   
     
summary.results <- rbind(summary.beta, summary.hyper)
summary.results <- summary.results[-which(X.indicator==2), ]
summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)

     
## Compile and return the results
  modelfit <- c(DIC, p.d, LMPL)
  names(modelfit) <- c("DIC", "p.d", "LMPL")

     if(length(which(X.indicator==2))==p)
     {
     samples <- list(lambda=mcmc(samples.lambda),  alpha=mcmc(samples.alpha), Z=mcmc(samples.Z), sigma2=mcmc(samples.sigma2), phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), rho=mcmc(samples.rho), gamma=mcmc(samples.gamma), delta=mcmc(samples.delta), fitted=mcmc(samples.fitted))
     }else
     {
     samples <- list(beta=mcmc(samples.beta.orig[ ,-which(X.indicator==2)]), lambda=mcmc(samples.lambda),  alpha=mcmc(samples.alpha), Z=mcmc(samples.Z), sigma2=mcmc(samples.sigma2), phi=mcmc(samples.phi), tau2=mcmc(samples.tau2), rho=mcmc(samples.rho), gamma=mcmc(samples.gamma), delta=mcmc(samples.delta), fitted=mcmc(samples.fitted))
     }
model.string <- c("Likelihood model - Poisson (log link function)", "\nLatent structure model - Clustering model with autoregressive CAR random effects\n")
     results <- list(formula=formula, samples=samples, fitted.values=fitted.values, residuals=residuals, stepchange=median.Z, modelfit=modelfit, summary.results=summary.results, model=model.string,  accept=accept.final)
  class(results) <- "carbayesST"
     if(verbose)
     {
     b<-proc.time()
     cat(" finished in ", round(b[3]-a[3], 1), "seconds")
     }else
     {}
  return(results)
}
