ST.ARCARadaptive <- function(formula, data=NULL, W, burnin=0, n.sample=1000, thin = 1, prior.mean.beta = NULL, prior.var.beta = NULL, prior.tau2 = NULL, prior.sigma2 = NULL, verbose = TRUE)
{ 
#### Check on the verbose option
     if(is.null(verbose)) verbose=TRUE     
     if(!is.logical(verbose)) stop("the verbose option is not logical.", call.=FALSE)

     if(verbose)
     {
     cat("Setting up the model\n")
     a<-proc.time()
     }else{}

     
  blocksize.beta <- 5
  blocksize.v <- 10
  z    <- which(W > 0, arr.ind = T)
  locs <- z[which(z[,1] < z[,2]), ]
  char.locs <- paste(locs[,1], ".", locs[,2], sep = "")
  n.edges <- nrow(locs)
  rho.fix = 0
  if(class(W) == "spam") W.spam <- W
  if(class(W) == "matrix") W.spam <- as.spam(W)
  if(!class(W) %in% c("matrix", "spam")) stop("W must be an object with class \"matrix\" or \"spam\"", call.=FALSE)  
  
  logit <- function(p) log(p/(1-p))
  inv_logit <- function(v) 1/(1+exp(-v))

  # interpret the formula
  frame <- try(suppressWarnings(model.frame(formula, data = data, na.action=na.pass)), silent=TRUE)
  if(class(frame)=="try-error") stop("the formula inputted contains an error, e.g the variables may be different lengths.", call.=FALSE)
  X <- try(suppressWarnings(model.matrix(object=attr(frame, "terms"), data=frame)), silent=TRUE)
  if(class(X)=="try-error") stop("the covariate matrix contains inappropriate values.", call.=FALSE)
  if(sum(is.na(X))>0) stop("the covariate matrix contains missing 'NA' values.", call.=FALSE)
  
  # set the initial parameters, get precision and beta from ls estimates
  p <- ncol(X)
  y <- model.response(frame)
  n.sites<-as.integer(nrow(W.spam))
  n.time<-as.integer(length(y)/n.sites)
  k<-as.integer(round(n.sites*n.time, 0))
     
  offset <- try(model.offset(frame), silent=TRUE)
  if(class(offset)=="try-error")   stop("the offset is not numeric.", call.=FALSE)
  if(is.null(offset))  offset <- rep(0,(n.time * n.sites))
  if(sum(is.na(offset))>0) stop("the offset has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(offset)) stop("the offset variable has non-numeric values.", call.=FALSE)
  lE <- offset
  E  <- exp(offset)
  

  
  ## Standardise the matrix
  X.standardised <- X
  X.sd <- apply(X, 2, sd)
  X.mean <- apply(X, 2, mean)
  X.indicator <- rep(NA, p)       # To determine which parameter estimates to transform back
  for(j in 1:p){
    if(length(table(X[ ,j]))>2){
      X.indicator[j] <- 1
      X.standardised[ ,j] <- (X[ ,j] - mean(X[ ,j])) / sd(X[ ,j])
    }else if(length(table(X[ ,j]))==1){
      X.indicator[j] <- 2
    }else{
      X.indicator[j] <- 0
    }
  } 
  
  # based on the blocksize.v provided, eigendecompose and create lists with relevent bits for v update
  if(is.numeric(blocksize.v)){
    ## Compute the blocking structure for v
    fromto <- seq(0, n.edges, by = blocksize.v)
    fromto[1] <- 0
    if(!n.edges %in% fromto) fromto <- c(fromto, n.edges)
    n.blocks <- length(fromto) - 1
    blockinds <- indVals <- vector("list", length = n.blocks)
    for(i in 1:n.blocks){
      diagrwcl        <- (fromto[i]+1):fromto[i+1]
      offdiagcl       <- (1:n.edges)[!(1:n.edges %in% diagrwcl)]
      indVals[[i]]    <- offdiagcl
      blockinds[[i]]  <- diagrwcl
    }
  } 
  
  # get triplet form of phi precision matrix from adjacency matrix, initial set of v's and adjacency collector
  v                  <- logit(rtrunc(n.edges, spec = "norm", mean = 0.999, sd = 0.001, a = 0, b=1))
  W_current<-W.spam
  W_current[locs][1:n.edges] <- inv_logit(v)
  W_current[locs[,2:1]][1:n.edges] <- inv_logit(v)
  Q_current          <- as.spam(diag(rowSums(W_current)) - W_current)
  tripList           <- vector("list", length = 2)
  tripList[[1]]      <- cbind(1:nrow(W_current), 1:nrow(W_current), rowSums(W_current) + 10^-7)
  tripList[[2]]      <- cbind(rbind(locs, locs[,2:1]), -rep(inv_logit(v), 2))
  tripCurrent        <- rbind(tripList[[1]], tripList[[2]])
  Q.space <- Q.space.prop <- spam(list(i = tripCurrent[,1], j = tripCurrent[,2], tripCurrent[,3]))
  # construct adjacency and precision matrices
  neighbours<-as.matrix(W_current)
  nneighbours<-rowSums(W.spam)
  Wstar <- diag(nneighbours)  - W.spam + 10^-7
  #   Q.space<-Q.space.prop<-as.spam(Wstar)
  alpha <- 0.9
  if(n.time > 1){
    Q.time <- as.spam(crossprod(diff(diag(n.time))))
    Q.time[1,1] <- Q.time[1,1] + 1
    Dg<-diag.spam(diag.spam(Q.time))
    R<-Q.time - Dg
    Dtime<-diag.spam( c(rep(1,nrow(Q.time)-1), 0))
    Dg<-Dg-Dtime
    Q.time <- Dg + Dtime*alpha^2+ R*alpha
    Q.time[n.time,n.time] <- 1
    Q<-Q.time %x% Q.space
    detTime<-determinant(Q.time, logarithm = T)
    detTime<-as.numeric(0.5*n.sites*(detTime$m)*(detTime$s))
  }  else {Q<-Q.space; Q.time <- 1; detTime = 1}
  
  # storage of parameters in the MCMC, starting values
  thin          <- thin - 1
  n.save        <- ifelse(thin == 0, (n.sample - burnin), (n.sample - burnin)/thin)
  accept.all    <- rep(0, 8)
  accept        <- accept.all
  samples.beta  <- array(NA, c(n.save, p))
  samples.phi   <- array(NA, c(n.save, n.sites * n.time))
  samples.tau2  <- samples.deviance <- samples.rho <- samples.vtau2 <- samples.alpha <- matrix(0, n.save, 1)
  samples.v     <- matrix(0, ncol = n.edges, nrow = c(n.save, n.sites*n.time))
  samples.fit   <- array(NA, c(n.save, n.sites * n.time))
  
  prior.max.tau <- 1000
  Xmat          <- as.matrix(X.standardised)
  glm_mod       <- glm(y ~-1+X.standardised, family = "poisson", offset = lE)
  beta_par      <- glm_mod$coefficients
  phi           <- log(y/E)
  XB            <- Xmat %*% beta_par
  tau           <- 0.001
  phi_tune      <- 0.5
  rho_tune      <- beta_tune <- 0.01
  W.tune        <- 1
  rho           <- ifelse(is.null(rho.fix), 0.99, rho.fix)
  tau_v         <- 100
  increment     <- 0
  
  # set spam check options to speed things up (a bit)
  spam.options( "cholsymmetrycheck" = FALSE)
  spam.options( "cholpivotcheck" = FALSE)
  spam.options( "safemode" = c(F, F, F))
  
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
  
  if(is.null(prior.sigma2)) prior.sigma2 <- c(0.001, 0.001)
  if(length(prior.sigma2)!=2) stop("the prior value for sigma2 is the wrong length.", call.=FALSE)    
  if(!is.numeric(prior.sigma2)) stop("the prior value for sigma2 is not numeric.", call.=FALSE)    
  if(sum(is.na(prior.sigma2))!=0) stop("the prior value for sigma2 has missing values.", call.=FALSE)    
     
     
  ## Compute the blocking structure for beta     
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
  proposal.sd.beta <- 0.01
  proposal.corr.beta <- solve(t(X.standardised) %*% X.standardised)
  chol.proposal.corr.beta <- chol(proposal.corr.beta)     
  
  
  # get a triplet form, cholesky and determinant for the initial ICAR precision matrix
  cholQphi      <- chol.spam(spam(list(i = tripCurrent[,1], j = tripCurrent[,2], tripCurrent[,3])))
  Qdetold       <- n.time*2*determinant(cholQphi, logarithm = T)$modulus
  # generate an initial phi vector and associated quadratic form with Q
  phiQphi       <- qformSPACETIME(Qtrip = tripCurrent, phi = phi, ntime = n.time, nsite = n.sites)
  Q             <- as.matrix(spam(list(i = tripCurrent[,1], j = tripCurrent[,2], tripCurrent[,3])))
  # get triplet form, cholesky, determinant and quadratic form with v for the edge precision matrix
  ridge         <- diag.spam(1, n.edges)
  Q_kI_trip     <- Reduce("cbind", triplet(ridge))
  v_15          <- v - 15
  vqform_current<- qform(Qtrip = Q_kI_trip, v_15)
  tau_v.shape   <- (n.edges/2) +  prior.tau2[1]
  tau_phi_shape <- (n.sites*n.time/2) + prior.sigma2[1]
  d             <- order(tripList[[2]][,1])
  
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

  
  lastblock<-(k - n.sites+1):k
  firstblock<-1:n.sites
  mult_by <- rep(seq(0, (n.time - 1 - 1)), each  = nrow(tripCurrent))
  trip_pretend_precision <- cbind(rep(tripCurrent[,1], (n.time - 1)) + mult_by*n.sites, 
                                  rep(tripCurrent[,2], (n.time - 1)) + mult_by*n.sites, 
                                  rep(tripCurrent[,3], (n.time - 1))) 
  perm          <- order(tripCurrent[,1], tripCurrent[,2])
  #   -------------------------------------------------------------------------------------------
  #   START THE MCMC SAMPLING
  #   -------------------------------------------------------------------------------------------
  
  for(j in 1:n.sample){
    # START ITERATING, ONLY SAVE thin^th ITERATION
    save.iter <- j > burnin && ((j %% thin == 0) | thin == 0)
    if(save.iter) increment <- increment+1
    
    # adjust the acceptance rate if required
    if(j %% 200 == 0){
      beta_tune   <- ifelse(accept[1] / accept[2] > 0.4,  2 * beta_tune, 0.5 * beta_tune)
      phi_tune    <- ifelse(accept[3] / accept[4] > 0.4,  2 * phi_tune, 0.5 * phi_tune)
      if(is.null(rho.fix)){ rho_tune <- ifelse(accept[5] / accept[6] > 0.4,  2 * rho_tune, 0.5 * rho_tune)}
      W.tune      <- ifelse(accept[7] / accept[8] > 0.4,  2 * W.tune, 0.5 * W.tune)
      accept.all  <- accept.all + accept
      accept      <- accept*0
    }
    
    # update ALPHA
    if(n.time > 1){
      phifirst         <- phi[-firstblock]
      philast          <- phi[-lastblock]
      trip_pretend_precision[,3] <- rep(tripCurrent[,3], (n.time - 1))                               
      philastQphilast  <- qform(trip_pretend_precision, philast)
      phifirstQphilast <- qform_asym(trip_pretend_precision, phifirst, philast) 
      mu_alpha         <- phifirstQphilast/philastQphilast
      mu_sigmasq       <- tau/philastQphilast
      alpha            <- rtrunc(n=1, spec="norm", a=10^-5, b=1 - 10^-5,  mean=mu_alpha, sd = sqrt(mu_sigmasq))
      Q.time           <- Dg + Dtime*alpha^2 + R*alpha
      Q.time[n.time,n.time] <- 1
      Q                <- Q.time %x% Q.space
      phiQphi          <- phi %*% (Q %*% phi)
      detTime          <- determinant(Q.time, logarithm = TRUE)
      detTime          <- (detTime$m)*(detTime$s)
    }
    
    # Gibbs update of tau_v
    if(vqform_current == 0) vqform_current <- 0.00001
    tau_scale <- vqform_current/2 + prior.tau2[2]
    tau_v     <- 1/rtrunc(n=1, spec="gamma", a=0.000001, b=Inf, shape=tau_v.shape, scale=(1/tau_scale))
    blocks    <- sample(1:n.blocks, n.blocks, replace = F)
    for(q in 1:n.blocks){
      #       update the blocks in a random order
      i <- blocks[q]
      vnew  <- v
      blockProposal <- rtrunc(n=length(blockinds[[i]]), spec="norm", a=-15, b=15,  mean=v[blockinds[[i]]], sd = W.tune)
      vnew[blockinds[[i]]]                 <- blockProposal
      tripStarList                         <- updatetripList(tripList, vold = v, vnew = vnew, nedges = n.edges)
      tripStar                             <- rbind(tripStarList[[1]], tripStarList[[2]])
      Q.space.prop                         <- Q.space
      Q.space.prop@entries                 <- tripStar[perm,3]
      # update the cholesky of the precision matrix & calculate the determinant
      Qstar          <- Q.time %x% Q.space.prop
      cholQphiNew    <- update(cholQphi, x = Q.space.prop) 
      detSpace       <- 2*determinant(cholQphiNew, logarithm = T)$modulus
      Qdet           <- (n.sites*detTime + n.time*detSpace)
      phiQphinew     <- phi %*% (Qstar %*% phi)
      v_15_prop      <- vnew - 15
      vqform_prop    <- qform(Qtrip = Q_kI_trip, v_15_prop)
      acceptance     <- exp(0.5*(Qdet - Qdetold) + (1/(2*tau))*(phiQphi - phiQphinew) + 0.5*(1/tau_v)*(vqform_current - vqform_prop))
      accept[8]      <- accept[8] + (1/n.blocks)
      if(runif(1)  <= acceptance){
        phiQphi          <- phiQphinew
        vqform_current   <- vqform_prop
        v                <- vnew
        tripCurrent      <- tripStar
        accept[7]        <- accept[7] + (1/n.blocks)
        tripList         <- tripStarList
        Qdetold          <- Qdet
        Q                <- Qstar
        cholQphi         <- cholQphiNew
        Q.space          <- Q.space.prop
      }
    }
    
    # update BETA
    proposal <- beta_par + (sqrt(proposal.sd.beta)* t(chol.proposal.corr.beta)) %*% rnorm(p)   
    proposal.beta <- beta_par
    offset.temp <- offset + as.numeric(phi)       
    
    for(r in 1:n.beta.block)
    {
      proposal.beta[beta.beg[r]:beta.fin[r]] <- proposal[beta.beg[r]:beta.fin[r]]
      prob <- poissonbetaupdate(X=X.standardised, nsites=n.sites, p=p, beta=beta_par, proposal=proposal.beta, 
                                offset=offset.temp, y=y, prior_meanbeta=prior.mean.beta, prior_varbeta=prior.var.beta)
      if(prob > runif(1)){
        beta_par[beta.beg[r]:beta.fin[r]] <- proposal.beta[beta.beg[r]:beta.fin[r]]
        accept[1] <- accept[1] + 1  
      }else{
        proposal.beta[beta.beg[r]:beta.fin[r]] <- beta_par[beta.beg[r]:beta.fin[r]]
      }
    }
    
    accept[2] <- accept[2] + n.beta.block    
    XB        <- X.standardised %*% beta_par
    
    
    # update PHI using one at a time M-H sampling
    tripListo     <- tripList[[2]][d,]
    tripListo[,3] <- abs(tripListo[,3])
    
    W.spam <- as.matrix(-Q.space + diag.spam(diag.spam(Q.space)))
    nneighbours <- rowSums(W.spam) 
    
    phi_update <- SPTICARphiVarb(W = as.matrix(W.spam), nsites = n.sites, ntimes = n.time, phiVarb = phi, 
                                 nneighbours = nneighbours, tau = tau, y = as.vector(y),  E = lE, 
                                 phiVarb_tune = phi_tune, 
                                 alpha = alpha, XB = as.vector(XB), beta_tune = beta_tune)
    
    phi       <- phi_update[[2]]
    phi       <- phi - mean(phi)
    accept[3] <- accept[3] + phi_update[[1]][2]
    accept[4] <- accept[4] + k
    
    # Gibbs update TAU using the gamma distribution
    Q <- Q.time %x% Q.space
    phiQphi <- phi %*% (Q %*% phi)
    tau_old<-tau
    tau_scale <- phiQphi/2 + prior.sigma2[2]
    tau     <- 1/rtrunc(n=1, spec="gamma", a=0.000001, b=Inf, shape=tau_phi_shape, scale=(1/tau_scale))
    # calculate the deviance
    fitted   <- exp(as.numeric(X.standardised %*% beta_par) + phi + lE)
    dev      <- - 2 * sum(y * log(fitted) -  fitted - lfactorial(y)) 
       
    # save samples if past burnin 
    if(save.iter){
      samples.beta[increment,]      <- beta_par
      samples.phi[increment,]       <- phi
      samples.fit[increment, ]      <- fitted
      samples.tau2[increment,]      <- tau
      samples.deviance[increment,]  <- dev
      samples.vtau2[increment,]     <- tau_v
      samples.rho[increment,]       <- rho
      samples.v[increment,]         <- v
      samples.alpha[increment,]     <- alpha
      samples.fit[increment,]       <- fitted
    }
    
    # print progress to the console
    if(j %in% percentage.points & verbose) setTxtProgressBar(progressBar, j/n.sample)
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
  accept.beta  <- 100 * accept.all[1] / accept.all[2]
  accept.phi   <- 100 * accept.all[3] / accept.all[4]
  accept.w     <- 100 * accept.all[7] / accept.all[8]
  accept.alpha <- 100
  accept.final <- c(accept.beta, accept.phi, accept.w)
  names(accept.final) <- c("beta", "phi", "w")
    
  # ## Compute information criterion (DIC, DIC3, WAIC)
  median.beta        <- apply(samples.beta, 2, median)
  regression.mat     <- matrix(X.standardised %*% median.beta, nrow = n.sites, ncol = n.time, byrow=FALSE)   
  median.phi         <- matrix(apply(samples.phi, 2, median), nrow = n.sites, ncol = n.time)
  offset.mat         <- matrix(offset, nrow = n.sites, ncol = n.time, byrow=FALSE) 
  fitted.median      <- as.numeric(exp(offset.mat + median.phi + regression.mat))
  deviance.fitted    <- -2 * sum(dpois(x=as.numeric(y), lambda=fitted.median, log=TRUE))
  p.d <- median(samples.deviance) - deviance.fitted
  DIC <- 2 * median(samples.deviance) - deviance.fitted     
  
  ## Compute the LMPL
  CPO <- rep(NA, (n.sites * n.time))
     for(j in 1:(n.sites * n.time))
     {
     CPO[j] <- 1/median((1 / dpois(x=y[j], lambda=samples.fit[ ,j])))    
     }
  LMPL <- sum(log(CPO))    
     
     
     ## Create the Fitted values
  fitted.values      <- apply(samples.fit, 2, median)
  residuals          <- as.numeric(y) - fitted.values
  
  #### transform the parameters back to the origianl covariate scale.
  samples.beta.orig <- samples.beta
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
    }
  }
  
  #### Create a summary object
  samples.beta.orig       <- mcmc(samples.beta.orig)
  summary.beta            <- t(apply(samples.beta.orig, 2, quantile, c(0.5, 0.025, 0.975))) 
  summary.beta            <- cbind(summary.beta, rep(n.save, p), rep(accept.beta,p))
  rownames(summary.beta)  <- colnames(X)
  colnames(summary.beta)  <- c("Median", "2.5%", "97.5%", "n.sample", "% accept")
  summary.hyper           <- array(NA, c(3, 5))     
  summary.hyper[1,1:3]    <- quantile(samples.tau2, c(0.5, 0.025, 0.975))
  summary.hyper[2,1:3]    <- quantile(samples.alpha, c(0.5, 0.025, 0.975))
  summary.hyper[3,1:3]    <- quantile(samples.vtau2, c(0.5, 0.025, 0.975))
  rownames(summary.hyper) <- c("sigma2", "alpha", "tau2")     
  summary.hyper[1, 4:5]   <- c(n.save, 100)     
  summary.hyper[2, 4:5]   <- c(n.save, accept.alpha)   
  summary.hyper[3, 4:5]   <- c(n.save, 100)       
  summary.results         <- rbind(summary.beta, summary.hyper)
  summary.results[ , 1:3] <- round(summary.results[ , 1:3], 4)
  summary.results[ , 4:5] <- round(summary.results[ , 4:5], 1)
  
  # convert v back to w, summarise and create a 'fitted' adjacency matrix
  samples.w <- inv_logit(samples.v)
  colnames(samples.w) <- char.locs
  get_prop_thresh <- function(v, thresh) as.numeric(!((sum(v < thresh)/length(v)) < 0.99))
  bdry99          <- apply(samples.w, 2, get_prop_thresh, thresh = 0.5)
  bdryMN          <- apply(samples.w, 2, mean)
  Wmn <- W99      <- matrix(NA, nrow = n.sites, ncol = n.sites)
  W99[locs]       <- bdry99
  Wmn[locs]       <- bdryMN
  

  ## Compile and return the results
  modelfit        <- c(DIC, p.d, LMPL)
  names(modelfit) <- c("DIC", "p.d", "LMPL")
  model.string    <- c("Likelihood model - Poisson (log link function)", 
                       "\nLatent structure model - Adaptive autoregressive CAR model\n")
  samples         <- list(beta = mcmc(samples.beta.orig), phi = mcmc(samples.phi),  
                          sigma2 = mcmc(samples.tau2), tau2 = mcmc(samples.vtau2), 
                         samples.w = samples.w, samples.alpha = mcmc(samples.alpha),
                         fitted = mcmc(samples.fit))
  stepchange <- list(Wmn = Wmn, W99 = W99)
  results         <- list(formula = formula, samples = samples, fitted.values = fitted.values, 
                          residuals = residuals, stepchange = stepchange, 
                          modelfit = modelfit, summary.results = summary.results, 
                          model = model.string,  accept = accept.final) 
  class(results) <- "carbayesST"
     if(verbose)
     {
     b<-proc.time()
     cat(" finished in ", round(b[3]-a[3], 1), "seconds")
     }else
     {}
  return(results)
}







