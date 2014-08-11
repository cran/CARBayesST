results.summarise <- function(samples, exceedences=NULL, quantiles=0.5)
{
#### Compute exceedence probabilities and posterior quantiles
#### for the matrix of samples
     if(class(samples) == "mcmc")     
     {
     N.all <- ncol(samples)
     n.sample <- nrow(samples)
          if(missing(exceedences))
          {
          N.exceedence <- 0     
          }else
          {
          N.exceedence <- length(exceedences)
          }
          
          if(missing(quantiles))
          {
          N.quantiles <- 0     
          }else
          {
          N.quantiles <- length(quantiles)
          }

        if(N.exceedence>0 & N.quantiles>0) 
        {
         posterior.exceedence <- array(NA, c(N.all,N.exceedence))
         posterior.quantile <- array(NA, c(N.all,N.quantiles))
         colnames(posterior.exceedence) <- exceedences
         colnames(posterior.quantile) <- quantiles
             
            for(j in 1:N.all)
            {
             mcmc <- samples[ ,j]
             posterior.quantile[j, ] <- quantile(mcmc, quantiles)
                 for(k in 1:N.exceedence)
                 {
                 posterior.exceedence[j, k] <- length(which(mcmc>exceedences[k]))/n.sample       
                 }
            }
          }else if(N.exceedence==0 & N.quantiles>0)
          {
          posterior.exceedence <- NULL
          posterior.quantile <- array(NA, c(N.all,N.quantiles))
          colnames(posterior.quantile) <- quantiles
             
            for(j in 1:N.all)
            {
             mcmc <- samples[ ,j]
             posterior.quantile[j, ] <- quantile(mcmc, quantiles)
            }    
          }else if(N.exceedence>0 & N.quantiles==0)
          {
          posterior.exceedence <- array(NA, c(N.all,N.exceedence))
          posterior.quantile <- NULL
          colnames(posterior.exceedence) <- exceedences
             
            for(j in 1:N.all)
            {
             mcmc <- samples[ ,j]
                  for(k in 1:N.exceedence)
                 {
                 posterior.exceedence[j, k] <- length(which(mcmc>exceedences[k]))/n.sample       
                 }
            }
          }else
          {
          posterior.quantile <- NULL
          posterior.exceedence <- NULL
          }
       
     results <- list(quantile=posterior.quantile, exceedence=posterior.exceedence)
     }else
     {
      stop("The samples object is not a mcmc type.", call.=FALSE)    
     }     

return(results)
}
 
