print.carbayesST <- function(x,...)
{
     if(class(x$stepchange)=="list")
     {
     #### Print out the model fitted
     cat("\n#################\n")
     cat("#### Model fitted\n")
     cat("#################\n")
     cat(x$model)
     cat("Regression equation - ")
     print(x$formula)

     #### Print out the results
     cat("\n############\n")
     cat("#### Results\n")
     cat("############\n")
     cat("Posterior quantiles for selected parameters and DIC\n\n")
     print(x$summary.results)
     cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "\n")
     cat("\nThe number of stepchanges identified in the random effect surface")
     cat("\nthat satisfy Prob(w_ij < 0.5|data) > 0.99 is \n")
     temp <- x$stepchange[[2]][!is.na(x$stepchange[[2]])]
     tab <- array(NA, c(1,2))
     tab[1, ] <- c(sum(temp), length(temp)- sum(temp))
     colnames(tab) <- c("stepchange", "no stepchange")
     print(tab)
     #print(sum(temp))
     #cat(" which is ")
     #print(round(100 * mean(temp),2))
     #cat("%\n")
     }else if(class(x$stepchange)=="matrix")
     {
     #### Print out the model fitted
     cat("\n#################\n")
     cat("#### Model fitted\n")
     cat("#################\n")
     cat(x$model)
     cat("Regression equation - ")
     print(x$formula)

     #### Print out the results
     cat("\n############\n")
     cat("#### Results\n")
     cat("############\n")
     cat("Posterior quantiles for selected parameters and DIC\n\n")
     print(x$summary.results)
     cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "\n")
     cat("\nNumber of clusters with the number of data points in each one")
     print(table(x$stepchange))
     }else
     {
     #### Print out the model fitted
     cat("\n#################\n")
     cat("#### Model fitted\n")
     cat("#################\n")
     cat(x$model)
     cat("Regression equation - ")
     print(x$formula)

     #### Print out the results
     cat("\n############\n")
     cat("#### Results\n")
     cat("############\n")
     cat("Posterior quantiles for selected parameters and DIC\n\n")
     print(x$summary.results)
     cat("\nDIC = ", x$modelfit[1], "     ", "p.d = ", x$modelfit[2], "\n")
     }
     
return(invisible(x))
}