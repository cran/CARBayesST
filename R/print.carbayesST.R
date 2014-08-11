print.carbayesST <- function(x,...)
{
     if(is.null(x$median.Z))
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
     cat("\nNumber of clusters with the number of data points in each one")
     print(table(x$median.Z))
     }
     
return(invisible(x))
}