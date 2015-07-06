\name{print.carbayesST}
\alias{print.carbayesST}


\title{
Print a summary of the fitted model to the screen.
}
\description{
This function takes a carbayesST object and returns a summary of the fitted model.
The summary includes posterior medians and 95 percent credible intervals, the
number of samples, the acceptance rate, the effective number of independent samples 
and the Geweke convergence diagnostic in the form of a Z-score.
}

\usage{
\method{print}{carbayesST}(x, \dots)
}



\arguments{
  \item{x}{
A carbayesST fitted model object
}
  \item{\dots}{
Ignored
}
}
\author{
Duncan Lee
}
