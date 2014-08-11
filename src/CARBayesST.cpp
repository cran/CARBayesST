#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector linpredcompute(NumericMatrix X, const int nsites, const int p, 
                          NumericVector beta, NumericVector offset)
{
//Create new objects
// Compute the linear predictor
NumericVector linpred(nsites);
double temp; 


//  Compute the linear predictor via a double for loop
     for(int j = 0; j < nsites; j++)
     {
     temp = 0;
      
          for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];     
          
     linpred[j] = temp + offset[j];  
     }


// Return the result
return linpred;
}

// [[Rcpp::export]]
double quadform(IntegerVector W_duplet1, IntegerVector W_duplet2, const int n_duplet, const int nsites, 
                NumericVector phi, NumericVector nneighbours, double diagonal, double offdiagonal)
{
  // Compute a quadratic form for the random effects
  // Create new objects 
  double tau2_posteriorscale;
  double tau2_quadform = 0, tau2_phisq = 0;
  int row, col;
  
  
  // Compute the off diagonal elements of the quadratic form
  for(int l = 0; l < n_duplet; l++)
  {
    row = W_duplet1[l] - 1;
    col = W_duplet2[l] - 1;
    tau2_quadform = tau2_quadform + phi[row] * phi[col]; 
  }
  
  
  // Compute the diagonal elements of the quadratic form          
  for(int l = 0; l < nsites; l++)
  {
    tau2_phisq =  tau2_phisq + pow(phi[l],2) * (diagonal * nneighbours[l] + 1 - diagonal);
  }
  
  
  // Compute the quadratic form
  tau2_posteriorscale = 0.5 * (tau2_phisq - offdiagonal * tau2_quadform);
  
  // Return the simulated value
  return tau2_posteriorscale;
}





// [[Rcpp::export]]
List poissoncarupdate(List W_list, const int nsites, NumericVector phi,double tau2, 
const NumericVector y, const double phi_tune, double rho_num, double rho_den, 
NumericVector offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priormean, priorvar, priorvardenom, propphi;
NumericVector phinew(nsites);
 
   
//  Update each random effect in turn
phinew = phi;
     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     IntegerVector neighbourvec = W_list[j];
     int m = neighbourvec.size();
     sumphi = 0;
          for(int l = 0; l < m; l++) sumphi += phinew[(neighbourvec[l]-1)];
      priorvardenom = (m * rho_den + (1-rho_den));
      priorvar = tau2 / priorvardenom;
      priormean = rho_num * sumphi / priorvardenom; 
      
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      oldlikebit = phinew[j] * y[j]  - exp(phinew[j]) * offset[j];
      newlikebit =  propphi * y[j]  - exp(propphi) * offset[j];
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew[j] = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
    }


List out(2);
out[0] = phinew;
out[1] = accept;
return out;
}




// [[Rcpp::export]]
List poissonindepupdate(const int nsites, NumericVector theta,double tau2, 
const NumericVector y, const double theta_tune, NumericVector offset)
{
// Update the spatially independent random effects 
//Create new objects
int accept=0;
double acceptance, proptheta;
double priorbit, oldlikebit, newlikebit;
NumericVector thetanew(nsites);
 
   
//  Update each random effect in turn
thetanew = theta;
     for(int j = 0; j < nsites; j++)
     {
      // propose a value  
      proptheta = rnorm(1, thetanew[j], sqrt(theta_tune))[0];
      
      // Accept or reject it
      priorbit = (0.5/tau2) * (pow(thetanew[j], 2) - pow(proptheta, 2));
      oldlikebit = thetanew[j] * y[j]  - exp(thetanew[j]) * offset[j];
      newlikebit =  proptheta * y[j]  - exp(proptheta) * offset[j];
      acceptance = exp(priorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          thetanew[j] = proptheta;
          accept = accept + 1;
          }
          else
          { 
          }
    }


List out(2);
out[0] = thetanew;
out[1] = accept;
return out;
}



// [[Rcpp::export]]
double poissonbetaupdate(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
               NumericVector proposal, NumericVector offset, NumericVector y, 
               NumericVector prior_meanbeta, NumericVector prior_varbeta)
{
// Compute the acceptance probability for beta
//Create new objects
double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit=0;
NumericVector lp_current(nsites), lp_proposal(nsites), fitted_current(nsites), fitted_proposal(nsites);


// Create the log likelihood acceptance probability component
lp_current = linpredcompute(X, nsites, p, beta, offset);
lp_proposal = linpredcompute(X, nsites, p, proposal, offset);     
     for(int j = 0; j < nsites; j++)     
     {
     fitted_current[j] = exp(lp_current[j]);
     fitted_proposal[j] = exp(lp_proposal[j]);
     oldlikebit = oldlikebit + y[j] * log(fitted_current[j]) - fitted_current[j];
     newlikebit = newlikebit + y[j] * log(fitted_proposal[j]) - fitted_proposal[j];
     }
likebit = newlikebit - oldlikebit;


// Create the prior acceptance component
     for(int j = 0; j < p; j++)     
     {
     priorbit = priorbit + 0.5 * pow((beta[j]-prior_meanbeta[j]),2) / prior_varbeta[j] - 0.5 * pow((proposal[j]-prior_meanbeta[j]),2) / prior_varbeta[j];
     }


// Compute the acceptance probability and return the value
acceptance = exp(likebit + priorbit);
return acceptance;
}




// [[Rcpp::export]]
NumericMatrix Zupdate(NumericMatrix Z, NumericMatrix Offset, NumericMatrix Y, const double alpha, 
NumericMatrix lambda, const int nsites, const int ntime, const int G, NumericVector SS,double Gstar)
{
//Create the new copy of Z
NumericMatrix Znew(nsites, ntime);

// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector lambdacurrent(G), prior3(G), prior4(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
     
// Update the elements in Z during time period 1
lambdacurrent = lambda(0 , _);

     for(int k = 0; k < nsites; k++) 
     {
     // Compute the Poisson likelihood component
     like1 = Y(k,0) * lambdacurrent - exp(lambdacurrent) * Offset(k,0);

     // Compute the prior part
     //prior1 = exp(-alpha * abs(SS - Znew(k,1)) - abs(SS - Gstar));
     prior1 = exp(-alpha * pow(SS - Znew(k,1),2) - pow(SS - Gstar,2));
     prior2 = log(prior1 / sum(prior1));    
     //prior3 = exp(- abs(SS - Gstar));
     prior3 = exp(- pow(SS - Gstar,2));
     prior4 = log(prior3 / sum(prior3));    

     // Compute the posterior probability
     posterior1 = like1 + prior2 + prior4;
     posterior2 = posterior1 - max(posterior1);
     posterior3 = exp(posterior2) / sum(exp(posterior2));
     u = runif(1)[0];          
     cs = posterior3[0];
     test = 1;
     counter = 1;
          while(test==1)
          {
               if(cs>u)
               {
               test = 0;     
               }else
               {
               counter = counter + 1;
               cs = cs + posterior3[(counter-1)];     
               }
          }
     Znew(k,0) = counter;
     }



//  Update the elements in Z during the remaining time periods except time ntime
     for(int j = 1; j < ntimeminus; j++)
     {
     lambdacurrent = lambda(j , _);     
          
          for(int k = 0; k < nsites; k++) 
          {
          // Compute the Poisson likelihood component
          like1 = Y(k,j) * lambdacurrent - exp(lambdacurrent) * Offset(k,j);
          
          // Compute the Markov prior component
          //prior1 = exp(-alpha * abs(SS - Znew(k,(j-1))) - abs(SS - Gstar));
          prior1 = exp(-alpha * pow(SS - Znew(k,(j-1)),2) - pow(SS - Gstar,2));
          prior2 = log(prior1 / sum(prior1));    
          //prior3 = exp(-alpha * abs(SS - Znew(k,(j+1))) - abs(SS - Gstar));
          prior3 = exp(-alpha * pow(SS - Znew(k,(j+1)),2) - pow(SS - Gstar,2));
          prior4 = log(prior3 / sum(prior3));    
          
          // Compute the posterior probabilities
          posterior1 = like1 + prior2 + prior4;
          posterior2 = posterior1 - max(posterior1);
          posterior3 = exp(posterior2) / sum(exp(posterior2));
          u=runif(1)[0];          
          cs = posterior3[0];
          test = 1;
          counter = 1;
               while(test==1)
               {
                    if(cs>u)
                    {
                    test = 0;     
                    }else
                    {
                    counter = counter + 1;
                    cs = cs + posterior3[(counter-1)];     
                    }
               }
          Znew(k,j) = counter;
          }
     }



// Update the element at time ntime
lambdacurrent = lambda(ntimeminus , _);
     for(int k = 0; k < nsites; k++) 
     {
     // Compute the Poisson likelihood component
     like1 = Y(k,ntimeminus) * lambdacurrent - exp(lambdacurrent) * Offset(k,ntimeminus);

     // Compute the prior part
     //prior1 = exp(-alpha * abs(SS - Znew(k,(ntime-2))) - abs(SS - Gstar));
     prior1 = exp(-alpha * pow(SS - Znew(k,(ntime-2)),2) - pow(SS - Gstar,2));
     prior2 = log(prior1 / sum(prior1));    

     // Compute the posterior probability
     posterior1 = like1 + prior2;
     posterior2 = posterior1 - max(posterior1);
     posterior3 = exp(posterior2) / sum(exp(posterior2));
     u = runif(1)[0];          
     cs = posterior3[0];
     test = 1;
     counter = 1;
          while(test==1)
          {
               if(cs>u)
               {
               test = 0;     
               }else
               {
               counter = counter + 1;
               cs = cs + posterior3[(counter-1)];     
               }
          }
     Znew(k,ntimeminus) = counter;
     }




// Return the updated values of Z 
return Znew;
}



// [[Rcpp::export]]
double alphaupdate(IntegerMatrix Z, const int nsites, NumericMatrix logratio, 
     const int ntime)
{
//Create the variable to save the log probability
double logprob=0;
int row, col;
    
//  Compute the probability
     for(int j = 1; j < ntime; j++)
     {
          for(int l = 0; l < nsites; l++) 
          {
           row = Z(l, j-1) - 1;
           col = Z(l, j) - 1;
           logprob = logprob + logratio(row,col);
          }
     }
     
// Return the log probability 
return logprob;
}




// [[Rcpp::export]]
List Xupdate(NumericMatrix Wspace, NumericMatrix Y, NumericMatrix offset,
     const int K, const int N, NumericMatrix X, double tau2, double proposalsdX)
{
//Create new objects
int accept=0;
NumericMatrix Xnew(K,N);
NumericVector Xcurrent(K), convcurrent(K), convadjust(K);
NumericVector Ytemp(K), offsettemp(K);
double Xprop;
double prob1, prob2, prob3, acceptance;
            
Xnew = X;
   
//  Update each random effect in turn
     for(int t = 0; t < N; t++)
     {
     // Specify quantities that do not change
     Ytemp = Y( _, t);
     offsettemp = offset( _, t);

     // Specify the weights for time t
     Xcurrent = Xnew( _, t);   
          for(int r = 0; r < K; r++)
          {
          convcurrent[r] = sum(Wspace(r, _) * Xcurrent);
          }

          for(int j = 0; j < K; j++)
          {
          // Propose a new value
          Xprop = rnorm(1, Xcurrent[j], proposalsdX)[0];
          convadjust = Wspace( _, j) * (Xprop - Xcurrent[j]);
         
          // Compute the acceptance probability
          prob1 = sum(Ytemp * convadjust);
          prob2 = sum(offsettemp * (exp(convcurrent) - exp(convcurrent  + convadjust))); 
          prob3 = (pow(Xcurrent[j],2) - pow(Xprop,2)) / (2 * tau2); 
          acceptance = exp(prob1 + prob2 + prob3);
               if(runif(1)[0] <= acceptance) 
               {
               Xnew(j,t) = Xprop;
               Xcurrent[j] = Xprop;
               convcurrent = convcurrent + convadjust;
               accept = accept + 1;
               }
               else
               { 
               }
          }
     }

// Return the results
List out(2);
out[0] = Xnew;
out[1] = accept;
return out;
}

