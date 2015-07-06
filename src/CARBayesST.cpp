#include <Rcpp.h>
using namespace Rcpp;

// This file contains the following functions:

// linpredcompute - computing the linear predictor for covariates.
// quadform - computing quadratic forms phi %*% Q %*% theta.
// poissoncarupdate - for updating spatial CAR effects based on a poisson likelihood.
// poissonindepupdate - for updating independent effects based on a poisson likelihood.
// poissonbetaupdate - for updating covariate effects based on a poisson likelihood.
// binomialbetaupdate - for updating covariate effects based on a binomial likelihood.
// binomialindepupdate - for updating independent effects based on a binomial likelihood.
// binomialcarupdate - for updating spatial CAR effects based on a binomial likelihood.
// gaussiancarupdate - for updating spatial CAR effects based on a gaussian likelihood.
// poissonarcarupdate - for updating spatio-temporal ARCAR effects based on a poisson likelihood.
// gammaquadformcompute - for computing the sum of quadratic forms for updating gamma in the ST.CARar model.
// tauquadformcompute - for computing the sum of quadratic forms for updating tau2 in the ST.CARar model.
// biomialarcarupdate - for updating spatio-temporal ARCAR effects based on a binomial likelihood.
// gaussianarcarupdate - for updating spatio-temporal ARCAR effects based on a gaussian likelihood.
// qform - computing a quadratic form from a triplet phi %*% Q %*% phi.
// qform_asym - computing a quadratic form from a triplet of the form phi %*% Q %*% theta.
// qformSPACETIME - computing a quadratic form in space nad time.
// SPTICARphiVarb - update space time ARCAR random effects from a varying non-binary W and a Poisson likelihood.
// updatetripList - update the triplet form of W based on a new estimate of W.
// SPTICARphiBinomial - update space time ARCAR random effects from a varying non-binary W and a binomial likelihood.
// SPTICARphiGaussian - update space time ARCAR random effects from a varying non-binary W and a gaussian likelihood.
// qform_difference_ST - this function works out the difference between two quadratic forms (x'Qx - y'Qy) where Q is a kronecker product of Q.space and Q.time
// qform_ST            - this function works out the quadratic forms phi'Qphi where Q is a kronecker product of Q.space and Q.time
// qform_ST_asym       - this function works out the quadratic forms phi1'Qphi2 where Q is a kronecker product of Q.space and Q.time
// update_Qtime        - this function updates the temporal precision matrix given an update of alpha, the temporal dependence par
// updatetriplets_rho  - updates the triplet form of Q.space given an update of the leroux parameter, rho
// updatetripList2     - updates the triplet form of Q.space given an update of the adjacency parameters, v
// Zupdatesqbin - updates the Z allocation parameters in the binomial cluster model.
// Zupdatesqpoi - updates the Z allocation parameters in the poisson cluster model.
// Zupdatesqgau - updates the Z allocation parameters in the gaussian cluster model.

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
double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                    NumericVector phi, NumericVector theta, double rho)
{
// Compute a quadratic form for the random effects
// Create new objects 
double tau2_posteriorscale;
double tau2_quadform = 0, tau2_phisq = 0;
int row, col;
   
   
// Compute the off diagonal elements of the quadratic form
     for(int l = 0; l < n_triplet; l++)
     {
     row = Wtriplet(l,0) - 1;
     col = Wtriplet(l,1) - 1;
     tau2_quadform = tau2_quadform + phi[(Wtriplet(l,0) - 1)] * theta[(Wtriplet(l,1) - 1)] * Wtriplet(l,2); 
     }
 
 
 // Compute the diagonal elements of the quadratic form          
     for(int l = 0; l < nsites; l++)
     {
     tau2_phisq = tau2_phisq + phi[l] * theta[l] * (rho * Wtripletsum[l] + 1 - rho);    
     }
           
     
// Compute the quadratic form
tau2_posteriorscale = 0.5 * (tau2_phisq - rho * tau2_quadform);

 
// Return the simulated value
return tau2_posteriorscale;
}



// [[Rcpp::export]]
List poissoncarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, const NumericMatrix y, const double phi_tune, 
     double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0,rowstart=0, rowend=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi, lpold, lpnew;
NumericVector phinew(nsites);

   
//  Update each random effect in turn
phinew = phi;

     for(int j = 0; j < nsites; j++)
    {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
     priormean = rho * sumphi / priorvardenom; 
     
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];

      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      
      oldlikebit = 0;
      newlikebit = 0;
        for(int i=0; i < ntime; i++)
        {
        lpold = mult_offset[i] * phinew[j] + offset(j, i);
        lpnew = mult_offset[i] * propphi + offset(j, i); 
        oldlikebit = oldlikebit + y(j,i) * lpold - exp(lpold);
        newlikebit = newlikebit + y(j,i) * lpnew - exp(lpnew);
        }       
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
double acceptance, proptheta, lpold, lpnew;
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
      lpold = thetanew[j] + offset[j];
      lpnew = proptheta + offset[j];
      oldlikebit = lpold * y[j]  - exp(lpold);
      newlikebit =  lpnew * y[j]  - exp(lpnew);
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
  NumericVector lp_current(nsites), lp_proposal(nsites);
  
  // Create the log likelihood acceptance probability component
  lp_current  = linpredcompute(X, nsites, p, beta, offset);
  lp_proposal = linpredcompute(X, nsites, p, proposal, offset);     
  for(int j = 0; j < nsites; j++){
    oldlikebit += y[j] * lp_current[j] -  exp(lp_current[j]);
    newlikebit += y[j] * lp_proposal[j] - exp(lp_proposal[j]);
  }
  likebit = newlikebit - oldlikebit;
  
  // Create the prior acceptance component
  for(int j = 0; j < p; j++) priorbit += 0.5 * pow((beta[j]-prior_meanbeta[j]),2) / prior_varbeta[j] - 0.5 * pow((proposal[j]-prior_meanbeta[j]),2) / prior_varbeta[j];
  
  // Compute the acceptance probability and return the value
  acceptance = exp(likebit + priorbit);
  return acceptance;
}




// [[Rcpp::export]]
double binomialbetaupdate(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
               NumericVector proposal, NumericVector offset, NumericVector y, 
               NumericVector failures, NumericVector prior_meanbeta,
               NumericVector prior_varbeta)
{
// Compute the acceptance probability for beta
//Create new objects
double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit=0;
NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);


// Create the log likelihood acceptance probability component
lp_current = linpredcompute(X, nsites, p, beta, offset);
lp_proposal = linpredcompute(X, nsites, p, proposal, offset);     
     for(int j = 0; j < nsites; j++)     
     {
     p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
     p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
     oldlikebit = oldlikebit + y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
     newlikebit = newlikebit + y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
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
List binomialindepupdate(const int nsites, NumericVector theta, double tau2, 
const NumericVector y, const NumericVector failures, const double theta_tune, NumericVector offset)
{
// Update the spatially independent random effects 
//Create new objects
int accept=0;
double acceptance, proptheta, lpold, lpnew, pold, pnew;
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
      lpold = thetanew[j] + offset[j];
      lpnew = proptheta + offset[j];
      pold = exp(lpold) / (1 + exp(lpold));
      pnew = exp(lpnew) / (1 + exp(lpnew));
      
      oldlikebit = y[j] * log(pold)  + failures[j] * log((1-pold));
      newlikebit =  y[j] * log(pnew)  + failures[j] * log((1-pnew));
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
List binomialcarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, const NumericMatrix y,  const NumericMatrix failures, const double phi_tune, 
     double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset)
{
// Update the spatially correlated random effects 
//Create new objects
int accept=0,rowstart=0, rowend=0;
double acceptance, sumphi;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
double priorvardenom, priormean, priorvar;
double propphi, lpold, lpnew, pold, pnew;
NumericVector phinew(nsites);

   
//  Update each random effect in turn
phinew = phi;

     for(int j = 0; j < nsites; j++)
    {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
     priormean = rho * sumphi / priorvardenom; 
     
      // propose a value  
      propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];

      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
      
      oldlikebit = 0;
      newlikebit = 0;
        for(int i=0; i < ntime; i++)
        {
        lpold = mult_offset[i] * phinew[j] + offset(j, i);
        lpnew = mult_offset[i] * propphi + offset(j, i); 
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));        
        
        oldlikebit = oldlikebit + y(j,i) * log(pold) + failures(j,i) * log((1-pold));
        newlikebit = newlikebit + y(j,i) * log(pnew) + failures(j,i) * log((1-pnew));
        }

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
NumericVector gaussiancarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, double nu2, const NumericVector offset, double rho, double ntime)
{
// Update the spatially correlated random effects 
//Create new objects
int rowstart=0, rowend=0;
double sumphi;
double fcmean, fcvar;
double priorvardenom, priormean, priorvar;
NumericVector phinew(nsites);

   
//  Update each random effect in turn
phinew = phi;

     for(int j = 0; j < nsites; j++)
    {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
     priormean = rho * sumphi / priorvardenom; 
     
      // compute the full conditional  
      fcvar = 1 / (1 / priorvar + ntime / nu2);
      fcmean = fcvar * (priormean / priorvar +  offset[j]/ nu2);
      phinew[j] = rnorm(1, fcmean, sqrt(fcvar))[0];
     }
        
return phinew;
}




// [[Rcpp::export]]
List poissonarcarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, const int ntime,
          NumericMatrix phi, double tau2, double gamma, double rho, 
          const NumericMatrix ymat, const double phi_tune, NumericMatrix offset,
          NumericVector denoffset)
{    
///////////////////////////////////////////    
// Specify variables needed in the function
///////////////////////////////////////////
double temp, priormean, priormeantemp1, priormeantemp2, priorvar, priorvardenom, acceptance;
double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew; 
NumericMatrix phinew(nsites,ntime);
phinew = phi;
int row, rowstart, rowend, accept=0;


//////////////////////////////////////////////
// Update the random effects at time 1 in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     priorvardenom = denoffset[j] * (1 + pow(gamma,2));
     priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * phinew(j,1);
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * phinew(row,1));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;

    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
    lpold = phinew(j,0) + offset(j, 0);
    lpnew = propphi + offset(j, 0); 
    oldlikebit = ymat(j,0) * lpold - exp(lpold);
    newlikebit = ymat(j,0) * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew(j,0) = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
     }
    
    
    
//////////////////////////////////////////////////////
// Update the random effects at times 2 to N-1 in turn
//////////////////////////////////////////////////////
     for(int t = 1; t < (ntime-1); t++)
     {
        for(int j = 0; j < nsites; j++)
        {
        // calculate prior mean and variance
        priorvardenom = denoffset[j] * (1 + pow(gamma,2));
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;
    
        // Propose a value and calculate the acceptance probability
        propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
        lpold = phinew(j,t) + offset(j, t);
        lpnew = propphi + offset(j, t); 
        oldlikebit = ymat(j,t) * lpold - exp(lpold);
        newlikebit = ymat(j,t) * lpnew - exp(lpnew);
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew(j,t) = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
        }
     }
    
    
    
//////////////////////////////////////////////
// Update the random effects at time N in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     priorvardenom = denoffset[j];
     priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
    lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
    lpnew = propphi + offset(j, (ntime-1)); 
    oldlikebit = ymat(j,(ntime-1)) * lpold - exp(lpold);
    newlikebit = ymat(j,(ntime-1)) * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew(j,(ntime-1)) = propphi;
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
List gammaquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, 
    const int n_triplet, const int nsites, const int ntime, NumericMatrix phi, double rho)
{    
NumericVector phi_t(nsites), phi_tminus1(nsites);
double num=0, den=0;

// Compute the sum of quadratic forms for updating gamma
    for(int t = 1; t < ntime; t++)
    {
    phi_t = phi(_,t);
    phi_tminus1 = phi(_,(t-1));    
    num = num + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_t, phi_tminus1, rho);
    den = den + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_tminus1, phi_tminus1, rho);    
    }


List out(2);
out[0] = num;
out[1] = den;

return out;
}





// [[Rcpp::export]]
double tauquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
    const int nsites, const int ntime, NumericMatrix phi, double rho, double gamma)
{    
NumericVector temp(nsites);
double num=0;

// Compute the sum of quadratic forms for updating tau
temp = phi(_,0);
num = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);

    for(int t = 1; t < ntime; t++)
    {
    temp = phi(_,t) - gamma * phi(_,(t-1));  
    num = num + quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
    }


return num;
}






// [[Rcpp::export]]
List binomialarcarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, const int ntime,
          NumericMatrix phi, double tau2, double gamma, double rho, 
          const NumericMatrix ymat, const NumericMatrix failuresmat,
          const double phi_tune, NumericMatrix offset,NumericVector denoffset)
{    
///////////////////////////////////////////    
// Specify variables needed in the function
///////////////////////////////////////////
double temp, priormean, priormeantemp1, priormeantemp2, priorvar, priorvardenom, acceptance;
double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, pold, pnew; 
NumericMatrix phinew(nsites,ntime);
phinew = phi;
int row, rowstart, rowend, accept=0;


//////////////////////////////////////////////
// Update the random effects at time 1 in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(gamma,2));
    priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,1));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * (phinew(row,1)));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;    
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
    lpold = phinew(j,0) + offset(j, 0);
    lpnew = propphi + offset(j, 0); 
    pold = exp(lpold) / (1 + exp(lpold));
    pnew = exp(lpnew) / (1 + exp(lpnew));        
    oldlikebit = ymat(j,0) * log(pold) + failuresmat(j,0) * log((1-pold));
    newlikebit = ymat(j,0) * log(pnew) + failuresmat(j,0) * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew(j,0) = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
     }
    
    
    
//////////////////////////////////////////////////////
// Update the random effects at times 2 to N-1 in turn
//////////////////////////////////////////////////////
     for(int t = 1; t < (ntime-1); t++)
     {
        for(int j = 0; j < nsites; j++)
        {
        // calculate prior mean and variance
        priorvardenom = denoffset[j] * (1 + pow(gamma,2));
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
    
        // Propose a value and calculate the acceptance probability
        propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
        lpold = phinew(j,t) + offset(j, t);
        lpnew = propphi + offset(j, t); 
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));        
        oldlikebit = ymat(j,t) * log(pold) + failuresmat(j,t) * log((1-pold));
        newlikebit = ymat(j,t) * log(pnew) + failuresmat(j,t) * log((1-pnew));
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew(j,t) = propphi;
          accept = accept + 1;
          }
          else
          { 
          }
        }
     }
    
    
    
//////////////////////////////////////////////
// Update the random effects at time N in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     priorvardenom = denoffset[j];
     priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
    lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
    lpnew = propphi + offset(j, (ntime-1)); 
    pold = exp(lpold) / (1 + exp(lpold));
    pnew = exp(lpnew) / (1 + exp(lpnew));        
    oldlikebit = ymat(j,(ntime-1)) * log(pold) + failuresmat(j,(ntime-1)) * log((1-pold));
    newlikebit = ymat(j,(ntime-1)) * log(pnew) + failuresmat(j,(ntime-1)) * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
          if(runif(1)[0] <= acceptance) 
          {
          phinew(j,(ntime-1)) = propphi;
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
NumericMatrix gaussianarcarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, const int ntime,
          NumericMatrix phi, double tau2, double nu2, double gamma, double rho, 
          NumericMatrix offset, NumericVector denoffset)
{    
///////////////////////////////////////////    
// Specify variables needed in the function
///////////////////////////////////////////
double temp, priormean, priormeantemp1, priormeantemp2, priorvar,priorvardenom;
double propphi, fcvar, fcmean; 
NumericMatrix phinew(nsites,ntime);
phinew = phi;
int row, rowstart, rowend;


//////////////////////////////////////////////
// Update the random effects at time 1 in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(gamma,2));
    priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,1));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * (phinew(row,1)));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
        
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,0)/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,0) = propphi;
    }
    
    
    
//////////////////////////////////////////////////////
// Update the random effects at times 2 to N-1 in turn
//////////////////////////////////////////////////////
     for(int t = 1; t < (ntime-1); t++)
     {
        for(int j = 0; j < nsites; j++)
        {
        // calculate prior mean and variance
        priorvardenom = denoffset[j] * (1 + pow(gamma,2));
        priorvar = tau2 / priorvardenom;
        priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
        priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
         
        // Compute the full conditional and update phi
        fcvar = 1 / (1 / priorvar + 1 / nu2);
        fcmean = fcvar * (priormean / priorvar +  offset(j,t)/ nu2);
        propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
        phinew(j,t) = propphi;
        }
     }
    
    
    
//////////////////////////////////////////////
// Update the random effects at time N in turn
//////////////////////////////////////////////
    for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     priorvardenom = denoffset[j];
     priorvar = tau2 / priorvardenom;
     priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     priormeantemp2 = 0;
          for(int l = rowstart; l < rowend; l++)
          {
          row = Wtriplet(l,1) - 1;
          temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
          priormeantemp2 = priormeantemp2 + temp; 
          }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
        
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,(ntime-1))/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,(ntime-1)) = propphi;
     }
    
return phinew;
}









// [[Rcpp::export]]
double qform(NumericMatrix Qtrip, NumericVector phi){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  for(int i = 0; i < nzero; i++){
    Qform += phi[Qtrip(i, 0) - 1] * Qtrip(i, 2) * phi[Qtrip(i, 1) - 1];
  }
  return(Qform);
}


// [[Rcpp::export]]
double qform_asym(NumericMatrix Qtrip, NumericVector phi1, NumericVector phi2){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  for(int i = 0; i < nzero; i++){
    Qform += phi1[Qtrip(i, 0) - 1] * Qtrip(i, 2) * phi2[Qtrip(i, 1) - 1];
  }
  return(Qform);
}

// [[Rcpp::export]]
double qformSPACETIME(NumericMatrix Qtrip, NumericVector phi, const int ntime, const int nsite){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  int spaceBlock = 0;
  for(int j = 0; j < ntime; j++){
    spaceBlock = j*nsite - 1;
    for(int i = 0; i < nzero; i++){
      Qform += phi[Qtrip(i, 0) + spaceBlock] * Qtrip(i, 2) * phi[Qtrip(i, 1) + spaceBlock];
    }
  }
  return(Qform);
}




// [[Rcpp::export]]
List SPTICARphiBinomial(NumericMatrix W, const int nsites, const int ntimes, 
                          NumericVector phi, NumericVector nneighbours, double tau,
                          const NumericVector y, double alpha, NumericVector XB, 
                          const double phiVarb_tune, NumericVector trials){
  
  // stuff associated with phiVarb update
  int    n;
  double theta_prop;
  double theta;
  double l_prob_new;
  double l_prob_old;
  double sumphiVarb;
  double priorvardenom;
  double priormean;
  double priorvar;
  double phi_new;
  double futuresumphiVarb;
  double pastsumphiVarb;
  double acceptance;
  NumericVector accept(4);
  double asqOne = 1 + pow(alpha, 2); 
  NumericVector phiVarb = clone(phi);
  
  //  outer loop for each random effect
  for(int i = 0; i < ntimes; i++) { 
    for(int j = 0; j < nsites; j++){
      sumphiVarb = 0;
      futuresumphiVarb = 0;
      pastsumphiVarb = 0;
      n = (i*nsites) + j;
      //  calulate the sum of the neighbours of j at time i
      for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
      // calculate prior variance
      priorvardenom = nneighbours[j];
      // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
        if(ntimes > 1) {
          if((0 < i) && (i < (ntimes-1))){
            priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
            //calculate the sum of the neighbours of j in the past and the futuren
            for(int l = 0; l < nsites; l++)  {
              futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
              pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
            }
            priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
            priorvar = tau/(priorvardenom*asqOne);
          } else if(i == 0) {
            priormean = phiVarb[((i+1)*nsites) + j];
            //calculate the sum of the neighbours of j in the future
            for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
            priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
            priorvar =  tau/(priorvardenom*asqOne);
          } else if(i == (ntimes-1)) {
            priormean = phiVarb[((i-1)*nsites) + j];
            //calculate the sum of the neighbours of j in the past
            for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
            priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
            priorvar = tau/(priorvardenom);
          }
        } else if(ntimes == 1){
          priorvar  = tau/priorvardenom;
          priormean = sumphiVarb/priorvardenom; 
        }
      
      // get the theta parameter in [0,1] of the binomial under proposal and old 
      phi_new      = rnorm(1, phiVarb[n], sqrt(priorvar*phiVarb_tune))[0];    
      theta        = 1 + exp(-XB[n] - phiVarb[n]);
      theta_prop   = 1 + exp(-XB[n] - phi_new);
      l_prob_old   = -y[n]*log(theta) + (trials[n] - y[n])*log(1 - (1/theta)) - (0.5/priorvar) * pow(phiVarb[n] - priormean, 2);
      l_prob_new   = -y[n]*log(theta_prop) + (trials[n] - y[n])*log(1 - (1/theta_prop)) - (0.5/priorvar) * pow(phi_new - priormean, 2);
      acceptance   = exp(l_prob_new - l_prob_old);
      if(acceptance >= 1){
        phiVarb[n] = phi_new;
        accept[1]++;
      } else {
        if(runif(1)[0] <= acceptance) {
          phiVarb[n] = phi_new;
          accept[1]++;
        }
      }
    }
  }
  List out(2);
  out[0] = accept;
  out[1] = phiVarb;
  return out;
}






// [[Rcpp::export]]
List SPTICARphiGaussian(NumericMatrix W, const int nsites, const int ntimes, 
                        NumericVector phi, NumericVector nneighbours, double tau, double lik_var,
                        const NumericVector y, double alpha, NumericVector XB){
  
  // stuff associated with phiVarb update
  int    n;
  double meanFullCon;
  double varFullCon;
  double sumphiVarb;
  double priorvardenom;
  double priormean;
  double priorvar;
  double futuresumphiVarb;
  double pastsumphiVarb;
  double asqOne = 1 + pow(alpha, 2); 
  NumericVector phiVarb = clone(phi);
  
  //  outer loop for each random effect
  for(int i = 0; i < ntimes; i++) { 
    for(int j = 0; j < nsites; j++){
      sumphiVarb = 0;
      futuresumphiVarb = 0;
      pastsumphiVarb = 0;
      n = (i*nsites) + j;
      //  calulate the sum of the neighbours of j at time i
      for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
      // calculate prior variance
      priorvardenom = nneighbours[j];
      // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
        if(ntimes > 1) {
          if((0 < i) && (i < (ntimes-1))){
            priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
            //calculate the sum of the neighbours of j in the past and the futuren
            for(int l = 0; l < nsites; l++)  {
              futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
              pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
            }
            priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
            priorvar = tau/(priorvardenom*asqOne);
          } else if(i == 0) {
            priormean = phiVarb[((i+1)*nsites) + j];
            //calculate the sum of the neighbours of j in the future
            for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
            priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
            priorvar =  tau/(priorvardenom*asqOne);
          } else if(i == (ntimes-1)) {
            priormean = phiVarb[((i-1)*nsites) + j];
            //calculate the sum of the neighbours of j in the past
            for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
            priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
            priorvar = tau/(priorvardenom);
          }
        } else if(ntimes == 1){
          priorvar  = tau/priorvardenom;
          priormean = sumphiVarb/priorvardenom; 
        }
      
      // get the mean and variance of the full conditional distribution
      varFullCon  = 1/((1/priorvar) + (1/lik_var));
      meanFullCon = ((priormean/priorvar) + (y[n] - XB[n])/lik_var)*varFullCon;
      phiVarb[n]  = rnorm(1, meanFullCon, sqrt(varFullCon))[0];    
    }
  }
  List out(2);
  out[0] = phiVarb;
  return out;
}


// [[Rcpp::export]]
List SPTICARphiVarb(NumericMatrix W, const int nsites, const int ntimes, 
                    NumericVector phiVarb, NumericVector nneighbours, double tau, 
                    const NumericVector y, const NumericVector E, 
                    const double phiVarb_tune, double alpha, NumericVector XB){
  
  //stuff associated with model objects
  int n;
  int k = ntimes*nsites;
  
  //General MCMC things
  NumericVector accept(4);
  double acceptance;
  
  // stuff associated with phiVarb update
  double newpriorbit;
  double oldpriorbit;
  double sumphiVarb;
  double gubbins;
  double priorvardenom, priormean;
  double priorvar;
  double propphiVarb;
  double futuresumphiVarb, pastsumphiVarb;
  double asqOne = 1 + pow(alpha, 2); 
  
  // stuff associated with tau update
  NumericVector arc(k);
  //Function rtrunc("rtrunc");
  
  //  outer loop for each random effect
  for(int i = 0; i < ntimes; i++) { 
    for(int j = 0; j < nsites; j++){
      sumphiVarb = 0;
      futuresumphiVarb = 0;
      pastsumphiVarb = 0;
      n = (i*nsites) + j;
      //  calulate the sum of the neighbours of j at time i
      for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
      // calculate prior variance
      priorvardenom = nneighbours[j];
      // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
      if(ntimes > 1) {
        if((0 < i) && (i < (ntimes-1))){
          priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
          //calculate the sum of the neighbours of j in the past and the futuren
          for(int l = 0; l < nsites; l++)  {
            futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
            pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
          }
          priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
          priorvar = tau/(priorvardenom*asqOne);
        } else if(i == 0) {
          priormean = phiVarb[((i+1)*nsites) + j];
          //calculate the sum of the neighbours of j in the future
          for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
          priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
          priorvar =  tau/(priorvardenom*asqOne);
        } else if(i == (ntimes-1)) {
          priormean = phiVarb[((i-1)*nsites) + j];
          //calculate the sum of the neighbours of j in the past
          for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
          priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
          priorvar = tau/(priorvardenom);
        }
      } else if(ntimes == 1){
        priorvar = tau/priorvardenom;
        priormean = 1*sumphiVarb/priorvardenom; 
      }
      
      // propose a value and accept or reject it 
      propphiVarb = rnorm(1, phiVarb[n], sqrt(priorvar*phiVarb_tune))[0];    
      newpriorbit = (0.5/priorvar) * pow(propphiVarb - priormean, 2); 
      oldpriorbit = (0.5/priorvar) * pow(phiVarb[n] - priormean, 2);
      gubbins = exp(E[n] + XB[n])*(exp(propphiVarb) - exp(phiVarb[n]));
      acceptance = exp(y[n]*(propphiVarb - phiVarb[n]) - newpriorbit + oldpriorbit - gubbins);
      arc[n] = acceptance;
      if(acceptance >= 1){
        phiVarb[n] = propphiVarb;
        accept[1]++;
      } else {
        if(runif(1)[0] <= acceptance) {
          phiVarb[n] = propphiVarb;
          accept[1]++;
        }
      }
    }
  }
  List out(3);
  out[0] = accept;
  out[1] = phiVarb;
  out[2] = arc;
  return out;
}


// [[Rcpp::export]]
double qform_difference_ST(NumericMatrix Qtrip, NumericMatrix Qtime, NumericVector phi, int nsites){ 
  int nrowSpace = Qtrip.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    if( !(Qtrip(i, 2) == 0) ){
      spRow  = Qtrip(i, 0);
      spCol  = Qtrip(i, 1);
      for(int j = 0; j < nrowTime; j++){
        tiRow  = Qtime(j, 0);
        tiCol  = Qtime(j, 1);
        stRow  = (tiRow - 1)*nsites + spRow;
        stCol  = (tiCol - 1)*nsites + spCol;
        Qform += phi[stCol - 1] * Qtrip(i, 2) * Qtime(j, 2) * phi[stRow - 1];
      }
    }
  }
  return(Qform);
}

// [[Rcpp::export]]
double qform_ST(NumericMatrix Qspace, NumericMatrix Qtime, NumericVector phi, int nsites){ 
  int nrowSpace = Qspace.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    spRow  = Qspace(i, 0);
    spCol  = Qspace(i, 1);
    for(int j = 0; j < nrowTime; j++){
      tiRow  = Qtime(j, 0);
      tiCol  = Qtime(j, 1);
      stRow  = (tiRow - 1)*nsites + spRow;
      stCol  = (tiCol - 1)*nsites + spCol;
      Qform += phi[stCol - 1] * Qspace(i, 2) * Qtime(j, 2) * phi[stRow - 1];
    }
  }
  return(Qform);
}

// [[Rcpp::export]]
double qform_ST_asym(NumericMatrix Qspace, NumericMatrix Qtime, NumericVector phi1, NumericVector phi2, int nsites){ 
  int nrowSpace = Qspace.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    spRow  = Qspace(i, 0);
    spCol  = Qspace(i, 1);
    for(int j = 0; j < nrowTime; j++){
      tiRow  = Qtime(j, 0);
      tiCol  = Qtime(j, 1);
      stRow  = (tiRow - 1)*nsites + spRow;
      stCol  = (tiCol - 1)*nsites + spCol;
      Qform += phi1[stCol - 1] * Qspace(i, 2) * Qtime(j, 2) * phi2[stRow - 1];
    }
  }
  return(Qform);
}

// [[Rcpp::export]]
NumericMatrix update_Qtime(NumericMatrix Qtime, double alpha, int rowNumberLastDiag){
  int nr = Qtime.nrow();
  double alphasq = alpha*alpha;
  NumericMatrix Qtime_new = clone(Qtime);
  for(int i = 0; i < nr; i++){
    if(Qtime(i, 0) == Qtime(i, 1))  Qtime_new(i, 2) = 1 + alphasq;
    if(!(Qtime(i, 0) == Qtime(i, 1)) ) Qtime_new(i, 2) = -alpha;
  }
  Qtime_new(rowNumberLastDiag,2) = 1;
  return Qtime_new;
}


// [[Rcpp::export]]
NumericMatrix updatetriplets_rho(NumericMatrix trips, int nsites, double rho_old, double rho_new, double fixedridge){    
  //   create a clone of the triplet matrix that will be used for output              
  NumericMatrix tripsnew = clone(trips);   
  int rows               = tripsnew.nrow();
  double one_rho_old     = 1 - rho_old;
  double one_rho_new     = 1 - rho_new;
  //   loop over all of the diagonal elements first.  To update the diagonal elements, 
  //   we first substract (1 - rho_old) and the fixed ridge and divide by rho_old.  Then multiply by rho_new
  //   and add 1 - rho_new.
  for(int i = 0; i < nsites; i++) {
    tripsnew(i, 2) = ((trips(i, 2) - one_rho_old - fixedridge)/rho_old)*rho_new + one_rho_new + fixedridge;
  }
  //   loop over all of the off-diagonal elements.  These are updated by simply 
  //   dividing by rho_old and multiplying by rho_new
  for(int j = nsites; j < rows; j++) {
    tripsnew(j, 2)     = (trips(j, 2)/rho_old)*rho_new;
  }
  //   output the updated triplets
  return tripsnew;
}



// [[Rcpp::export]]
List updatetripList2(NumericMatrix trips, NumericVector vold, 
                     NumericVector vnew, const int nedges, int nsites, IntegerVector block, int block_length, double rho, 
                     double fixedridge){                   
  //   create a clone of the triplet matrix                   
  NumericMatrix temporary = clone(trips); 
  NumericMatrix difference = clone(trips);
  for(int l = 0; l < temporary.nrow(); l++) difference(l, 2) = 0;
  //   stuff needed for intermediate calculations
  double oldoffdiag_binary, newoffdiag_binary;
  double newoffdiag_logit,  oldoffdiag_logit;
  int    diag_inds_1,       diag_inds_2, i;
  double olddiag_binary_1,  olddiag_binary_2;
  double newdiag_binary_1,  newdiag_binary_2;
  //   loop over all edges, stop only at elements where vold and vnew differ
  //   then perform and update of the corresponding triplets
  for(int j = 0; j < block_length; j++) {  
    i = block[j] - 1;
    //       this is the old off diagonal element (-ve to make +ve)
    oldoffdiag_binary = -temporary(i + nsites, 2)/rho;
    //       old and new v elements
    oldoffdiag_logit  = vold[i];
    newoffdiag_logit  = vnew[i];
    //       convert new v onto [0,1] scale
    newoffdiag_binary = -(1/(1 + exp(-newoffdiag_logit)));
    //       replace triplets with new v
    difference(i + nsites, 2) = rho*(temporary(i + nsites, 2) - newoffdiag_binary);
    temporary(i + nsites, 2)  = rho*newoffdiag_binary;
    difference(i + nsites + nedges, 2) = rho*(temporary(i + nsites + nedges, 2) - newoffdiag_binary);
    temporary(i + nsites + nedges, 2)  = rho*(newoffdiag_binary);
    
    //       now need to find x and y coords of these offdiags to update diags
    diag_inds_1 = temporary(i + nsites, 0) - 1;
    diag_inds_2 = temporary(i + nsites, 1) - 1;
    //       get the old binary elements
    olddiag_binary_1 = (temporary(diag_inds_1, 2) - fixedridge - (1 - rho))/rho;
    olddiag_binary_2 = (temporary(diag_inds_2, 2) - fixedridge - (1 - rho))/rho;
    //       calculate and replace with new ones
    newdiag_binary_1 = (olddiag_binary_1 - oldoffdiag_binary)  - newoffdiag_binary;
    newdiag_binary_2 = (olddiag_binary_2 - oldoffdiag_binary)  - newoffdiag_binary;
    temporary(diag_inds_1, 2)  = (newdiag_binary_1*rho) + fixedridge + (1 - rho);
    temporary(diag_inds_2, 2)  = (newdiag_binary_2*rho) + fixedridge + (1 - rho);
    difference(diag_inds_1, 2) = trips(diag_inds_1, 2) -  temporary(diag_inds_1, 2);
    difference(diag_inds_2, 2) = trips(diag_inds_2, 2) -  temporary(diag_inds_2, 2);
  }
  //   output to a list where the first element is the updated diagonal triplets
  //   second element is the offdiagonal triplets
  List out(2);
  out[0] = temporary;
  out[1] = difference;
  return out;
}





// [[Rcpp::export]]
NumericMatrix Zupdatesqbin(NumericMatrix Z, NumericMatrix Offset, NumericMatrix Y, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar, NumericMatrix failures)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G), lp(G), prob(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   



// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    lp = Offset(k,0) + lambda;
    prob = exp(lp) / (1 + exp(lp));
    like1 = Y(k,0) * log(prob) + failures(k,0) * log((1 - prob));
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
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
    Z(k,0) = counter;
    }




//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        lp = Offset(k,j) + lambda;
        prob = exp(lp) / (1 + exp(lp));
        like1 = Y(k,j) * log(prob) + failures(k,j) * log((1 - prob));
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
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
        Z(k,j) = counter;
    }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    lp = Offset(k,ntimeminus) + lambda;
    prob = exp(lp) / (1 + exp(lp));
    like1 = Y(k,ntimeminus) * log(prob) + failures(k,ntimeminus) * log((1 - prob));
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
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
    Z(k,ntimeminus) = counter;
    }


return Z;
}








// [[Rcpp::export]]
NumericMatrix Zupdatesqpoi(NumericMatrix Z, NumericMatrix Offset, NumericMatrix Y, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   


// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = Y(k,0) * lambda - exp(lambda) * Offset(k,0);
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
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
    Z(k,0) = counter;
    }


//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        like1 = Y(k,j) * lambda - exp(lambda) * Offset(k,j);
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
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
        Z(k,j) = counter;
        }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = Y(k,ntimeminus) * lambda - exp(lambda) * Offset(k,ntimeminus);
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
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
    Z(k,ntimeminus) = counter;
    }


return Z;
}






// [[Rcpp::export]]
NumericMatrix Zupdatesqgau(NumericMatrix Z, NumericMatrix Offset, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar, const double nu2)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   



// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = -pow((Offset(k,0) - lambda),2) / (2*nu2);
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
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
    Z(k,0) = counter;
    }


//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        like1 = -pow((Offset(k,j) - lambda),2) / (2*nu2);
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
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
        Z(k,j) = counter;
        }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = -pow((Offset(k,ntimeminus) - lambda),2) / (2*nu2);
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
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
    Z(k,ntimeminus) = counter;
    }




return Z;
}
