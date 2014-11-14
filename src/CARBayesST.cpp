#include <Rcpp.h>
using namespace Rcpp;

// This file contains the following functions:

// linpredcompute - computing the linear predictor for covariates.
// quadform - computing quadratic forms phi %*% Q %*% theta.
// poissoncarupdate - for updating spatial CAR effects based on a poisson likelihood.
// poissonindepupdate - for updating independent effects based on a poisson likelihood.
// poissonbetaupdate - for updating covariate effects based on a poisson likelihood.
// poissonarcarupdate - for updating spatio-temporal ARCAR effects based on a poisson likelihood.
// Zupdate - for updating the allocation variables Z in the cluster models.
// norm - for computing the normalising constant for the parameters controlling the Z prior.
// qform - computing a quadratic form from a triplet phi %*% Q %*% phi.
// qform_asym - computing a quadratic form from a triplet of the form phi %*% Q %*% theta.
// qformSPACETIME - computing a quadratic form in space nad time.
// SPTICARphiVarb - update space time ARCAR random effects from a varying non-binary W.
// updatetripList - update the triplet form of W based on a new estimate of W.


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
                NumericVector phi, NumericVector theta, NumericVector nneighbours, 
                double diagonal, double offdiagonal)
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
    tau2_quadform = tau2_quadform + phi[row] * theta[col]; 
  }
  
  
  // Compute the diagonal elements of the quadratic form          
  for(int l = 0; l < nsites; l++)
  {
    tau2_phisq =  tau2_phisq + phi[l] * theta[l] * (diagonal * nneighbours[l] + 1 - diagonal);
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
List poissonarcarupdate(List W_list, const int nsites, const int ntime,
          NumericMatrix phi, double tau2, double gamma, double rho, 
          const NumericMatrix ymat, const double phi_tune, NumericMatrix offset,
          NumericVector denoffset)
{    
// Update the spatially correlated random effects 
double sumphi1, sumphi2, sumphi3, priormeantemp1, priormeantemp2, priormeantemp3, priorvartemp1;
double priorvar, priormean, propphi, acceptance;
double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
NumericMatrix phinew(nsites,ntime);
phinew = phi;
int accept = 0;

// Update the random effects at time 1 in turn
     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     IntegerVector neighbourvec = W_list[j];
     int m = neighbourvec.size();
     sumphi1 = 0;
     sumphi2 = 0;
          for(int l = 0; l < m; l++)
          {
          sumphi1 += phinew((neighbourvec[l]-1),0);
          sumphi2 += phinew((neighbourvec[l]-1),1);
          }      
      priormeantemp1 = (1 + pow(gamma,2)) * rho * sumphi1;
      priormeantemp2 = gamma * (denoffset[j] * phi(j,1)  - rho * sumphi2);
      priorvartemp1 = denoffset[j] * (1 + pow(gamma,2));
      priormean = (priormeantemp1 + priormeantemp2) / priorvartemp1;
      priorvar = tau2 / priorvartemp1;      
      
      // propose a value  
      propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
      oldlikebit = phinew(j,0) * ymat(j,0)  - exp(phinew(j,0)) * offset(j,0);
      newlikebit =  propphi * ymat(j,0)  - exp(propphi) * offset(j,0);
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




// Update the random effects at times 2 to ntime-1 in turn
     for(int t = 1; t < (ntime-1); t++)
     {
          for(int j = 0; j < nsites; j++)
          {
          // calculate prior mean and variance
          IntegerVector neighbourvec = W_list[j];
          int m = neighbourvec.size();
          sumphi1 = 0;
          sumphi2 = 0;
          sumphi3 = 0;
               for(int l = 0; l < m; l++)
               {
               sumphi1 += phinew((neighbourvec[l]-1),t);
               sumphi2 += phinew((neighbourvec[l]-1),(t-1));
               sumphi3 += phinew((neighbourvec[l]-1),(t+1));
               }      
          priormeantemp1 = (1 + pow(gamma,2)) * rho * sumphi1;
          priormeantemp2 = gamma * (denoffset[j] * phi(j,(t-1))  - rho * sumphi2);
          priormeantemp3 = gamma * (denoffset[j] * phi(j,(t+1))  - rho * sumphi3);
          priorvartemp1 = denoffset[j] * (1 + pow(gamma,2));
          priormean = (priormeantemp1 + priormeantemp2 + priormeantemp3) / priorvartemp1;
          priorvar = tau2 / priorvartemp1;      
      
          // propose a value  
          propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      
          // Accept or reject it
          newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
          oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
          oldlikebit = phinew(j,t) * ymat(j,t)  - exp(phinew(j,t)) * offset(j,t);
          newlikebit =  propphi * ymat(j,t)  - exp(propphi) * offset(j,t);
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
          


// Update the random effects at time ntime in turn
     for(int j = 0; j < nsites; j++)
     {
     // calculate prior mean and variance
     IntegerVector neighbourvec = W_list[j];
     int m = neighbourvec.size();
     sumphi1 = 0;
     sumphi2 = 0;
          for(int l = 0; l < m; l++)
          {
          sumphi1 += phinew((neighbourvec[l]-1),(ntime-1));
          sumphi2 += phinew((neighbourvec[l]-1),(ntime-2));
          }      
      priormeantemp1 = rho * sumphi1;
      priormeantemp2 = gamma * (denoffset[j] * phi(j,(ntime-2))  - rho * sumphi2);
      priorvartemp1 = denoffset[j];
      priormean = (priormeantemp1 + priormeantemp2) / priorvartemp1;
      priorvar = tau2 / priorvartemp1;      
      
      // propose a value  
      propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
      
      // Accept or reject it
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
      oldlikebit = phinew(j,(ntime-1)) * ymat(j,(ntime-1))  - exp(phinew(j,(ntime-1))) * offset(j,(ntime-1));
      newlikebit =  propphi * ymat(j,(ntime-1))  - exp(propphi) * offset(j,(ntime-1));
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
NumericMatrix Zupdate(NumericMatrix Z, NumericMatrix Offset, NumericMatrix Y, const double alpha, 
NumericMatrix lambda, const int nsites, const int ntime, const int G, NumericVector SS,double Gstar, double delta)
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
     prior1 = exp(-alpha * pow(SS - Znew(k,1),2) - delta * pow(SS - Gstar,2));
     prior2 = log(prior1 / sum(prior1));    
     prior3 = exp(- delta * pow(SS - Gstar,2));
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
          prior1 = exp(-alpha * pow(SS - Znew(k,(j-1)),2) - delta * pow(SS - Gstar,2));
          prior2 = log(prior1 / sum(prior1));    
          prior3 = exp(-alpha * pow(SS - Znew(k,(j+1)),2) - delta * pow(SS - Gstar,2));
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
     prior1 = exp(-alpha * pow(SS - Znew(k,(ntime-2)),2) - delta * pow(SS - Gstar,2));
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
NumericVector norm(NumericVector Z, double G, double Gstar, double alpha, 
     double delta, NumericVector SS, const int K, const int Nall)
{
// Elements to undertake the posterior sampling
double c1=0, c2=0, temp=0;


// Compute the normalising constant
c1 = K * log(sum(exp(-delta * pow(SS - Gstar,2))));
     
     for(int k = 0; k < Nall; k++) 
     {     
     temp =  log(sum(exp(-alpha * pow(SS-Z[k],2) -delta * pow(SS - Gstar,2))));
     c2 = c2 + temp;
     }
     
// Return the normalising constant
NumericVector c3(2);
c3[0] = c1;
c3[1] = c2;
return c3;
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
List SPTICARphiVarb(NumericMatrix W, const int nsites, const int ntimes, 
                          NumericVector phiVarb, NumericVector nneighbours, double tau, 
                          const NumericVector y, const NumericVector E, 
                          const double phiVarb_tune, double alpha, NumericVector XB,
                          const double beta_tune){
  
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
List updatetripList(List trips, NumericVector vold, 
                    NumericVector vnew, const int nedges){    
   //   create a clone of the triplet matrix                   
   List temporary = clone(trips);    
   //   seperate into the diagonal elements and off diagonal elements
   NumericMatrix outdiagtrip = temporary[0];
   NumericMatrix outoffdiagtrip = temporary[1];
   //   stuff needed for intermediate calculations
   double oldoffdiag_binary, newoffdiag_binary;
   double newoffdiag_logit,  oldoffdiag_logit;
   int    diag_inds_1,       diag_inds_2;
   double olddiag_binary_1,  olddiag_binary_2;
   //   loop over all edges, stop only at elements where vold and vnew differ
   //   then perform and update of the corresponding triplets
   for(int i = 0; i < nedges; i++) {  
     if( !(vold[i] == vnew[i]) ){
       //       this is the old off diagonal element (-ve to make +ve)
       oldoffdiag_binary = -outoffdiagtrip(i, 2);
       //       old and new v elements
       oldoffdiag_logit  = vold[i];
       newoffdiag_logit  = vnew[i];
       //       convert new v onto [0,1] scale
       newoffdiag_binary = -(1/(1 + exp(-newoffdiag_logit)));
       //       replace triplets with new v
       outoffdiagtrip(i, 2) = newoffdiag_binary;
       outoffdiagtrip(i + nedges, 2) = newoffdiag_binary;
       //       now need to find x and y coords of these offdiags to update diags
       diag_inds_1 = outoffdiagtrip(i, 0) - 1;
       diag_inds_2 = outoffdiagtrip(i, 1) - 1;
       //       get the old binary elements
       olddiag_binary_1 = outdiagtrip(diag_inds_1, 2);
       olddiag_binary_2 = outdiagtrip(diag_inds_2, 2);
       //       calculate and replace with new ones
       outdiagtrip(diag_inds_1, 2) = (olddiag_binary_1 - oldoffdiag_binary)  - newoffdiag_binary;
       outdiagtrip(diag_inds_2, 2) = (olddiag_binary_2 - oldoffdiag_binary)  - newoffdiag_binary; 
     }
   }
   //   output to a list where the first element is the updated diagonal triplets
   //   second element is the offdiagonal triplets
   List out(2);
   out[0] = outdiagtrip;
   out[1] = outoffdiagtrip;
   return out;
}

