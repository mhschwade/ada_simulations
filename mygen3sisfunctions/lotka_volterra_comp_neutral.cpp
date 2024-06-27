#include <math.h>
#include <Rcpp.h>
#include "./mygen3sisfunctions.h" // Library with cpp functions to use with gen3sis.
using namespace Rcpp;

/* Last modification on November, 21th of 2023. */

/******************************************************************************/
/* Function to solve a Lotka Volterra Competition EDO system using RK2 method.*/
/* This LVC model is modified to include an intrinsic growth rate depending in 
environmental conditions (2 conditions). */

/************** Guide of the parameter list (function input). *****************/
// lotka_volterra_comp: a vector with the final abundances resulting of the 
//                      solution of Lotka-Volterra Competition model (modified).
// Nt: vector with the initial abundances of species populations
// Theta: vector with the local traits of the species. Traits are environmental 
//        optimal for conditions 1 and 2.
// pars[0]: bmax - maximum birth rate.
// pars[1]: Kr - carrying capacity by recruitment. Kr = K/rmax.
// pars[2]: alpha - interspecific competition coefficient.
// pars[3]: theta1_var2 - 2*sigma_theta1^2. sigma_theta1 is the niche breath of environmental condition 1.
// pars[4]: theta2_var2 - 2*sigma_theta2^2. sigma_theta2 is the niche breath of environmental condition 2.
// pars[5]: theta1 - environmental condition 1.
// pars[6]: theta2 - environmental condition 2.
// dt: the time step used to solve the system of EDOs using RK2.
// times: the size of temporal interval between the begin and the end of the 
//        EDOs solution. It's the temporal ecological timestepof the simulation.
// [[Rcpp::export]]
NumericVector lotka_volterra_comp_neutral(NumericVector Nt, NumericMatrix Theta, NumericVector pars, double dt=0.01, double times=1){
  
  int Ntimes;
  int nSp=Nt.length();
  //double Ntotal;
  NumericVector r(nSp);
  NumericVector k1(nSp), k2(nSp);
  
  pars.names() = CharacterVector({"bmax", "Kr", "alpha", "theta1_var2", "theta2_var2", "theta1", "theta2"});
  
  Ntimes = times/dt; // Integer.
  
  //Ntotal=0;
  for(int i=0; i<nSp; i++){
    //Ntotal += Nt[i];
    
    /*
    wTheta2 = ((Theta(i, 0) - pars["theta1"])*(Theta(i, 0) - pars["theta1"])/pars["theta1_var2"] + 
              (Theta(i, 1) - pars["theta2"])*(Theta(i, 1) - pars["theta2"])/pars["theta2_var2"]);

    r[i] = pars["bmax"]*exp(-wTheta2) - 1;
    */

    // A neutral definition of r. The intrinsic groth doesn't depend on environment.
    r[i] = pars["bmax"] - 1;


    //printf("sp%d\t %.2f\t %.2f\t", i, Theta(i,0), r[i]);
    //printf("sp: %d\t Theta1: %.2f\t wTheta2: %.2f\t r:%.2f\t", i, Theta(i,0), wTheta2, r[i]);
  }
  //printf("\n");
  for(int t=0; t<Ntimes; t++){
    
    k1 = dN_LVC(Nt, r, pars["Kr"], pars["alpha"])*dt;
    k2 = dN_LVC(Nt+k1, r, pars["Kr"], pars["alpha"])*dt;
    
    Nt += (k1 + k2)/2;
  }
  
  return Nt;
}
