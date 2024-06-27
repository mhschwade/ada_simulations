#include <math.h>
#include <Rcpp.h>
#include "./mygen3sisfunctions.h" // Library with cpp functions to use with gen3sis.
using namespace Rcpp;

/***********************************************************************************************
*	Author: Marcelo H. Schwade
*	Contact: marceloh.schwade@gmail.com 
*
*	Last modification on April, 9th of 2024. 
***********************************************************************************************/

/******************************************************************************/
/* Function to solve a Lotka Volterra Competition EDO system using RK2 method.*/
/* This LVC model is a neutral model, without differences between species or individuals. */

/************** Guide of the parameter list (function input). *****************/
// lotka_volterra_comp_neutral2: a vector with the final abundances resulting of the 
//                      solution of Lotka-Volterra Competition model (modified).
// Nt: vector with the initial abundances of species populations
// Theta: vector with the local traits of the species. Traits are environmental 
//        optimal for conditions 1 and 2.
// pars[0]: bmax - maximum birth rate.
// pars[1]: Kr - carrying capacity by recruitment. Kr = K/rmax.
// pars[2]: alpha - interspecific competition coefficient.
// dt: the time step used to solve the system of EDOs using RK2.
// times: the size of temporal interval between the begin and the end of the 
//        EDOs solution. It's the temporal ecological timestepof the simulation.
// [[Rcpp::export]]
NumericVector lotka_volterra_comp_neutral2(NumericVector Nt, NumericVector pars, double dt=0.01, double times=1){
  
  int Ntimes;
  int nSp=Nt.length();
  NumericVector r (nSp, (pars[0] - 1));   // Define the intrinsec growth rate as a same fix value for all species.
  NumericVector k1(nSp), k2(nSp);
  
  pars.names() = CharacterVector({"bmax", "Kr", "alpha"});
  
  Ntimes = times/dt; // Integer.
  
  /*
  for(int i=0; i<nSp; i++){
    r[i] = pars["bmax"] - 1;
  }
  */


  for(int t=0; t<Ntimes; t++){
    
    k1 = dN_LVC(Nt, r, pars["Kr"], pars["alpha"])*dt;
    k2 = dN_LVC(Nt+k1, r, pars["Kr"], pars["alpha"])*dt;
    
    Nt += (k1 + k2)/2;
    
  }
  
  return Nt;
}
