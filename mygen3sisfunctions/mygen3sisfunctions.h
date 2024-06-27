#include <Rcpp.h>
using namespace Rcpp;
//#include <stdio.h>
//#include <math.h>

/* This is a library to reference a set of functions in cpp developed to gen3sis
 simulations. */


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

// Function to calculate the derivative of Lotka-Volterra Competition model.
NumericVector dN_LVC(NumericVector Nt, NumericVector r, double Kr, double alpha);

// Function to solve a system of EDOs of Lotka-Volterra Competition model using RK2.
/*NumericVector lotka_volterra_comp(NumericVector Nt, 
                                  NumericMatrix Theta,  
                                  NumericVector pars, 
                                  double dt, 
                                  double times);
*/
/******************************************************************************/
/** Function to calculate the derivative of a Lotka-Volterra Competion model. */
// [[Rcpp: export]]
NumericVector dN_LVC(NumericVector Nt, NumericVector r, double Kr, double alpha){
  int nSp=Nt.length(), i;
  double Ntotal;
  NumericVector dN(nSp);
  
  Ntotal = 0;
  for(i=0; i<nSp; i++){
    Ntotal += Nt[i];
  }
  
  for(i=0; i<nSp; i++){
    dN[i] = Nt[i]*(r[i] - (Nt[i] + alpha*(Ntotal - Nt[i]))/Kr);
  }
  
  return(dN);
}
