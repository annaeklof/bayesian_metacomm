#include <Rcpp.h>
using namespace Rcpp;

// calculate local extinction probs (deltas) of a focal species, in each patch
// Input:
// - Arow: one row of the food web's adjacency matrix, corresponding to the
//         focal species
// - pmat: matrix of patch occupancies; pmat(i,k) is species i's occupancy of
//         patch k
// - dmat: matrix of extinction probabilities; dmat(i,k) is species i's extinction
//         probability in patch k
// - pi_sp: baseline extinction probabilities of the focal species in each patch
// - alpha: first parameter of the Beta distribution
// - beta: second parameter of the Beta distribution
// - nreps: number of iterations for the Monte Carlo sim of the Bayesian network
// Output:
// - vector of extinction probabilities for focal species in each patch
// [[Rcpp::export]]
NumericVector getdelta(NumericVector Arow, NumericMatrix pmat,
                          NumericMatrix dmat, NumericVector pi_sp,
                          double alpha, double beta, int reps) {
  int S=Arow.size(); // number of species
  int N=pmat.ncol(); // number of patches
  int rep, prodsum, tot, i, k;
  double frac, pext;
  NumericVector extinct(S), rands_e(S), rands_r(reps), result(N);
  if (sum(Arow)==0) { // if species is basal:
    result=pi_sp; // extiction probs are the baseline values
  } else { // if species is nonbasal:
    for (k=0; k<N; k++) { // for each patch:
      result[k]=0.0; // zero out local extinction prob
      rands_r=Rcpp::runif(reps); // generate (reps) random numbers for future use
      for (rep=0; rep<reps; rep++) { // for (reps) repetitions:
        prodsum=0; // cumulative sum of product of Arow[i] and extinct[i]
        tot=0; // sum of Arow's entries (row sum of A[sp,])
        rands_e=Rcpp::runif(S); // generate S random numbers, for future use
        for (i=0; i<S; i++) { // for each species:
          // extinct if random value is less than 1-(1-delta)*p in local patch
          extinct[i]=(rands_e[i]<(1-(1-dmat(i,k))*pmat(i,k)) ? 1 : 0);
          prodsum+=Arow[i]*extinct[i]; // is ith species extinct?
          tot+=Arow[i]; // accumulate sum of entries of Arow
        } // end of for i
        frac=(double)prodsum/((double)tot); // fraction of extinct species
        // extinction probability of focal species (Bayesian NW formula):
        pext=pi_sp[k]+(1-pi_sp[k])*R::pbeta(frac, alpha, beta, true, false);
        // extinct if random number is less than pext
        result[k]+=(rands_r[rep]<pext ? 1.0 : 0.0);
      } // end of for rep
      result[k]/=(double)reps; // divide result by number of repetitions
    } // end of for k
  } // end of else
  return(result);
}

// obtain equilibrium patch occupancies numerically
// Input:
// - p: vector of initial patch occupancies of focal species
// - M: dispersal matrix of focal species
// - d: vector of extinction probabilities of focal species
// Output:
// - vector of updated patch occupancies for focal species in each patch
// [[Rcpp::export]]
NumericVector itersol(NumericVector p, NumericMatrix M, NumericVector d) {
  int N=p.size(); // number of patches
  int converge=0; // variable to check whether numerical iteration has converged
  int iter, k, l;
  NumericVector p_next(N), result=p, E_inverse(N), g(N);
  for (k=0; k<N; k++) E_inverse[k]=-1/log(1-d[k]); // the inverse extinction rates
  while (!converge) { // if solution hasn't converged yet:
    p_next=result; // initialize next occupancy values to the current ones
    for (iter=0; iter<10; iter++) { // for ten iterations:
      for (k=0; k<N; k++) { // calculate principal map g: colonization/extinction
        g[k]=0; // initialize ith value of principal map
        for (l=0; l<N; l++) g[k]=g[k]+M(k,l)*p_next[l]; // colonization (M %*% p)
      } // end of for k
      g=g*E_inverse; // colonization*(inverse extinction)
      p_next=g/(g+1); // iteration map: g/(g+1)
    } // end of for iter
    if (sum(pow(result-p_next, 2))<1.0e-8) converge=1; // has result converged?
    result=p_next; // in either case, update result with current estimate of p
  } // end of while
  return(result);
}
