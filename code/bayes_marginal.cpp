#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int getMarginal(NumericVector Arow, NumericVector deltas, double be,
                double alpha, double beta, int reps) {
  int l = Arow.size(), rep, prodsum, tot, i, marginal = 0;
  double randnum, frac, pext;
  NumericVector extinct(l), rands_e(l), rands_r = Rcpp::runif(reps);
  for (rep = 0; rep < reps; rep++) {
    prodsum = 0;
    tot = 0;
    rands_e = Rcpp::runif(l);
    for (i = 0; i < l; i++) {
      extinct[i] = (rands_e[i] < deltas[i] ? 1 : 0);
      prodsum += Arow[i] * extinct[i];
      tot += Arow[i];
    }
    frac = (double)prodsum / ((double)tot);
    pext = be + (1 - be) * R::pbeta(frac, alpha, beta, true, false);
    marginal += (rands_r[rep] < pext ? 1 : 0);
  }
  return(marginal);
}

// /*** R
// getMarginal(c(1, 1, 0), rep(0.5, 3), 0.5, 1, 1, 10000)
// */
