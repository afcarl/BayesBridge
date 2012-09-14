
#include "Matrix/Matrix.h"
#include "RNG/RNG.h"
#include "BridgeRegression.h"

int main(int argc, char** argv)
{
  // The data.
  Matrix X; X.read("X.lars", true);
  Matrix y; y.read("Y.lars", true);

  // True data.
  Matrix beta_data; beta_data.read("beta.data", true);

  BridgeRegression br(X, y);

  uint M = 20000;
  uint P = X.cols();

  double alpha = 0.5;
  double tau_seed = 0.001;

  // Least squares.
  Matrix XX(X, X, 'T', 'N');
  Matrix ls(X, y, 'T', 'N');
  symsolve(XX, ls);

  cout << "LS:\n" << ls;

  // EM

  Matrix beta_EM(P);

  int iter_direct = br.EM(beta_EM, tau_seed, 1.0, alpha, 10e8, 10e-9, 1000);
  cout << "iters (dir): " << iter_direct << "\n";
  cout << "beta (EM):\n" << beta_EM;

  // int iter_cg     = br.EM(beta_EM, 1.0, alpha, 10e9, 10e-9, 10, true);
  // cout << "iters (cg): " << iter_cg << "\n";
  // cout << "beta:\n" << beta_EM;

  // MCMC

  // Prior parameters.
  double sig2_shape = 0.0;
  double sig2_scale = 0.0;
  double nu_shape   = 2.0;
  double nu_rate    = 2.0;

  Matrix beta (P, 1, M, 0.5); beta[0].copy(ls);
  Matrix u    (P, 1, M, 0.0);
  Matrix omega(P, 1, M, 1.0);
  Matrix sig2 (1, 1, M, 0.5);
  Matrix tau  (1, 1, M, tau_seed);

  Matrix a(P);
  Matrix b(P); b.fill(10e9);

  RNG r;

  // MCMC
  for(uint i = 1; i < M; i++){

    br.sample_u(u[i], beta[i-1], omega[i-1], tau[i-1](0), alpha, r); // Sample u.

    // Set a,b.
    for(uint j = 0; j < P; j++){
      b(j) = (1.0 - u[i](j)) * exp( log(omega[i-1](j)) / alpha) * tau[i-1](0);
      a(j) = exp( alpha * log( fabs(beta[i-1](j)) / ( (1 - u[i](j)) * tau[i-1](0)) ) );
    }

    br.sample_omega(omega[i], a, alpha, r);                       // Sample omega.
    br.sample_beta(beta[i], beta[i-1], b, sig2[i-1](0), r, 1); // Sample beta.
    br.sample_sig2(sig2[i], beta[i], sig2_shape, sig2_scale, r);  // Sample sig2.
    //br.sample_tau(tau[i], beta[i], alpha, nu_shape, nu_rate, r);  // Sample tau.

   }

  // Write to file.
  u.write("u.post", true);
  omega.write("omega.post", true);
  beta.write("beta.post", true);
  sig2.write("sig2.post", true);
  //tau.write("tau.post", true);

  return 0;
}

//////////////////////////////////////////////////////////////////////
			    // APPENDIX //
//////////////////////////////////////////////////////////////////////

/*

  Running the MCMC you see that the results are _very_ sensitive to
  the values of sig2 and tau.  In particular, if one assumes that sig2
  and tau are known they they get posteriors with orders of magnitude
  less variance.

  Also, one must be very careful of the autocorrelation induced by
  sampling beta un-jointly.  When I ran the simulation for just beta
  and sigma, which reduces to the classic case once you set b_j to be
  very large, I found that I need to run a huge number of samples to
  get my results to agree with what I was supposed to get.  The
  numbers that were off most were those that had the highest
  autocorrelation in their time series.

 */
