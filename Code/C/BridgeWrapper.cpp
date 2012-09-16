#include "BridgeWrapper.hpp"
#include <time.h>
#include <iostream>
#include <math.h>
#include <exception>

using std::cout;
using std::cerr;

#ifdef USE_R
#include "R.h"
#include "Rmath.h"
#include <R_ext/Utils.h>
#endif

// #include <algorithm>

// Test function.
// double add(double *sum, double *a, double *b)
// {
//   *sum = *a + *b;
// }

//////////////////////////////////////////////////////////////////////
		    // EXPECTATION MAXIMIZATION //
//////////////////////////////////////////////////////////////////////

int EM(Matrix & beta, MatrixFrame &y, MatrixFrame &X,
	double ratio, double alpha, double lambda_max,
	double tol, int max_iter, bool use_cg)
{
  BridgeRegression br(X, y);
  int iter;

  try{
    iter =  br.EM(beta, 1.0, ratio, alpha, lambda_max, tol, max_iter, use_cg);
  }
  catch (std::exception& e) {
    printf("Error: %s\n", e.what());
    printf("Aborting EM.\n");
  }

  return iter;
}

//////////////////////////////////////////////////////////////////////
		       // BRIDGE REGRESSION //
//////////////////////////////////////////////////////////////////////

//--------------------------------------------------------------------
// Bridge Regression when sig2 and tau are known.
double bridge_regression(MatrixFrame & beta,
		       MatrixFrame & u,
		       MatrixFrame & omega,
		       const MatrixFrame & y,
		       const MatrixFrame & X,
		       double sig2_known,
		       double tau_known,
		       double alpha,
		       uint burn)
{

  BridgeRegression br(X, y);

  // uint P = X.cols();
  // uint N = X.rows();
  uint M = beta.mats();

  // Initialize beta.
  Matrix ls;
  br.least_squares(ls);
  beta[0].copy(ls);
  u[0].fill(0.0);

  // Matrix u    (P, 1, M, 0.0);
  // Matrix omega(P, 1, M, 1.0);
  Matrix sig2 (1, 1, M, sig2_known);
  Matrix tau  (1, 1, M, tau_known);

  // We must make sure that we have seeded u, omega, beta with values
  // that won't produce negative values for b.
  // beta[0].fill(tau_known); omega[0].fill(2.0);

  RNG r;

  // Details.
  cout << "Bridge Regression..." << "\n";
  cout << "Known: sig2, tau, alpha" << "\n";
  cout << "Burn-in: " << burn << ", Num. Samples: " << M << "\n";

  // Keep track of time.
  clock_t start, end;

  try {

    start = clock();

    // I think I want to to go:
    // Set u = 0.
    // Set beta = LS.
    // Repeat.
    // sample w
    // sample beta
    // sample u.

    // Burn-In.
    for(uint i = 0; i < burn; i++){
      br.sample_omega(omega[0], beta[0], u[0], tau[0](0), alpha, r);
      br.sample_u(u[0], beta[0], omega[0], tau[0](0), alpha, r);
      br.sample_beta(beta[0], beta[0], u[0], omega[0], sig2[0](0), tau[0](0), alpha, r);
      #ifdef USE_R
      if (i % 10 == 0) R_CheckUserInterrupt();
      #endif
    }

    end = clock();

    double total_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
    printf("Expect approx. %g sec. for %i samples.\n", total_time * M / burn, M);
    
    start = clock();

    // MCMC
    for(uint i = 1; i < M; i++){
      br.sample_omega(omega[i], beta[i-1], u[i-1], tau[i-1](0), alpha, r);
      br.sample_u(u[i], beta[i-1], omega[i], tau[i](0), alpha, r);
      br.sample_beta(beta[i], beta[i-1], u[i], omega[i], sig2[i](0), tau[i](0), alpha, r);
      #ifdef USE_R
      if (i % 10 == 0) R_CheckUserInterrupt();
      #endif
    }

    end = clock();

    cout << "Sampling complete: " << (double)(end - start) / CLOCKS_PER_SEC
	 << " sec. for " << M << " iterations." << "\n";

  }
  catch (std::exception& e) {
    printf("Error: %s\n", e.what());
    printf("Aborting Gibbs sampler.\n");
  }

  return (double)(end - start) / CLOCKS_PER_SEC;

} // bridge_regression

//--------------------------------------------------------------------
double bridge_regression(MatrixFrame & beta,
		       MatrixFrame & u,
		       MatrixFrame & omega,
		       MatrixFrame & sig2,
		       const MatrixFrame & y,
		       const MatrixFrame & X,
		       double tau_known,
		       double alpha,
		       double sig2_shape,
		       double sig2_scale,
		       uint burn,
		       int beta_iter)
{
  BridgeRegression br(X, y);

  // uint P = X.cols();
  // uint N = X.rows();
  uint M = beta.mats();

  // Initialize beta.
  Matrix ls;
  br.least_squares(ls);
  beta[0].copy(ls);
  u[0].fill(0.0);

  //Matrix u    (P, 1, M, 0.0);
  //Matrix omega(P, 1, M, 1.0);
  Matrix tau  (1, 1, M, tau_known);

  // beta[0].fill(tau_known); omega[0].fill(2.0);

  RNG r;

  // Details.
  cout << "Bridge Regression..." << "\n";
  cout << "Known: tau, alpha" << "\n";
  cout << "Burn-in: " << burn << ", Num. Samples: " << M << "\n";

  // Keep track of time.
  clock_t start, end;

  try{

    start = clock();

    br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.

    // Burn-In.
    for(uint i = 0; i < burn+1; i++){
      br.sample_omega(omega[0], beta[0], u[0], tau[0](0), alpha, r);
      br.sample_u(u[0], beta[0], omega[0], tau[0](0), alpha, r);
      br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.
      br.sample_beta(beta[0], beta[0], u[0], omega[0], 
		     sig2[0](0), tau[0](0), alpha, r, beta_iter);
      #ifdef USE_R
      if (i % 10 == 0) R_CheckUserInterrupt();
      #endif
    }

    end = clock();

    double total_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
    printf("Expect approx. %g sec. for %i samples.\n", total_time * M / burn, M);

    start = clock();

    // MCMC
    for(uint i = 1; i < M; i++){
      br.sample_omega(omega[i], beta[i-1], u[i-1], tau[i-1](0), alpha, r);
      br.sample_u(u[i], beta[i-1], omega[i], tau[i](0), alpha, r);
      br.sample_sig2(sig2[i], beta[i-1], sig2_shape, sig2_scale, r);  // Sample sig2.
      br.sample_beta(beta[i], beta[i-1], u[i], omega[i], 
		     sig2[i](0), tau[i](0), alpha, r, beta_iter);
      #ifdef USE_R
      if (i % 10 == 0) R_CheckUserInterrupt();
      #endif
    }

    end = clock();

    cout << "Sampling complete: " << (double)(end - start) / CLOCKS_PER_SEC
	 << " sec. for " << M << " iterations." << "\n";

  }
  catch (std::exception& e) {
    printf("Error: %s\n", e.what());
    printf("Aborting Gibbs sampler.\n");
  }

  return (double)(end - start) / CLOCKS_PER_SEC;

} // bridge_regression

//--------------------------------------------------------------------
double bridge_regression(MatrixFrame & beta,
			 MatrixFrame & u,
			 MatrixFrame & omega,
			 MatrixFrame & sig2,
			 MatrixFrame & tau,
			 const MatrixFrame & y,
			 const MatrixFrame & X,
			 double alpha,
			 double sig2_shape,
			 double sig2_scale,
			 double nu_shape,
			 double nu_rate,
			 uint burn)
{

  BridgeRegression br(X, y);

  // uint P = X.cols();
  // uint N = X.rows();
  uint M = beta.mats();

  // Initialize beta.
  Matrix ls;
  br.least_squares(ls);
  beta[0].copy(ls);
  u[0].fill(0.0);

  // Matrix u    (P, 1, M, 0.0);
  // Matrix omega(P, 1, M, 1.0);

  RNG r;

  // Details.
  cout << "Bridge Regression..." << "\n";
  cout << "Known: alpha" << "\n";
  cout << "Burn-in: " << burn << ", Num. Samples: " << M << "\n";

  // Keep track of time.
  clock_t start, end;

  try {

    br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.
    br.sample_tau(tau[0], beta[0], alpha, nu_shape, nu_rate, r);  // Sample tau.

    start = clock();

    // THE ORDER OF SAMPLING IS EXTREMELY IMPORTANT - SEE APPENDIX.

    // Burn-In.
    for(uint i = 0; i < burn; i++){
      br.sample_tau(tau[0], beta[0], alpha, nu_shape, nu_rate, r);  // Sample tau.
      br.sample_omega(omega[0], beta[0], u[0], tau[0](0), alpha, r);
      br.sample_u(u[0], beta[0], omega[0], tau[0](0), alpha, r);
      br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.
      br.sample_beta(beta[0], beta[0], u[0], omega[0], sig2[0](0), tau[0](0), alpha, r);
      #ifdef USE_R
      if (i % 10 == 0) R_CheckUserInterrupt();
      #endif
    }

    end = clock();

    double total_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
    printf("Expect approx. %g sec. for %i samples.\n", total_time * M / burn, M);

    start = clock();

    // MCMC
    for(uint i = 1; i < M; i++){
      br.sample_tau(tau[i], beta[i-1], alpha, nu_shape, nu_rate, r);  // Sample tau.
      br.sample_omega(omega[i], beta[i-1], u[i-1], tau[i](0), alpha, r);
      br.sample_u(u[i], beta[i-1], omega[i], tau[i](0), alpha, r);
      br.sample_sig2(sig2[i], beta[i-1], sig2_shape, sig2_scale, r);  // Sample sig2.
      br.sample_beta(beta[i], beta[i-1], u[i], omega[i], sig2[i](0), tau[i](0), alpha, r);
      #ifdef USE_R
      if (i % 10 == 0) R_CheckUserInterrupt();
      #endif
    }

    end = clock();

    cout << "Sampling complete: "  << (double)(end - start) / CLOCKS_PER_SEC 
	 << " sec. for " << M << " iterations." << "\n";

  }
  catch (std::exception& e) {
    printf("Error: %s\n", e.what());
    printf("Aborting Gibbs sampler.\n");
  }

  return (double)(end - start) / CLOCKS_PER_SEC;

} // bridge_regression

//--------------------------------------------------------------------
double bridge_regression(MatrixFrame & beta,
			 MatrixFrame & u,
			 MatrixFrame & omega,
			 MatrixFrame & sig2,
			 MatrixFrame & tau,
			 const MatrixFrame & y,
			 const MatrixFrame & X,
			 double alpha,
			 double sig2_shape,
			 double sig2_scale,
			 double nu_shape,
			 double nu_rate,
			 double true_sig2, // 
			 double true_tau , // Stored in tau  if true.
			 uint burn)
{
  BridgeRegression br(X, y);

  // uint P = X.cols();
  // uint N = X.rows();
  uint M = beta.mats();
  bool know_sig2 = true_sig2 > 0;
  bool know_tau  = true_tau  > 0;

  RNG r;

  // Details.
  printf("Bridge Regression: known alpha=%g", alpha);
  if (know_sig2) printf(", sig2=%g", true_sig2);
  if (know_tau ) printf(", tau=%g", true_tau);
  printf("\nBurn-in: %i, Num. Samples: %i\n", burn, M);

  // Initialize.
  Matrix ls;
  br.least_squares(ls);
  beta[0].copy(ls);

  u[0].fill(0.0);

  if ( know_sig2 ) sig2.fill(true_sig2);
  if ( know_tau  ) tau.fill(true_tau);

  // Keep track of time.
  clock_t start, end;

  try {

    if (!know_sig2) br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.
    if (!know_tau ) br.sample_tau(tau[0], beta[0], alpha, nu_shape, nu_rate, r);  // Sample tau.

    start = clock();

    // Burn-In.
    for(uint i = 0; i < burn; i++){
      if (!know_tau) br.sample_tau(tau[0], beta[0], alpha, nu_shape, nu_rate, r);  // Sample tau.
      br.sample_omega(omega[0], beta[0], u[0], tau[0](0), alpha, r);
      br.sample_u(u[0], beta[0], omega[0], tau[0](0), alpha, r);
      if (!know_sig2) br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.
      br.sample_beta(beta[0], beta[0], u[0], omega[0], sig2[0](0), tau[0](0), alpha, r);
      #ifdef USE_R
      if (i % 10 == 0) R_CheckUserInterrupt();
      #endif
    }

    end = clock();

    double total_time = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
    printf("Expect approx. %g sec. for %i samples.\n", total_time * M / burn, M);

    start = clock();

    // MCMC
    for(uint i = 1; i < M; i++){
      if (!know_tau) br.sample_tau(tau[i], beta[i-1], alpha, nu_shape, nu_rate, r);  // Sample tau.
      br.sample_omega(omega[i], beta[i-1], u[i-1], tau[i](0), alpha, r);
      br.sample_u(u[i], beta[i-1], omega[i], tau[i](0), alpha, r);
      if (!know_sig2) br.sample_sig2(sig2[i], beta[i-1], sig2_shape, sig2_scale, r);  // Sample sig2.
      br.sample_beta(beta[i], beta[i-1], u[i], omega[i], sig2[i](0), tau[i](0), alpha, r);
      #ifdef USE_R
      if (i % 10 == 0) R_CheckUserInterrupt();
      #endif
    }

    end = clock();

    cout << "Sampling complete: "  << (double)(end - start) / CLOCKS_PER_SEC 
	 << " sec. for " << M << " iterations." << "\n";

  }
  catch (std::exception& e) {
    printf("Error: %s\n", e.what());
    printf("Aborting Gibbs sampler.\n");
  }

  return (double)(end - start) / CLOCKS_PER_SEC;

} // bridge_regression

//////////////////////////////////////////////////////////////////////
		   // BRIDGE REGRESSION WRAPPERS //
//////////////////////////////////////////////////////////////////////

void bridge_EM(double *betap,
	       const double *yp,
	       const double *Xp,
	       const double *ratio,
	       const double *alpha,
	       const int *P,
	       const int *N,
	       const double *lambda_max,
	       const double *tol,
	             int *max_iter,
	       const bool *use_cg)
{
  Matrix beta(*P, 1);

  Matrix y(yp, *N, 1);
  Matrix X(Xp, *N, *P);

  if (*use_cg) cout << "Using conjugate gradient method.\n";

  int total_iter = EM(beta, y, X, *ratio, *alpha, *lambda_max, *tol, *max_iter, *use_cg);
  *max_iter = total_iter;

  MatrixFrame beta_mf(betap, *P);
  beta_mf.copy(beta);
}

//--------------------------------------------------------------------
void bridge_reg_know_sig2(double *betap,
			  double *up,
			  double *omegap,
			  const double *yp,
			  const double *Xp,
			  const double *sig2_known,
			  const double *tau_known,
			  const double *alpha,
			  const int *P,
			  const int *N,
			  const int *M,
			  const int *burn,
			  double *runtime)
{
  Matrix y(yp, *N,  1);
  Matrix X(Xp, *N, *P);

  Matrix beta (*P, 1, *M);
  Matrix u    (*P, 1, *M, 0.0);
  Matrix omega(*P, 1, *M, 1.0);

  #ifdef USE_R
  GetRNGstate();
  #endif

  *runtime = bridge_regression(beta, u, omega, y, X, 
			       *sig2_known, *tau_known, *alpha, (uint)(*burn));

  //  bridge_regression(beta, y, X, *sig2_known, *tau_known, *alpha, (uint)(*burn));

  #ifdef USE_R
  PutRNGstate();
  #endif

  MatrixFrame beta_mf (betap , *P, 1, *M);
  MatrixFrame u_mf    (up    , *P, 1, *M);
  MatrixFrame omega_mf(omegap, *P, 1, *M);
  beta_mf.copy(beta);
  u_mf.copy(u);
  omega_mf.copy(omega);

} // bridge_regression

//--------------------------------------------------------------------
void bridge_reg_know_tau(double *betap,
			 double *up,
			 double *omegap,
			 double *sig2p,
			 const double *yp,
			 const double *Xp,
			 const double *tau_known,
			 const double *alpha,
			 const double *sig2_shape,
			 const double *sig2_scale,
			 const int *P,
			 const int *N,
			 const int *M,
			 const int *burn,
			 double *runtime,
			 const int *beta_iter)
{
  Matrix y(yp, *N,  1);
  Matrix X(Xp, *N, *P);

  Matrix beta (*P, 1, *M);
  Matrix u    (*P, 1, *M, 0.0);
  Matrix omega(*P, 1, *M, 1.0);
  Matrix sig2 ( 1, 1, *M);

  #ifdef USE_R
  GetRNGstate();
  #endif

  *runtime = bridge_regression(beta, u, omega, sig2, y, X,
			       *tau_known, *alpha,
			       *sig2_shape, *sig2_scale,
			       (uint)*burn,
			       *beta_iter);

  #ifdef USE_R
  PutRNGstate();
  #endif

  // std::copy(&beta(0), &beta(0) + beta.vol(), betap);

  MatrixFrame beta_mf (betap, *P, 1, *M);
  MatrixFrame u_mf    (up    , *P, 1, *M);
  MatrixFrame omega_mf(omegap, *P, 1, *M);
  MatrixFrame sig2_mf (sig2p,  1, 1, *M);

  beta_mf.copy(beta);
  u_mf.copy(u);
  omega_mf.copy(omega);
  sig2_mf.copy(sig2);

} // bridge_regression

//--------------------------------------------------------------------
void bridge_regression(double *betap,
		       double *up,
		       double *omegap,
		       double *sig2p,
		       double *taup,
		       const double *yp,
		       const double *Xp,
		       const double *alpha,
		       const double *sig2_shape,
		       const double *sig2_scale,
		       const double *nu_shape,
		       const double *nu_rate,
		       const int *P,
		       const int *N,
		       const int *M,
		       const int *burn,
		       double *runtime)
{
  Matrix y(yp, *N,  1);
  Matrix X(Xp, *N, *P);

  Matrix beta(*P, 1, *M);
  Matrix u    (*P, 1, *M, 0.0);
  Matrix omega(*P, 1, *M, 1.0);
  Matrix sig2( 1, 1, *M);
  Matrix tau ( 1, 1, *M);

  #ifdef USE_R
  GetRNGstate();
  #endif

  *runtime = bridge_regression(beta, u, omega, sig2, tau, y, X,
			       *alpha,
			       *sig2_shape, *sig2_scale,
			       *nu_shape, *nu_rate,
			       *burn);

  #ifdef USE_R
  PutRNGstate();
  #endif

  MatrixFrame beta_mf (betap, *P, 1 , *M);
  MatrixFrame u_mf    (up    , *P, 1, *M);
  MatrixFrame omega_mf(omegap, *P, 1, *M);
  MatrixFrame sig2_mf (sig2p, 1 , 1 , *M);
  MatrixFrame tau_mf  (taup , 1 , 1 , *M);

  beta_mf.copy(beta);
  u_mf.copy(u);
  omega_mf.copy(omega);
  sig2_mf.copy(sig2);
  tau_mf.copy (tau );

} // bridge_regression

void bridge_regression_all(double *betap,
		       double *up,
		       double *omegap,
		       double *sig2p,
		       double *taup,
		       const double *yp,
		       const double *Xp,
		       const double *alpha,
		       const double *sig2_shape,
		       const double *sig2_scale,
		       const double *nu_shape,
		       const double *nu_rate,
		       const double *true_sig2,
		       const double *true_tau,
		       const int *P,
		       const int *N,
		       const int *M,
		       const int *burn,
		       double *runtime)
{
  Matrix y(yp, *N,  1);
  Matrix X(Xp, *N, *P);

  Matrix beta(*P, 1, *M);
  Matrix u    (*P, 1, *M, 0.0);
  Matrix omega(*P, 1, *M, 1.0);
  Matrix sig2( 1, 1, *M);
  Matrix tau ( 1, 1, *M);

  #ifdef USE_R
  GetRNGstate();
  #endif

  *runtime = bridge_regression(beta, u, omega, sig2, tau, y, X,
			       *alpha,
			       *sig2_shape, *sig2_scale,
			       *nu_shape, *nu_rate,
			       *true_sig2, *true_tau,
			       *burn);

  #ifdef USE_R
  PutRNGstate();
  #endif

  MatrixFrame beta_mf (betap, *P, 1 , *M);
  MatrixFrame u_mf    (up    , *P, 1, *M);
  MatrixFrame omega_mf(omegap, *P, 1, *M);
  MatrixFrame sig2_mf (sig2p, 1 , 1 , *M);
  MatrixFrame tau_mf  (taup , 1 , 1 , *M);

  beta_mf.copy(beta);
  u_mf.copy(u);
  omega_mf.copy(omega);
  sig2_mf.copy(sig2);
  tau_mf.copy (tau );

} // bridge_regression


/*
void passnothing(void)
{
  // The data.
  Matrix X; X.read("X.lars", true);
  Matrix y; y.read("Y.lars", true);

  // True data.
  Matrix beta_data; beta_data.read("beta.data", true);

  uint M = 10000;
  uint P = X.cols();

  double alpha = 0.5;
  double tau_seed = 0.01;

  // Least squares.
  Matrix XX(X, X, 'T', 'N');
  Matrix ls(X, y, 'T', 'N');
  symsolve(XX, ls);

  cout << "LS:\n" << ls;

  Matrix beta(P      , (uint)1, M);
  Matrix sig2((uint)1, (uint)1, M);

  // bridge_regression(beta, y, X, 10.0, 1.0, 0.5, 100);
  bridge_regression(beta, sig2, y, X, 10.0, 0.5, 0.0, 0.0, 500);

  // beta.write("beta.post");

  // // EM
  // Matrix beta_EM(P);
  // int iter_direct = br.EM(beta_EM, tau_seed, 1.0, alpha, 10e8, 10e-9, 1000);
  // cout << "iters (dir): " << iter_direct << "\n";
  // cout << "beta (EM):\n" << beta_EM;

}

// To check speed issues.
void donothing(double* qbeta,
	       double* qsig2,
	       double* qtau,
	       double* qy,
	       double* qX,
	       double* qalpha,
	       double* qsig2_shape,
	       double* qsig2_scale,
	       double* qnu_shape,
	       double* qnu_rate,
	       int* qP,
	       int* qN,
	       int* qM,
	       int* qburn)
{
  // The data.
  Matrix X; X.read("X.lars", true);
  Matrix y; y.read("Y.lars", true);

  // True data.
  Matrix beta_data; beta_data.read("beta.data", true);

  uint M = 10000;
  uint P = X.cols();

  double alpha = 0.5;
  double tau_seed = 0.01;

  // Least squares.
  Matrix XX(X, X, 'T', 'N');
  Matrix ls(X, y, 'T', 'N');
  symsolve(XX, ls);

  cout << "LS:\n" << ls;

  Matrix beta(P      , (uint)1, M);
  Matrix sig2((uint)1, (uint)1, M);

  // bridge_regression(beta, y, X, 10.0, 1.0, 0.5, 100);
  bridge_regression(beta, sig2, y, X, 10.0, 0.5, 0.0, 0.0, 500);

  // beta.write("beta.post");

  // // EM
  // Matrix beta_EM(P);
  // int iter_direct = br.EM(beta_EM, tau_seed, 1.0, alpha, 10e8, 10e-9, 1000);
  // cout << "iters (dir): " << iter_direct << "\n";
  // cout << "beta (EM):\n" << beta_EM;

}
*/

//--------------------------------------------------------------------
void bridge_regression_test(MatrixFrame & beta,
			    MatrixFrame & u,
			    MatrixFrame & omega,
			    MatrixFrame & tau,
			    const MatrixFrame & y,
			    const MatrixFrame & X,
			    double sig2_known,
			    double alpha,
			    double nu_shape,
			    double nu_rate,
			    uint burn)
{
  BridgeRegression br(X, y);

  // uint P = X.cols();
  // uint N = X.rows();
  uint M = beta.mats();

  // Initialize beta.
  Matrix ls;
  br.least_squares(ls);
  beta[0].copy(ls);

  //Matrix u    (P, 1, M, 0.0);
  //Matrix omega(P, 1, M, 1.0);
  Matrix sig2  (1, 1, M, sig2_known);

  // beta[0].fill(tau_known); omega[0].fill(2.0);

  RNG r;

  // Details.
  cout << "Bridge Regression..." << "\n";
  cout << "Known: tau, alpha" << "\n";
  cout << "Burn-in: " << burn << ", Num. Samples: " << M << "\n";

  // Keep track of time.
  clock_t start, end;

  start = clock();

  // br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.
  omega[0].fill(2.0);
  tau[0](0) = 2.0;
  beta[0].fill(1.0);

  // Burn-In.
  for(uint i = 0; i < burn; i++){
    br.sample_omega(omega[0], beta[0], u[0], tau[0](0), alpha, r);
    br.sample_tau(tau[0], beta[0], alpha, nu_shape, nu_rate, r);  // Sample tau.
    br.sample_beta(beta[0], beta[0], u[0], omega[0], sig2[0](0), tau[0](0), alpha, r);
    br.sample_u(u[0], beta[0], omega[0], tau[0](0), alpha, r);
    // cout << u[0] << " " << omega[0] << " " << tau[0] << "\n";
  }

  end = clock();
  
  double total_time = (double)(end - start) / CLOCKS_PER_SEC;
  printf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
  printf("Expect approx. %g sec. for %i samples.\n", total_time * M / burn, M);

  start = clock();

  // MCMC
  for(uint i = 1; i < M; i++){
    br.sample_omega(omega[i], beta[i-1], u[i-1], tau[i-1](0), alpha, r);
    br.sample_tau(tau[i], beta[i-1], alpha, nu_shape, nu_rate, r);  // Sample tau.
    br.sample_beta(beta[i], beta[i-1], u[i-1], omega[i], sig2[i](0), tau[i](0), alpha, r);
    br.sample_u(u[i], beta[i], omega[i], tau[i](0), alpha, r);
  }

  end = clock();

  cout << "Sampling complete: " << (double)(end - start) / CLOCKS_PER_SEC
       << " sec. for " << M << " iterations." << "\n";

} // bridge_regression

void rtnorm_left(double *x, double *left, double *mu, double *sig, int *num)
{
  RNG r;

  #ifdef USE_R
  GetRNGstate();
  #endif

  for(int i=0; i < *num; ++i){
    #ifdef USE_R
    if (i%1==0) R_CheckUserInterrupt();
    #endif

    x[i] = r.tnorm(left[i], mu[i], sig[i]);
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
}

void rtnorm_both(double *x, double *left, double* right, double *mu, double *sig, int *num)
{
  RNG r;

  #ifdef USE_R
  GetRNGstate();
  #endif

  for(int i=0; i < *num; ++i){
    #ifdef USE_R
    if (i%1==0) R_CheckUserInterrupt();
    #endif

    x[i] = r.tnorm(left[i], right[i], mu[i], sig[i]);
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
}

//--------------------------------------------------------------------
double bridge_regression_stable(MatrixFrame & beta,
			      MatrixFrame & lambda,
			      MatrixFrame & sig2,
			      const MatrixFrame & y,
			      const MatrixFrame & X,
			      double tau_known,
			      double alpha,
			      double sig2_shape,
			      double sig2_scale,
			      uint burn)
{
  BridgeRegression br(X, y);

  // uint P = X.cols();
  // uint N = X.rows();
  uint M = beta.mats();

  // Initialize beta.
  Matrix ls;
  br.least_squares(ls);
  beta[0].copy(ls);

  //Matrix u    (P, 1, M, 0.0);
  //Matrix omega(P, 1, M, 1.0);
  Matrix tau  (1, 1, M, tau_known);

  // beta[0].fill(tau_known); omega[0].fill(2.0);

  RNG r;

  // Details.
  cout << "Bridge Regression..." << "\n";
  cout << "Known: tau, alpha" << "\n";
  cout << "Burn-in: " << burn << ", Num. Samples: " << M << "\n";

  // Keep track of time.
  clock_t start, end;

  start = clock();

  // br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.

  // Burn-In.
  for(uint i = 0; i < burn+1; i++){
    #ifdef USE_R
    R_CheckUserInterrupt();
    #endif
    br.sample_lambda(lambda[0], beta[0], alpha, tau[0](0), r);
    br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.
    br.sample_beta_stable(beta[0], lambda[0], alpha, sig2[0](0), tau[0](0), r);
  }

  end = clock();

  double total_time = (double)(end - start) / CLOCKS_PER_SEC;
  printf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
  printf("Expect approx. %g sec. for %i samples.\n", total_time * M / burn, M);

  start = clock();

  // MCMC
  for(uint i = 1; i < M; i++){
    #ifdef USE_R
    R_CheckUserInterrupt();
    #endif
    br.sample_lambda(lambda[i], beta[i-1], alpha, tau[i-1](0), r);
    br.sample_sig2(sig2[i], beta[i-1], sig2_shape, sig2_scale, r);  // Sample sig2.
    br.sample_beta_stable(beta[i], lambda[i], alpha, sig2[i](0), tau[i-1](0), r);
  }

  end = clock();

  cout << "Sampling complete: " << (double)(end - start) / CLOCKS_PER_SEC
       << " sec. for " << M << " iterations." << "\n";

  return (double)(end - start) / CLOCKS_PER_SEC;

} // bridge_regression

//--------------------------------------------------------------------
double bridge_regression_stable(MatrixFrame & beta,
				MatrixFrame & lambda,
				MatrixFrame & sig2,
				MatrixFrame & tau,
				const MatrixFrame & y,
				const MatrixFrame & X,
				double alpha,
				double sig2_shape,
				double sig2_scale,
				double nu_shape,
				double nu_rate,
				uint burn)
{
  BridgeRegression br(X, y);

  // uint P = X.cols();
  // uint N = X.rows();
  uint M = beta.mats();

  // Initialize beta.
  Matrix ls;
  br.least_squares(ls);
  beta[0].copy(ls);

  //Matrix u    (P, 1, M, 0.0);
  //Matrix omega(P, 1, M, 1.0);
  //Matrix tau  (1, 1, M, tau_known);

  // beta[0].fill(tau_known); omega[0].fill(2.0);

  RNG r;

  // Details.
  cout << "Bridge Regression..." << "\n";
  cout << "Known: tau, alpha" << "\n";
  cout << "Burn-in: " << burn << ", Num. Samples: " << M << "\n";

  // Keep track of time.
  clock_t start, end;

  start = clock();

  // br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.

  // Burn-In.
  for(uint i = 0; i < burn+1; i++){
    #ifdef USE_R
    R_CheckUserInterrupt();
    #endif
    br.sample_tau(tau[0], beta[0], alpha, nu_shape, nu_rate, r);  // Sample tau.
    br.sample_lambda(lambda[0], beta[0], alpha, tau[0](0), r);
    br.sample_sig2(sig2[0], beta[0], sig2_shape, sig2_scale, r);  // Sample sig2.
    br.sample_beta_stable(beta[0], lambda[0], alpha, sig2[0](0), tau[0](0), r);
  }

  end = clock();

  double total_time = (double)(end - start) / CLOCKS_PER_SEC;
  printf("Burn-in complete: %g sec. for %i iterations.\n", total_time, burn);
  printf("Expect approx. %g sec. for %i samples.\n", total_time * M / burn, M);

  start = clock();

  // MCMC
  for(uint i = 1; i < M; i++){
    #ifdef USE_R
    R_CheckUserInterrupt();
    #endif
    br.sample_tau(tau[i], beta[i-1], alpha, nu_shape, nu_rate, r);  // Sample tau.
    br.sample_lambda(lambda[i], beta[i-1], alpha, tau[i](0), r);
    br.sample_sig2(sig2[i], beta[i-1], sig2_shape, sig2_scale, r);  // Sample sig2.
    br.sample_beta_stable(beta[i], lambda[i], alpha, sig2[i](0), tau[i-1](0), r);
  }

  end = clock();

  cout << "Sampling complete: " << (double)(end - start) / CLOCKS_PER_SEC
       << " sec. for " << M << " iterations." << "\n";

  return (double)(end - start) / CLOCKS_PER_SEC;

} // bridge_regression

//--------------------------------------------------------------------
void bridge_reg_know_tau_stable(double *betap,
				double *lambdap,
				double *sig2p,
				const double *yp,
				const double *Xp,
				const double *tau_known,
				const double *alpha,
				const double *sig2_shape,
				const double *sig2_scale,
				const int *P,
				const int *N,
				const int *M,
				const int *burn,
				double *runtime)
{
  Matrix y(yp, *N,  1);
  Matrix X(Xp, *N, *P);

  Matrix beta  (*P, 1, *M);
  Matrix lambda(*P, 1, *M, 1.0);
  Matrix sig2  ( 1, 1, *M);

  #ifdef USE_R
  GetRNGstate();
  #endif

  *runtime = bridge_regression_stable(beta, lambda, sig2, y, X,
				      *tau_known, *alpha,
				      *sig2_shape, *sig2_scale,
				      (uint)*burn);

  #ifdef USE_R
  PutRNGstate();
  #endif

  // std::copy(&beta(0), &beta(0) + beta.vol(), betap);

  MatrixFrame beta_mf (betap, *P, 1, *M);
  MatrixFrame lambda_mf(lambdap, *P, 1, *M);
  MatrixFrame sig2_mf (sig2p,  1, 1, *M);

  beta_mf.copy(beta);
  lambda_mf.copy(lambda);
  sig2_mf.copy(sig2);

} // bridge_regression

//--------------------------------------------------------------------
void bridge_reg_stable(double *betap,
		       double *lambdap,
		       double *sig2p,
		       double *taup,
		       const double *yp,
		       const double *Xp,
		       const double *alpha,
		       const double *sig2_shape,
		       const double *sig2_scale,
		       const double *nu_shape,
		       const double *nu_rate,
		       const int *P,
		       const int *N,
		       const int *M,
		       const int *burn,
		       double *runtime)
{
  Matrix y(yp, *N,  1);
  Matrix X(Xp, *N, *P);

  Matrix beta  (*P, 1, *M);
  Matrix lambda(*P, 1, *M, 1.0);
  Matrix sig2  ( 1, 1, *M);
  Matrix tau   ( 1, 1, *M);

  #ifdef USE_R
  GetRNGstate();
  #endif

  *runtime = bridge_regression_stable(beta, lambda, sig2, tau,
				      y, X,
				      *alpha,
				      *sig2_shape, *sig2_scale,
				      *nu_shape, *nu_rate,
				      (uint)*burn);

  #ifdef USE_R
  PutRNGstate();
  #endif

  // std::copy(&beta(0), &beta(0) + beta.vol(), betap);

  MatrixFrame beta_mf (betap, *P, 1, *M);
  MatrixFrame lambda_mf(lambdap, *P, 1, *M);
  MatrixFrame sig2_mf (sig2p,  1, 1, *M);
  MatrixFrame tau_mf  (taup,   1, 1, *M);

  beta_mf.copy(beta);
  lambda_mf.copy(lambda);
  sig2_mf.copy(sig2);
  tau_mf.copy(tau);

} // bridge_regression


//////////////////////////////////////////////////////////////////////
// END OF CODE //
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
			    // APPENDIX //
//////////////////////////////////////////////////////////////////////

/*
  It appears that the order in which I draw u, beta, omega, and tau
  makes a difference.  I need to figure out why this is.

  When I drew thing by omega, beta, u, tau I got a different answer
  than when I drew things by u, omega, tau, beta.  James uses the
  latter.

  sig2 is independent of u and omega given beta.  So its order doesn't
  matter.

  Order clearly matters because we want to sample u immediately after
  we have sampled beta.  Otherwise we might end up with
  1-tau|beta_j|omega_j^{-1/alpha} being negative.

  I think this has something to do with how we are thinking about
  likelihood.  The posterior for tau is calculated as if we have
  integrated out the auxillery variables u and omega.  But I haven't
  sampled omega and u "jointly" if I sample tau in between.

  I think we need to sample beta, u, omega, tau.  We want u to follow
  beta and, apparently, we need u and omega to be sampled using the
  same tau.  I suppose this makes sense.  Or the same beta and tau for
  that matter.  We can't switch beta and tau here because we want to
  sample u immediately after sampling beta.  u and omega need to be
  sampled using the same tau because they were constructed using the
  same tau.

  We didn't calculate the posterior distribution of tau using (beta,
  omega, u) joint likelihood.  Presumably we could do this.  I'd need
  to write routine for some sort of truncated distribution.  I think
  it would just be \1{|beta_j| / ((1-u_j) * omega_j^{1/alpha})<=tau}
  p(tau).
 */
