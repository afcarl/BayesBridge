//////////////////////////////////////////////////////////////////////
/*
 Herein we implement Gibbs sampling for Bridge Regression a la Polson
 and Scott.  For a detailed description of the specifics of the setup
 and algorithm please see their paper The Bayesian Bridge
 (http://arxiv.org/abs/1109.2279).

 One starts with the basic regression:

   y = X beta + ep, ep ~ N(0, sig2 I).

 To regularize the selection of beta one may chose a variety of
 priors.  In the Bridge regression the prior lives within the familiy
 of exponential prior distributions described by

   p(beta | alpha, tau) \propto exp( - | beta_j / tau |^alpha ).

 GIBBS SAMPLING:

 The challenge is to compute the posterior distribution efficiently.
 Polson and Scott use a Normal mixure of Bartlett-Fejer kernels to
 carry out this procedure.

 See BridgeRegression.h for the conditional posteriors used when Gibbs
 Sampling.

 R WRAPPER:

 We also provide functions that may be called from R so that one can
 speed up their Gibbs sampling or prevent copying large matrices when
 doing an expectation maximization.

 In what follows:

   - y is N x 1.
   - X is N x P.
   - beta is P x 1.

 */
//////////////////////////////////////////////////////////////////////

#ifndef __BRIDGE__
#define __BRIDGE__

#include "Matrix.h"
#include "RNG.h"
#include "BridgeRegression.h"

//////////////////////////////////////////////////////////////////////
		    // EXPECTATION MAXIMIZATION //
//////////////////////////////////////////////////////////////////////

int EM(Matrix & beta, MatrixFrame &y, MatrixFrame &X,
	double ratio, double alpha, double lambda_max,
	double tol, int max_iter, bool use_cg=false);

//////////////////////////////////////////////////////////////////////
		       // BRIDGE REGRESSION //
//////////////////////////////////////////////////////////////////////

void bridge_regression(MatrixFrame & beta,
		       MatrixFrame & u,
		       MatrixFrame & omega,
		       const MatrixFrame & y,
		       const MatrixFrame & X,
		       double sig2_known,
		       double tau_known,
		       double alpha,
		       uint burn_in);

void bridge_regression(MatrixFrame & beta,
		       MatrixFrame & u,
		       MatrixFrame & omega,
		       MatrixFrame & sig2,
		       const MatrixFrame & y,
		       const MatrixFrame & X,
		       double tau,
		       double alpha,
		       double sig2_shape,
		       double sig2_scale,
		       uint burn_in);

void bridge_regression(MatrixFrame & beta,
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
		       uint burn_in);

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
			    uint burn_in);

void bridge_regression_stable(MatrixFrame & beta,
			      MatrixFrame & lambda,
			      MatrixFrame & sig2,
			      const MatrixFrame & y,
			      const MatrixFrame & X,
			      double tau_known,
			      double alpha,
			      double sig2_shape,
			      double sig2_scale,
			      uint burn_in);

//////////////////////////////////////////////////////////////////////
			    // WRAPPERS //
//////////////////////////////////////////////////////////////////////


extern "C"
{
  void bridge_EM(double *beta,
		 const double *y,
		 const double *X,
		 const double *ratio,
		 const double *alpha,
		 const int *P,
		 const int *N,
		 const double *lambda_max,
		 const double *tol,
		       int *max_iter,
		 const bool *use_cg);

  void bridge_reg_know_sig2(double *beta,
			    double *up,
			    double *omegap,
			    const double *y,
			    const double *X,
			    const double *sig2_known,
			    const double *tau_known,
			    const double *alpha,
			    const int *P,
			    const int *N,
			    const int *M,
			    const int *burn_in);

  void bridge_reg_know_tau(double *beta,
			   double *up,
			   double *omegap,
  			   double *sig2,
  			   const double *y,
  			   const double *X,
  			   const double *tau_known,
  			   const double *alpha,
  			   const double *sig2_shape,
  			   const double *sig2_scale,
  			   const int *P,
  			   const int *N,
  			   const int *M,
  			   const int *burn_in);

  void bridge_regression(double *beta,
			 double *up,
			 double *omegap,
  			 double *sig2,
  			 double *tau,
  			 const double *y,
  			 const double *X,
  			 const double *alpha,
  			 const double *sig2_shape,
  			 const double *sig2_scale,
  			 const double *nu_shape,
  			 const double *nu_rate,
  			 const int *P,
  			 const int *N,
  			 const int *M,
  			 const int *burn_in);

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
				const int *burn_in);

}

#endif

//////////////////////////////////////////////////////////////////////
			  // END OF CODE //
//////////////////////////////////////////////////////////////////////

