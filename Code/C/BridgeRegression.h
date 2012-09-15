// -*- c++-mode; -*-
//////////////////////////////////////////////////////////////////////
/*
 Herein we provide conditional posterior sampling and expectation
 maximization for Bridge Regression a la Polson and Scott.  For a
 detailed description of the specifics of the setup and algorithm
 please see their paper The Bayesian Bridge
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

 EXPECTATION MAXIMIZATION:

 They also present an expectation maximization algorithm to calculate
 the posterior mode of beta.

 */
//////////////////////////////////////////////////////////////////////

#ifndef __BRIDGEREGRESSION__
#define __BRIDGEREGRESSION__

#include "Matrix.h"
#include "RNG.hpp"
#include <iostream>
#include "retstable.c"
#include <math.h>
#include <Eigen/Core>
#include <Eigen/SVD>

using std::cout;
using std::cerr;

#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )

//////////////////////////////////////////////////////////////////////
			// CLASS DEFINITION //
//////////////////////////////////////////////////////////////////////

class BridgeRegression
{

 protected:

  // Dimension of beta.
  uint P;

  // Consider two cases: N > P and N < P.

  // Stored values, which are reused.
  Matrix y;
  Matrix X;
  Matrix XX;
  Matrix Xy;
  Matrix XX_sub;

  Matrix RT;
  Matrix RTInv;

  // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M;
  Eigen::MatrixXd tV;
  Eigen::MatrixXd a;
  Eigen::MatrixXd d;

 public:

  // Constructors:
  BridgeRegression();
  BridgeRegression(const MF& X_, const MF& y_);

  // Least squares solution.
  void least_squares(Matrix & ls);

  // For sampling beta
  void rtnorm_gibbs(MF beta, MF bmean, MF Prec, double sig2, MF b, RNG& r);
  void rtnorm_gibbs_wrapper(MF beta, double sig2, MF b, RNG& r);

  void rtnorm_gibbs(double *betap, 
		    double *ap, double *tVp, double *dp, 
		    double* bp, double* sig2p, 
		    int *Pp, RNG& r);

  void sample_beta(MF beta, const MF& beta_prev, 
		   const MF& u, const MF& omega, 
		   double sig2, double tau, double alpha, 
		   RNG& r, int niter=1);

  // For sampling everything else.
  void sample_u(MF u, const MF& beta, const MF& omega, double tau, double alpha, RNG& r);
  void sample_omega(MF omega, const MF& beta, const MF& u, double tau, double alpha, RNG& r);
  void sample_sig2(MF sig2, const MF& beta, double sig2_shape, double sig2_scale, RNG& r);
  void sample_tau(MF tau, const MF& beta, double alpha, double nu_shape, double nu_rate, RNG& r);

  void sample_lambda(MF lambda, MF beta, double alpha, double tau, RNG& r);
  void sample_beta_stable(MF beta, MF lambda, double alpha, double sig2, double tau, RNG& r);

  // Expectation Maximization.
  int EM(Matrix& beta, double sig, double tau, double alpha,
	 double lambda_max, double tol, int max_iter, bool use_cg=false);

};

// BR is BridgeRegression
#ifndef BR
typedef BridgeRegression BR;
#endif

//////////////////////////////////////////////////////////////////////
			  // Constructors //
//////////////////////////////////////////////////////////////////////

BR::BridgeRegression()
{
  std::cerr << "Warning: Default constructor called." << std::endl;
} // BridgeRegression

//--------------------------------------------------------------------
BR::BridgeRegression(const MF& X_, const MF& y_) 
  : P(X_.cols())
  , y(y_)
  , X(X_)
  , XX(P, P)
  , Xy(P)
{
  // Check conformity.
  if (X.rows()!=y.rows())
    cerr << "Error: X and y do not conform." << std::endl;
  //std:cerr << "Error: X and y do not conform." << std::endl;

  gemm(XX, X, X, 'T', 'N');
  gemm(Xy, X, y, 'T', 'N');

  // Need to deal with P = 1.
  if (P > 1) {
    XX_sub.resize(1, P-1, P);

    Matrix ss("N", P-1);
    for(uint j = 0; j < P; j++){
      XX_sub[j].copy(XX, j, ss);
      if (j < P-1) ss(j) = j;
    }
  }

  symsqrt(RT, XX);
  syminvsqrt(RTInv, XX);

  int n = X_.rows();
  int p = X_.cols();
  Eigen::Map<Eigen::MatrixXd> Xmap(&X(0), n, p);
  Eigen::Map<Eigen::MatrixXd> ymap(&y(0), n, 1);
  // Eigen::JacobiSVD<Eigen::MatrixXd> svd(Xmap, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // We need to deal with the underdetermined case.
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Xmap, Eigen::ComputeThinU | Eigen::ComputeFullV);
  Eigen::MatrixXd A = svd.matrixU() * svd.singularValues().asDiagonal();

  tV = svd.matrixV().transpose();
  a.resize(p, 1);
  d.resize(p, 1);
  int ddim = (p <= n) ? p : n;
  a.block(0,0,ddim,1) = A.transpose() * ymap;
  d.block(0,0,ddim,1) = svd.singularValues();

  // cout << d;
  // cout << "d:" << d.rows() << " " << d.cols() << "\n";

} // Bridge Regression

//////////////////////////////////////////////////////////////////////
			// Helper Functions //
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
			 // Least Squares //
//////////////////////////////////////////////////////////////////////

void BR::least_squares(Matrix & ls)
{
  try {
    ls.clone(Xy);
    symsolve(XX, ls);
  }
  catch (std::exception& e){
    ls.clone(Matrix(P)); 
    ls.fill(0.0);
    printf("Warning: cannot calculate least squares estimate; X'X is singular.\n");
    printf("Warning: setting least squares estimate to 0.0.\n");
  }
}

//////////////////////////////////////////////////////////////////////
		       // Posterior Samplers //
//////////////////////////////////////////////////////////////////////

void BR::sample_u(MF u, const MF& beta, const MF& omega, double tau, double alpha, RNG& r)
{
  double right;
  for(uint j = 0; j < P; j++){
    right = 1 - fabs(beta(j)) / tau * exp( -1.0 * log(omega(j)) / alpha );
    // u(j) = right > 0 ? r.flat(0, right) : 0;
    u(j) = r.flat(0, right);
    // COMMENT COMMENT
    if (u(j) < 0) {
      cout << "Warning: sampled negative value for u." << "\n";
      cout << beta(j) << " " << omega(j) << " " << tau << " " << right << "\n";
    }
  }
} // sample_u

//--------------------------------------------------------------------
void BR::sample_omega(MF omega, const MF& beta, const MF& u, double tau, double alpha, RNG& r)
{
  double prob;
  for(uint j = 0; j < P; j++){
    double a_j = exp( alpha * log( fabs(beta(j)) / ( (1 - u(j)) * tau) ) );
    // prob = ( 1 - alpha * (1 + a(j)) ) / (1 - alpha * a(j));
    prob = alpha / (1 + alpha * a_j);
    if (r.unif() > prob){
      omega(j) = r.gamma_rate(1.0, 1.0);
    }
    else{
      omega(j) = r.gamma_rate(2.0, 1.0);
    }
    omega(j) += a_j;
  }
} // sample_omega

//--------------------------------------------------------------------

void BR::rtnorm_gibbs(MF beta, MF bmean, MF Prec, double sig2, MF b, RNG& r)
{
  // Matrix RT, RTInv;
  // symsqrt(RT, Prec);
  // syminvsqrt(RTInv, Prec);

  Matrix RTInvZSub(P, P);
  RTInvZSub.fill(0.0);

  Matrix m(RT, bmean);
  Matrix z(RT, beta);
  double v = sig2; 

  // cout << "beta: " << beta;
  // cout << "b: " << b;
  // cout << "Prec:\n" << Prec;
  // cout << "RT:\n" << RT;
  // cout << "RTInv:\n" << RT;

  Matrix ss("N", P-1);

  for (uint i=0; i<P; i++) {

    Matrix left(P); left.fill(0.0);
    Matrix right(P); right.fill(0.0);

    for (uint j=0; j<P; j++) {

      for (uint k=0; k<P-1; k++)
	RTInvZSub(j,i) += RTInv(j,ss(k)) * z(ss(k));

      if (RTInv(j,i) > 0) {
	left(j)  = (-1.0*b(j) - RTInvZSub(j,i)) / fabs(RTInv(j,i));
	right(j) = (     b(j) - RTInvZSub(j,i)) / fabs(RTInv(j,i));
      }
      else {
	left(j)  = -1.0*(b(j) - RTInvZSub(j,i)) / fabs(RTInv(j,i));
	right(j) = (     b(j) + RTInvZSub(j,i)) / fabs(RTInv(j,i));
      }

    }

    // cout << "left:" << left;
    // cout << "right:" << right;
    // cout << i << ": z: " << z;
      
    double lmax = maxAll(left);
    double rmin = minAll(right);
    
    try {
      z(i) = r.tnorm(lmax, rmin, m(i), sqrt(v));
    }
    catch (std::exception& e) {
      cout << "left:" << left;
      cout << "right:" << right;
      cout << i << ": z: " << z;
      cout << "beta: " << beta;
      cout << "b: " << b;
      throw e;
    }

    if (i < P-1) ss(i) = i;

    // cout << "lmax: " << lmax << ", rmin: " << rmin 
    // 	 << ", z_i: " << z(i) << ", m_i: " << m(i) << "\n";
  
  }

  gemm(beta, RTInv, z);

}

////////////////////////////////////////////////////////////////////////////////

void BR::rtnorm_gibbs_wrapper(MF beta, double sig2, MF b, RNG& r)
{
  int Pint = P;
  rtnorm_gibbs(&beta(0), &a(0), &tV(0), &d(0), &b(0), &sig2, &Pint, r);
}

void BR::rtnorm_gibbs(double *betap, 
		      double *ap, double *tVp, double *dp, 
		      double* bp, double* sig2p, 
		      int *Pp, RNG& r)
{
  // Anything with a "p" suffix is a pointer.

  int P = *Pp;
  double sig = sqrt(*sig2p);
  MatrixFrame beta(betap, P);
  MatrixFrame tV(tVp, P, P);
  Matrix z(tV, beta);
  double *zp = &z(0);
  // Matrix vj(P);

  for (int i=0; i<P; i++) {
    double lmax = -INFINITY;
    double rmin =  INFINITY;

    for (int j=0; j<P; j++) {
      double vji = tVp[i+j*P];
      MatrixFrame vj(&tVp[j*P], P);
      // double rji = dot(tV.col(j), z) - vji * zp[i];
      double rji = dot(vj, z) - vji * zp[i];
      double dif = bp[j] - rji;
      double sum = bp[j] + rji;
      double left  = (vji > 0 ? -sum : -dif) / fabs(vji);
      double right = (vji > 0 ?  dif :  sum) / fabs(vji);
      lmax = MAX(lmax, left );
      rmin = MIN(rmin, right);
    }

    // double dx = rmin - lmax;
    double mean = ap[i] / (dp[i] * dp[i]);
    double sd   = sig / dp[i];

    // I need to be careful here.  It may be the case that dp is almost zero or negative!
    if (dp[i] > 1e-8){
      zp[i] = r.tnorm(lmax, rmin, mean, sd);
    }
    else {
      // double lw = lmax < rmin ? lmax : rmin;
      // double up = lmax > rmin ? lmax : rmin;
      // if (lw!=lmax) printf("Problem with lmax,rmin: %g, %g \n", lmax, rmin);
      // zp[i] = lw + (up-lw) * r.unif();
      zp[i] = r.flat(lmax, rmin);
    }

  }

  gemm(beta, tV, z, 'T', 'N');
}

////////////////////////////////////////////////////////////////////////////////

// There are multiple ways one may calculate the conditional distributions of
// beta_j | beta_{-j}.  Initially, I considered the joint distribution and then
// used regression theory to calculate the conditional.  This is a bad idea--you
// could have singular precisions.  It is better to calculate beta_j based upon
// likelihood.

#ifdef NOROTATE

void BR::sample_beta(MF beta, const MF& beta_prev, const MF& u, const MF& omega, double sig2, double tau, double alpha, RNG& r, int niter)
{
  Matrix beta_sub(P-1);
  Matrix XXbeta(1); 

  beta.copy(beta_prev);

  for(int i=0; i<niter; i++){

    Matrix ss("N", P-1);

    for(uint j = 0; j < P; j++){

     XXbeta(0) = 0.0;
      
      // If P > 1.
      beta_sub.copy(beta, ss, 0);
      gemm(XXbeta, XX_sub[j], beta_sub);
      // Else keep XXbeta(0) = 0.

      // mean and variance
      double m = ( Xy(j) - XXbeta(0) ) / XX(j,j);
      double v = sig2 / XX(j,j);

      // Calculate b(j).
      double b_j = (1.0 - u(j)) * exp( log(omega(j)) / alpha) * tau;

      beta(j) = r.tnorm(-1.0*b_j, b_j, m, sqrt(v));

      // COMMENT COMMENT
      if (fabs(beta(j)) > b_j) {
	cout << "b(j) problem: ";
	cout << b_j << " " << beta(j) << "\n";
      }
      if (j < P-1) ss(j) = j;
    }

  }

} // sample_beta

#else //------------------------------------------------------------------------------

void BR::sample_beta(MF beta, const MF& beta_prev, const MF& u, const MF& omega, double sig2, double tau, double alpha, RNG& r, int niter)
{
  Matrix b(P);
  for(uint j=0; j<P; j++)
    b(j) = (1.0 - u(j)) * exp( log(omega(j)) / alpha) * tau;

  beta.copy(beta_prev);

  // Could precompute.
  // Matrix bhat(P);
  // least_squares(bhat);

  for(int i=0; i<niter; i++) {
    // rtnorm_gibbs(beta, bhat, XX, sig2, b, r);
    rtnorm_gibbs_wrapper(beta, sig2, b, r);
  }
} // sample_beta

#endif

//--------------------------------------------------------------------
void BR::sample_sig2(MF sig2, const MF& beta, double sig2_shape, double sig2_scale, RNG& r)
{
  double shape = sig2_shape + 0.5 * y.rows();

  Matrix temp(y);
  gemm(temp, X, beta, 'N', 'N', -1.0, 1.0);

  double rss = dot(temp, temp);

  double scale = sig2_scale + 0.5 * rss;

  // cout << shape << " " << scale << "\n";

  sig2(0) = r.igamma(shape, scale);
} // sample_sig2

//--------------------------------------------------------------------
void BR::sample_tau(MF tau, const MF& beta, double alpha, double nu_shape, double nu_rate, RNG& r)
{
  double shape = nu_shape + ((double)P) / alpha;

  double rate = nu_rate;
  for(uint j = 0; j < P; j++){
    rate += exp( alpha * log(fabs(beta(j))) );
  }

  double nu = r.gamma_rate(shape, rate);

  tau(0) = exp(-1.0 * log(nu) / alpha);
} // sample_tau

//------------------------------------------------------------------------------
void BR::sample_lambda(MF lambda, MF beta, double alpha, double tau, RNG& r)
{
  for (int j=0; j<(int)P; j++)
    lambda(j) = 2 * retstable_LD(beta(j)*beta(j) / (tau*tau), 0.5 * alpha, r);
}

//------------------------------------------------------------------------------
void BR::sample_beta_stable(MF beta, MF lambda, double alpha, double sig2, double tau, RNG& r)
{
  Matrix VInv(XX);
  for(uint i=0; i<P; i++)
    VInv(i,i) += lambda(i) * sig2 / (tau * tau);
  // cout << "VInv:\n" << VInv;

  Matrix V;
  syminv(VInv, V);
  // cout << "V:\n" << V;

  Matrix L;
  chol(L, V);
  hprodeq(L, sqrt(sig2));
  // cout << "L:\n" << L;

  // The mean
  gemm(beta, V, Xy, 'N', 'N');
  // cout << "Beta:\n" << beta;

  Matrix ndraw(P);
  r.norm(ndraw, 0.0, 1.0);
  // cout << "ndraw:\n" << ndraw;

  gemm(beta, L, ndraw, 'N', 'N', 1.0, 1.0);
}

//////////////////////////////////////////////////////////////////////
			   // EM DIRECT //
//////////////////////////////////////////////////////////////////////

// Solves system directly or using conjugate gradient method.

// Returns the order of the number of "solves" needed to find a solution.  To
// find x in Ax = b you need p solves where p is the dimension of b.  Using the
// conjugate gradient method you cand find an okay x in fewer than p iterations.
// This algorithm boils down to solving a linear system several times.

int BR::EM(Matrix& beta, double sig, double tau, double alpha,
	   double lambda_max, double tol, int max_iter, bool use_cg)
{

  Matrix lambda(P);
  Matrix ss("W", P);

  Matrix XX(X, X, 'T', 'N'); // Already exists
  Matrix b(X, y, 'T', 'N');  // Already exists

  double dist = tol + 1.0;
  int    iter = 0;
  int    p    = P;
  int    N    = X.rows();

  int total_iter = p;  // We do a symsolve below.

  // Regarding comment below.  I think James let tau^* = tau/sigma and then
  // dropped the *.  Having both parameters tau and sig2 is redundant since only
  // the ratio tau/sigma matters.  However, it is somewhat confusing to read the
  // paper since they don't say this.
  double c1   = alpha * exp( (2-alpha) * (log(tau) - log(sig)) );
  double c2   = exp( -2 * (log(tau) - log(sig)) );

  Matrix EM_X(X);
  Matrix EM_b(b);
  Matrix old_beta(b);

  // Do one maximization step.
  Matrix A(XX);
  Matrix new_beta(b);
  symsolve(A, new_beta);

  while (dist > tol && iter < max_iter){

    // Expectation Step.
    uint num = 0;
    for(uint j = 0; j < (uint)p; j++){
      lambda(j) = c1 * exp( (alpha - 2) * log( fabs(new_beta(j)) ) );
      if (lambda(j) < lambda_max) {
	  ss(num) = ss(j);
	  lambda(num) = lambda(j);
	  old_beta(num) = new_beta(j);
	  ++num;
      }
    }

    // Delete entries that are too large.
    if (num < (uint)p) {
      if (num == 0){
	beta.clone(Matrix(P));
	return iter;
      }
      p = num;
      ss.resize(p);
      lambda.resize(p);
      old_beta.resize(p);
      EM_X.clone(X, Matrix("W", N), ss);
      XX.resize(p, p);
      gemm(XX, EM_X, EM_X, 'T', 'N');
      EM_b.clone(b, ss, 0);
    }

    // Maximization step.
    A.clone(XX);
    for(uint j = 0; j < (uint)p; j++)
      A(j,j) += c2 * lambda(j);
    if (!use_cg) { // Solve the system using LAPACK - by Cholesky I think.
      new_beta.clone(EM_b);
      symsolve(A, new_beta);
      total_iter = total_iter + p;
    }
    else{ // Solve the system using conjugate gradient method.
      new_beta.clone(old_beta);
      int cg_iter = cg(new_beta, A, EM_b, tol, p);
      total_iter = total_iter + cg_iter;
    }

    // Regarding the conjugate gradient method.  There is the following
    // potential problem: I found that things weren't converging when the
    // tolerance wasn't strict enough, i.e. if the conjugate gradient method
    // stops without getting close to the solution to Ax = b.  Suppose when I
    // don't set the tolerance low enough I calculate x1.  On the expectation
    // step x1 may not move lambda that much.  Then on the next maximization
    // step we are going to be solving essentially the same problem since the
    // operator A won't have changed that much.  But the tolerance will be the
    // same so we will be calculating x2 which is almost the same as x1, both of
    // which are the wrong solution.

    // UPDATE: This may have had something to do with the tolerance in the cg
    // algorithm.  Previously the tolernace wasn't an absolute thing but
    // relative to delta_0.

    // One possible solution.  At the last step, run the conjugate gradient
    // method p steps to get the best possible solution.

    // Calculate distance and increment iter.
    Matrix diff = new_beta - old_beta;
    dist = sqrt( dot(diff, diff) );
    ++iter;
  }

  beta.clone(Matrix(P));
  for(uint j = 0; j < (uint)p; j++){
    beta(ss(j)) = new_beta(j);
  }

  return total_iter;
}

//////////////////////////////////////////////////////////////////////
			  // END OF FILE //
//////////////////////////////////////////////////////////////////////

#endif

// void BR::preprocess(const MF& var)
// {
//   Matrix idx;

//   Matrix ss("N", P-1);

//   Matrix Sigma_12((uint)1, P-1);
//   Matrix Sigma_22(    P-1, P-1);

//   for(uint j = 0; j < P; j++){

//     Sigma_12.copy(var,  j, ss);
//     Sigma_22.copy(var, ss, ss);

//     // Regression Matrix.
//     Matrix Prec_22("I", P-1);
//     posv(Sigma_22, Prec_22, 'L');
//     gemm(regmat[j], Sigma_12, Prec_22);

//     // Mean Term.
//     Matrix m2(P-1);
//     m2.copy(betahat, ss, 0);
//     mean_term[j](0) = betahat(j);
//     gemm(mean_term[j], regmat[j], m2, 'N', 'N', -1.0, 1.0);

//     // Conditional Variance / SD.
//     Matrix temp(1);
//     gemm(temp, regmat[j], Sigma_12, 'N', 'T');
//     condsd(j) = sqrt(var(j,j) - temp(0));

//     if (j < P-1) ss(j) = j;
//   }

//   // cout << "mean_term:\n" << MatrixFrame(&mean_term(0), P);
//   // cout << "betahat:\n" << betahat;

// }
