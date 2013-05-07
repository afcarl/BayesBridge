#include "Matrix.h"
#include "RNG.hpp"
#include "BridgeWrapper.h"
#include "BridgeRegression.h"

// #include "HmcSampler.h"

#include <iostream>
#include <string>
#include <GetPot>

using std::cout;
using std::string;

int main(int argc, char** argv)
{

  GetPot cl(argc, argv);

  string y_file = cl.follow("", "-y");
  string X_file = cl.follow("", "-X");
  string fpre   = cl.follow("prefix", "--name");

  if (y_file.empty() || X_file.empty()) {
    printf("Usage: test_gibbs -y y.db -X X.db [--stable] [--ortho] [--samp 100000] [--burn 10000]\n");
    return 0;
  }

  bool do_stb = cl.search("--stable");
  bool ortho  = cl.search("--ortho");

  uint M = cl.follow(100000, "--samp");
  uint burn = cl.follow(10000, "--burn");

  // int d = 3;
  // int s = 1;
  // HmcSampler hmc(d, s);

  // The data.
  // Matrix X; X.read("X.lars", true);
  // Matrix X; X.read("Q.lars", true);
  // Matrix y; y.read("Y.lars", true);

  Matrix X; X.load(X_file, true);
  Matrix y; y.load(y_file, true);  

  // // True data.
  // Matrix beta_data; beta_data.read("beta.data", true);

  uint P = X.cols();
  uint N = X.rows();

  printf("N: %i, P: %i\n", N, P);

  double my_alpha = 0.5;

  // Least squares.
  Matrix XX; mult(XX, X, X, 'T', 'N');
  Matrix ls; mult(ls, X, y, 'T', 'N');

  symsolve(XX, ls);
  // cout << "LS:\n" << ls << "\n";

  Matrix beta(P      , (uint)1, M);

  Matrix u    (P, 1, M, 0.0);
  Matrix omega(P, 1, M, 2.0);
  Matrix shape(P, 1, M, 0.0);

  Matrix sig2((uint)1, (uint)1, M);
  Matrix tau ((uint)1, (uint)1, M);
  Matrix alpha((uint)1, (uint)1, M);

  double rt= 0;

  if (!do_stb && !ortho) 
    rt = bridge_regression(beta, u, omega, shape, sig2, tau, alpha, y, X, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, my_alpha, burn);
  if (!do_stb && ortho)
    rt = bridge_regression_ortho(beta, u, omega, shape, sig2, tau, alpha, y, X, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, my_alpha, burn);
  
  Matrix lambda(N, (uint)1, M);

  if (do_stb && !ortho) 
    rt = bridge_regression_stable(beta, lambda, sig2, tau, alpha, y, X, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, my_alpha, burn);
  if (do_stb && ortho)
    rt = bridge_regression_stable_ortho(beta, lambda, sig2, tau, alpha, y, X, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, 0.0, 0.0, my_alpha, burn);

  printf("Regression done.\n");
  // beta.write(fpre.append("beta.post"));

  beta.reshape(P, M);
  cout << "LS PMean" << "\n";
  cout << ls.cbind(rowMeans(beta)) << "\n";
  
  Matrix runtime(1);
  runtime(0) = rt;

  // u.write("u.post");
  // omega.write("omega.post");
  // tau.write("tau.post");

  // RNG r;
  // Matrix blah(100000);
  // for(uint i = 0; i < blah.vol(); i++)
  //   blah(i) = r.tnorm(2.0, 3.0);
  // blah.write("tn.data");

  // // EM
  // BridgeRegression br(X, y);
  // Matrix beta_EM(P);
  // double tau_seed = 2.0;

  // int iter_direct = br.EM(beta_EM, 50.0, tau_seed, alpha, 10e8, 10e-9, 1000, false);
  // cout << "Not CG...\n";
  // cout << "iters (dir): " << iter_direct << "\n";
  // cout << "beta (EM):\n" << beta_EM;

  // iter_direct = br.EM(beta_EM, 50.0, tau_seed, alpha, 10e8, 10e-9, 1000, true);
  // cout << "CG...\n";
  // cout << "iters (dir): " << iter_direct << "\n";
  // cout << "beta (EM):\n" << beta_EM;

  // // rtnorm
  // Matrix A(2, 2);
  // A << "2 1 0 1";
  // Matrix AA(A, A, 'T', 'N');
  // Matrix beta(2);
  // beta << "0 0";
  // Matrix bmean(2);
  // bmean << "1 0";
  // Matrix b(2);
  // b << "2 1";
  // double sig2 = 1.0;

  // int nsamp = 1;
  // Matrix samp(2, 1, nsamp);
  // samp.fill(0.0);

  // RNG r;
  // BridgeRegression br2(AA, beta);

  // Matrix mn(2); mn.fill(2);

  // for(int i=0; i<nsamp; i++) {
  //   br2.rtnorm_gibbs(beta, bmean, AA, sig2, b, r);
  //   samp[i].copy(beta);
  //   mn(0) += beta(0);
  //   mn(1) += beta(1);
  // }
  
  // cout << "mean: " << hdiveq(mn, (double)nsamp);

  // samp.write("samp.txt");

}
