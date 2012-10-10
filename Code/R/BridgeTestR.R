## To test our R files.

## Generates summary statistics.
sum.stat <- function(gbs, thin=1)
{
  n = nrow(gbs$beta);
  p = ncol(gbs$beta);
  gbs$beta = gbs$beta[seq(1, n, thin),];
  if (p==1) {gbs$beta = as.matrix(gbs$beta);}
  sstat = matrix(nrow=p, ncol=9);

  sstat[,1] = apply(gbs$beta, 2, mean);
  sstat[,2] = apply(gbs$beta, 2, sd);
  sstat[,3] = apply(gbs$beta, 2, effectiveSize);
  sstat[,4] = sstat[,3] / gbs$runtime;
  sstat[,5] = apply(gbs$beta, 2, ESS);
  sstat[,6] = sstat[,5] / gbs$runtime;
  sstat[,7] = sstat[,1] / sstat[,2];
  sstat[,8] = apply(gbs$beta, 2, function(x){quantile(x, 0.1)});
  sstat[,9] = apply(gbs$beta, 2, function(x){quantile(x, 0.9)});

  colnames(sstat) = c("mean", "sd", "ESS", "ESS.sec", "myESS", "myESS.sec", "t", "Q.1", "Q.9");

  sstat
}

################################################################################

if (TRUE) {

  source("~/RPackage/BayesBridge/Code/R/BridgeTMix.R")
  source("~/RPackage/BayesBridge/Code/R/BridgeNMix.R")
  source("~/RPackage/BayesLogit/Code/R/Efficiency.R")

  library("coda")

  sig2.shape = 0.0
  sig2.scale = 0.0
  nu.shape = 2.0
  nu.rate  = 2.0

  load("~/RPackage/BayesBridge/Code/C/diabetes.RData")
  cov.name = colnames(diabetes$x);
  y = diabetes$y;
  X = diabetes$x;

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }

  LS = solve(t(X) %*% X, t(X) %*% y);

  nsamp = 5000
  burn  = 1000
  tau = 0.0
  alpha = 0.5

  out.tri = bridge.tmix.R(y, X, nsamp, alpha=0.5, sig2.shape, sig2.scale, nu.shape, nu.rate,
    burn=burn, sig2=0.0, tau=tau, verbose=500)

  out.nrm = bridge.nmix.R(y, X, nsamp, alpha=0.5, sig2.shape, sig2.scale, nu.shape, nu.rate,
    burn=burn, sig2=0.0, tau=tau, verbose=500)

  stat.tri = sum.stat(out.tri)
  stat.nrm = sum.stat(out.nrm)

  stat.tri
  stat.nrm

}
