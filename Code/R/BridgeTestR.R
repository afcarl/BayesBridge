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

if (FALSE) {

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

  out.tri = bridge.tmix.R(y, X, nsamp, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate,
    burn=burn, sig2=0.0, tau=tau, verbose=500)

  out.nrm = bridge.nmix.R(y, X, nsamp, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate,
    burn=burn, sig2=0.0, tau=tau, verbose=500)

  stat.tri = sum.stat(out.tri)
  stat.nrm = sum.stat(out.nrm)

  stat.tri
  stat.nrm

  ##---------------------------------------------------------------------------
  ## Looking for multimodality.

  P = ncol(X);
  bk = 60;
  
  for (i in 1:P) {
    beta.1 = out.tri$beta[out.tri$shape[,i]==1,i]
    beta.2 = out.tri$beta[out.tri$shape[,i]==2,i]
    h1 = hist(beta.1, breaks=bk)
    h2 = hist(beta.2, breaks=bk)
    ymax=max(c(h1$counts, h2$counts));
    xmin = min(c(h1$breaks, h2$breaks));
    xmax = max(c(h1$breaks, h2$breaks))
    plot(h1, col="#FF000088", ylim=c(0,ymax))
    plot(h2, col="#0000FF66", add=TRUE)
    readline("<ENTER>");
  }
  
}
