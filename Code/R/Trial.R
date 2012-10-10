## Here we try to test whether the MCMCs are, in fact, producing the same results.

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

if (is.loaded("Bridge.so")) dyn.unload("Bridge.so");
if (!is.loaded("Bridge.so")) dyn.load("Bridge.so");
source("~/RPackage/BayesBridge/Code/C/BridgeWrapper.R");

source("BridgeTMix.R")
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

nsamp = 10000
burn  = 2000
## tau = 41
tau = 0.0
alpha = 0.5
ortho=FALSE
multi = 1;
 
tri.list = list()
stb.list = list()
## R.list   = list()
i = 0

idc = (i+1):(i+2); for (i in idc) {

  gb.tri = bridge.reg.tri(y, X, multi*nsamp, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate, burn=multi*burn, tau=tau, ortho=ortho)
  gb.stb = bridge.reg.stb(y, X, nsamp, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate, burn=burn, tau=tau, ortho=ortho)
  ## gb.R   = bridge.tmix.R(y, X, nsamp, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate, burn=burn, tau=tau, verbose=500);
  
  ss.tri = sum.stat(gb.tri)
  ss.stb = sum.stat(gb.stb)
  ## ss.R   = sum.stat(gb.R)

  of.interest = c(1,2,3,7,8,9);
  tri.list[[i]] = ss.tri[,of.interest]
  stb.list[[i]] = ss.stb[,of.interest]
  ## R.list[[i]]   = ss.R  [,of.interest]

}

nr = nrow(tri.list[[1]]);
nc = ncol(tri.list[[1]]);
nl = length(tri.list);

tri.array = array(dim=c(nr, nc, nl));
colnames(tri.array) = colnames(tri.list[[1]]);
stb.array = tri.array
## R.array = tri.array
for (j in 1:nl) {
  tri.array[,,j] = tri.list[[j]]
  stb.array[,,j] = stb.list[[j]]
  ## R.array[,,j]   = R.list[[j]]
}

apply(tri.array[,1,], 1, mean)
apply(stb.array[,1,], 1, mean)

################################################################################
                                 ## APPENDIX ##
################################################################################

## Note to self: above is a decent strategy.  You can keep track of many
## collecting summary statistics in an on-going fashion.  This is useful because
## it may be the case that it takes a long time for the MCMC to converge.  This
## appears to be the case in the above example.  For only 5 MCMCs of 8 * 10,000
## and 10,000 samples each it is unclear that they are doing the same thing.
## However, after 20 runs the mean of the sample means seem to be the same,
## which was not obvious after only 5 runs.

## I run the triangle method for 8 times longer to try and make each run result
## in roughly the same number of effective samples.  I want to make sure that
## the minimum effective sample sizes are roughly the same.

## Remember: the samples can look different because one is mixing much better.
