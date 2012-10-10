## This doesn't seem to work perfectly: detach("package:BayesBridge")

library("coda")
## source("~/RPackage/BayesLogit/Code/R/Efficiency.R")

################################################################################
                                  ## Setup ##
################################################################################

run <- list("EFRON" = FALSE, # Efron's diabetes data
            "BH"    = FALSE, # Boston Housing
            "BHI"   = FALSE, # Boston Housing with interactions
            "NIR"   = FALSE) # NIR

## Ortogonalizing matrices by QR.
oth <- list("EFRON" = FALSE,
            "BH"    = FALSE,
            "BHI"   = FALSE)

## RUN INFO
nsamp = 10000
burn = 2000
alpha = 0.5
ntrials = 2
tau = 0 ## Set to <= 0 for unknown tau.
betaburn = 0
use.hmc = FALSE

save.it  = FALSE ## Write output to file
plot.it  = FALSE ## Plot histograms
print.it = TRUE  ## Print summary.

################################################################################
                           ## Load .so or package ##
################################################################################

use.library=TRUE
if (use.library) {
  ## library("BayesBridge", lib.loc="~/RPackage/BayesBridge/Code/BBPackage/Test/");
  library("BayesBridge");
} else {
  if (is.loaded("Bridge.so")) dyn.unload("Bridge.so");
  if (!is.loaded("Bridge.so")) dyn.load("Bridge.so");
  source("~/RPackage/BayesBridge/Code/C/BridgeWrapper.R");
}

################################################################################
                     ## Make sure everything is working ##
################################################################################

## TRIAL RUN ##

if (FALSE) {
  library("BayesBridge")
  data(diabetes, package="BayesBridge");
}

if (FALSE) {
  dyn.unload("Bridge.so")
  if (!is.loaded("Bridge.so")) dyn.load("~/RPackage/BayesBridge/Code/R/Bridge.so");
  source("~/RPackage/BayesBridge/Code/C/BridgeWrapper.R");
}

if (FALSE) {

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

  gb1 = bridge.reg.tri(y, X, nsamp, 0.5, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, burn, ortho=FALSE, betaburn=betaburn, use.hmc=use.hmc);
  gb2 = bridge.reg.stb(y, X, nsamp, 0.5, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, burn, ortho=FALSE);

  sstat.1 = sum.stat(gb1);
  sstat.2 = sum.stat(gb2);

  sstat.1
  sstat.2

  ## Make sure things look right.
  for(i in 1:10){
    par(mfrow=c(1,2));
    hist(gb1$beta[,i], breaks=100, prob=TRUE);
    ## hist(gb2$beta[,i], breaks=100, prob=TRUE);
    hist(gb2$beta[,i], breaks=100, prob=TRUE);
    readline("Press Enter...");
  }

}

################################################################################

################################################################################

## Generates summary statistics.
sum.stat <- function(gbs, thin=1)
{
  n = nrow(gbs$beta);
  p = ncol(gbs$beta);
  gbs$beta = gbs$beta[seq(1, n, thin),];
  if (p==1) {gbs$beta = as.matrix(gbs$beta);}
  sstat = matrix(nrow=p, ncol=7);

  sstat[,1] = apply(gbs$beta, 2, mean);
  sstat[,2] = apply(gbs$beta, 2, sd);
  sstat[,3] = apply(gbs$beta, 2, effectiveSize);
  sstat[,4] = sstat[,3] / gbs$runtime;
  ## sstat[,5] = apply(gbs$beta, 2, ESS);
  ## sstat[,6] = sstat[,5] / gbs$runtime;
  sstat[,5] = sstat[,1] / sstat[,2];
  sstat[,6] = apply(gbs$beta, 2, function(x){quantile(x, 0.1)});
  sstat[,7] = apply(gbs$beta, 2, function(x){quantile(x, 0.9)});

  ## colnames(sstat) = c("mean", "sd", "ESS", "ESS.sec", "myESS", "myESS.sec", "t", "Q.1", "Q.9");
  colnames(sstat) = c("mean", "sd", "ESS", "ESS.sec", "t", "Q.1", "Q.9");

  sstat
}

################################################################################
                             ## Compare Routines ##
################################################################################

compare.it <- function(y, X, nsamp=10000,
                       alpha=0.5,
                       sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,
                       tau=0.0, 
                       burn=0, ntrials=1, ortho=FALSE, betaburn=0, use.hmc=FALSE)
{
  ## Returns a list of summary statistics, the last Gibbs MCMC, and a
  ## table with the average statistics over several runs.
  
  tri.info = matrix(nrow=ntrials, ncol=5);
  colnames(tri.info) = c("ave.ESS", "sd.ESS", "ave.ESS.sec", "sd.ESS.sec", "runtime");
  stb.info = tri.info;

  tri.stat = list();
  stb.stat = list();

  for (i in 1:ntrials){
    
    gb.tri = bridge.reg.tri(y, X, nsamp, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate,
      0.0, tau, burn, ortho=ortho, betaburn=betaburn, use.hmc=use.hmc);
    gb.stb = bridge.reg.stb(y, X, nsamp, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate,
      0.0, tau, burn, ortho=ortho)
    
    sstat.tri = sum.stat(gb.tri);
    sstat.stb = sum.stat(gb.stb);

    tri.stat[[i]] = sstat.tri;
    stb.stat[[i]] = sstat.stb;

    tri.info[i,"ave.ESS"] = mean(sstat.tri[,"ESS"]);
    tri.info[i,"sd.ESS" ] = sd(sstat.tri[,"ESS"]);
    tri.info[i,"ave.ESS.sec"] = mean(sstat.tri[,"ESS.sec"]);
    tri.info[i,"sd.ESS.sec"]  = sd(sstat.tri[,"ESS.sec"]);
    tri.info[i,"runtime"] = gb.tri$runtime

    stb.info[i,"ave.ESS"] = mean(sstat.stb[,"ESS"]);
    stb.info[i,"sd.ESS" ] = sd(sstat.stb[,"ESS"]);
    stb.info[i,"ave.ESS.sec"] = mean(sstat.stb[,"ESS.sec"]);
    stb.info[i,"sd.ESS.sec"]  = sd(sstat.stb[,"ESS.sec"]);
    stb.info[i,"runtime"] = gb.stb$runtime
  }

  OUT <- list("tri.info"=tri.info, "stb.info"=stb.info,
              "tri.list"=tri.stat, "stb.list"=stb.stat,
              "tri.stat"=simplify2array(tri.stat), "stb.stat"=simplify2array(stb.stat),
              "tri.gb"=gb.tri, "stb.gb"=gb.stb);

  OUT
}

##------------------------------------------------------------------------------

plot.info <- function(info, cnames, P=ncol(info$tri.gb$beta))
{
  ## sig2
  par(mfrow=c(1,2));
  hist(info$tri.gb$sig2, breaks=100, prob=TRUE, main=paste("tri: sig2"));
  hist(info$stb.gb$sig2, breaks=100, prob=TRUE, main=paste("stb: sig2"));
  cat("sig2: \n");
  cat("tri:", mean(info$tri.gb$sig2), sd(info$tri.gb$sig2), "\n");
  cat("stb:", mean(info$stb.gb$sig2), sd(info$stb.gb$sig2), "\n");
  readline("Press Enter...");

  ## tau
  if (!is.null(info$tri.gb$tau)) {
    par(mfrow=c(1,2));
    hist(info$tri.gb$tau, breaks=100, prob=TRUE, main=paste("tri: tau"));
    hist(info$stb.gb$tau, breaks=100, prob=TRUE, main=paste("stb: tau"));
    cat("tau: \n");
    cat("tri:", mean(info$tri.gb$tau), sd(info$tri.gb$tau), "\n");
    cat("stb:", mean(info$stb.gb$tau), sd(info$stb.gb$tau), "\n");
    readline("Press Enter...");
  }

  ## beta
  for(i in 1:P){
    par(mfrow=c(1,2));
    hist(info$tri.gb$beta[,i], breaks=100, prob=TRUE, main=paste("tri:", cnames[i]));
    hist(info$stb.gb$beta[,i], breaks=100, prob=TRUE, main=paste("stb:", cnames[i]));
    cat("Variable", i, "\n");
    cat("tri:", mean(info$tri.gb$beta[,i]), sd(info$tri.gb$beta[,i]), "\n");
    cat("stb:", mean(info$stb.gb$beta[,i]), sd(info$stb.gb$beta[,i]), "\n");
    readline("Press Enter...");
  }
}

##------------------------------------------------------------------------------

run.it <- function(y, X, nsamp=1000,  burn=100, 
                   alpha=0.5, sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,
                   ntrials=1, tau=0.0,
                   save.it=FALSE, print.it=FALSE, plot.it=FALSE, name="somerun", ortho=FALSE, betaburn=0, use.hmc=FALSE)
{
  ## Runs the comparison and plots/prints/saves the resulting data.
  
  info = compare.it(y, X, nsamp=nsamp,
    alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, nu.shape, nu.rate, tau=tau,
    burn=burn, ntrials=ntrials, ortho=ortho, betaburn=betaburn, use.hmc=use.hmc)

  ## Make sure things look right.
  if (plot.it) {
    cnames = colnames(X);
    plot.info(info, cnames);
  }

  if (print.it) {
    print("tri.stat:")
    print(info$tri.stat[[ntrials]])
    print("stb.stat:")
    print(info$stb.stat[[ntrials]])
  }

  cat("Info size:", object.size(info) / 2^20, "Mb\n")
  print("tri:")
  print(info$tri.info)
  print("stb:")
  print(info$stb.info)

  filename = paste("info", name, "RData", sep=".");
  if (save.it) save(info, file=filename, compress=TRUE);

  info
}

################################################################################
                              ## SUITE OF TESTS ##
################################################################################

##------------------------------------------------------------------------------
## Efron's diabetes data.

if (run$EFRON) {

  if (tau>0) tau = 41
  
  ## Load data.
  ## load("~/RPackage/BayesBridge/Code/C/diabetes.RData")
  data("diabetes", package="BayesBridge");
  cov.name = colnames(diabetes$x);
  y = diabetes$y;
  X = diabetes$x;
  cnames = cov.name
  
  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }

  info <- run.it(y, X, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=print.it, plot.it=plot.it, name="efron",
                 ortho=FALSE, betaburn=betaburn, use.hmc=use.hmc)

}

##------------------------------------------------------------------------------
## Boston Housing

if (run$BH) {

  if (tau>0) tau   = 0.15
  
  data("BostonHousing", package="mlbench")

  ## Setup
  y = BostonHousing$medv

  ## No interactions or squared terms.
  X = model.matrix(medv ~ ., data=BostonHousing);

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }
  X = X[,-1]
  
  cnames = colnames(X);

  info <- run.it(y, X, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=print.it, plot.it=plot.it, name="BH",
                 ortho=FALSE, betaburn=betaburn, use.hmc=use.hmc)  

}

##------------------------------------------------------------------------------
## Boston Housing with interactions

if (run$BHI) {

  if (tau>0) tau = 0.15
  
  data("BostonHousing", package="mlbench")

  ## Setup
  y = BostonHousing$medv

  ## No interactions or squared terms.
  X = model.matrix(medv ~ ., data=BostonHousing);
  cnames = colnames(X);

  ## Include interations and squared terms.
  idc  = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13);
  nidc = length(idc)

  X.base  = model.matrix(medv ~ . * ., data=BostonHousing);

  X = X.base
  cnames = colnames(X.base);
  for (i in 1:nidc){
    nm = paste(colnames(BostonHousing)[idc[i]], 2, sep="");
    cnames = c(cnames, nm);
    X = cbind(X, BostonHousing[,idc[i]] * BostonHousing[,idc[i]]);
  }

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }
  X = X[,-1]

  cnames = colnames(X);

  info <- run.it(y, X, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=FALSE, plot.it=plot.it, name="BHI",
                 ortho=FALSE, betaburn=betaburn, use.hmc=use.hmc)
  
}

##------------------------------------------------------------------------------
## Ozone

if (FALSE) {

  data("Ozone", package="mlbench")
  ## Has missing data.
  
}


##------------------------------------------------------------------------------
## NIR Glucose - has more predictors than observations

if (run$NIR) {

  ## RUN INFO
  nsamp = 10000
  burn  = 2000
  alpha = 0.5
  tau   = 100
  ntrials = 1

  ## Bad situation: 166 rows, 236 columns in design matrix.
  
  data("NIR", package="chemometrics")

  y = NIR$yGlcEtOH[,1]
  X = cbind(1.0, NIR$xNIR);
  X = as.matrix(X);

  info <- run.it(y, X, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=TRUE, plot.it=plot.it, name="NIR",
                 ortho=FALSE, betaburn=betaburn, use.hmc=use.hmc)

}

################################################################################
                         ## ORTHOGONAL DESIGN MATRIX ##
################################################################################

##------------------------------------------------------------------------------
## Efron - orthogonalized

if (oth$EFRON) {

  if (tau>0) tau = 41;
  
  ## Load data.
  ## load("~/RPackage/BayesBridge/Code/C/diabetes.RData")
  data("diabetes", package="BayesBridge");
  y = diabetes$y;
  X = diabetes$x;

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }

  ## Ortongonalize
  Q = qr.Q(qr(X))
  colnames(Q) = c(1:ncol(Q))

  info <- run.it(y, Q, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=TRUE, plot.it=plot.it, name="efron-ortho", ortho=TRUE)
  
}

##------------------------------------------------------------------------------
## Boston Housing - orthogonalized

if (oth$BH) {

  if (tau>0) tau = 0.15
  
  data("BostonHousing", package="mlbench")

  ## Setup
  y = BostonHousing$medv

  ## No interactions or squared terms.
  X = model.matrix(medv ~ ., data=BostonHousing);

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }
  X = X[,-1]

  ## Ortongonalize
  Q = qr.Q(qr(X))
  colnames(Q) = c(1:ncol(Q))

  info <- run.it(y, Q, nsamp=nsamp, burn=burn,
                        alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                        save.it=save.it, print.it=TRUE, plot.it=plot.it, name="BH-ortho", ortho=TRUE)  

}

##------------------------------------------------------------------------------
## Boston Housing with interactions

if (oth$BHI) {

  if (tau>0) tau=0.15
  
  data("BostonHousing", package="mlbench")

  ## Setup
  y = BostonHousing$medv

  ## No interactions or squared terms.
  X = model.matrix(medv ~ ., data=BostonHousing);
  cnames = colnames(X);

  ## Include interations and squared terms.
  idc  = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13);
  nidc = length(idc)

  X.base  = model.matrix(medv ~ . * ., data=BostonHousing);

  X = X.base
  cnames = colnames(X.base);
  for (i in 1:nidc){
    nm = paste(colnames(BostonHousing)[idc[i]], 2, sep="");
    cnames = c(cnames, nm);
    X = cbind(X, BostonHousing[,idc[i]] * BostonHousing[,idc[i]]);
  }

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }
  X = X[,-1]

  ## Ortongonalize
  Q = qr.Q(qr(X))
  colnames(Q) = c(1:ncol(Q))

  info <- run.it(y, Q, nsamp=nsamp, burn=burn,
                        alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                        save.it=save.it, print.it=FALSE, plot.it=plot.it, name="BHI-ortho", ortho=TRUE)
  
}

################################################################################
                                ## RUN STUFF ##
################################################################################

if (FALSE) {

  ## Known Tau
  temp.gb.tri = bridge.reg.tri(y, X, nsamp=nsamp,
    alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, tau=tau, burn=burn);
  temp.gb.stb = bridge.reg.stb(y, X, nsamp=nsamp,
    alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, tau=tau, burn=burn);

  cnames=colnames(X)
  
}

################################################################################
                                 ## APPENDIX ##
################################################################################

## Fall, 2011

## When using R for the random number generation the routines RUN FASTER when
## compiling the package, as opposed to compiling a .so file by hand and then
## using that.  Not sure why that is.

## It appears that the number of samples I take also affects the relative speed
## of each routine.  If you try nsamp = 1e6 vs nsmap=1e5 there is a difference
## in ESS/sec.  There is not a difference in ESS.

## It think something must be off on my bridge with triangles routine.  The
## posterior means are not in complete agreement with the stable method.

## Spring, 2012

## It appears that with the Efron data the triangular mixture is competitive.
## However, with the Boston Housing data the triangular mixture is not
## competitive.

## Both the Boston Housing data with interactions and squared terms and the NIR
## data (which has P > N), have odd posterior behavior when computed using the
## mixture of triangles.  At times it appears that the posterior has been stuck
## in several local modes.  The normal mixture by stables has well behaved
## posteriors.  This is reflected, I think, in the ESS.  If you look at the ESS
## of specific coordinates of beta, the odd looking posterior distributions will
## have extremely small ESS.  Like in the single digits out of 10000 or 100000
## draws.  This is interesting to me because it provides an example of when MCMC
## does not work well.

## Fall, 2012

## It now appears that the best we can do is use the rotation method for
## sampling the truncated normal.

## Note to self: it is helpful to have a fully C routine for posterior
## estimation.  In that case you can find out where the code is the slowest.  I
## was able to go from code that was slower than the stable method to a code
## that was faster by analyzing that data.

## Originally, we thought about maybe using orthogonal design matrices to "win".
## THe problem is that, in that case, one would have an orthogonal posterior
## variance using the stable method, which would make it run fast enough to beat
## the triangle method.

## I'm not how sure much we get hurt by including the slice variables.
## Unfortunately, I don't see any easy way to get rid of them.  If examines the
## conditional posterior $(beta | \omgea, sig2, tau)$, then he or she will see
## that it is a truncated normal tilted by some complicated peicewise
## polynomial--not easy to simulate.

## The reality is that the stable method works well.  You see how that method is
## advantageous if you look at the NIR data.  It appears that the stable method
## mixes much better.  But this is a degenerate situation as you are pretty much
## getting nothing but prior.  There are a few instances where it seems like
## there may be some weight not on zero.
