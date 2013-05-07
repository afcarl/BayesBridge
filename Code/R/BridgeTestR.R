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
  library("pgnorm")

  sig2.shape = 0.0
  sig2.scale = 0.0
  nu.shape = 2.0
  nu.rate  = 2.0
  alpha.a  = 1.0
  alpha.b  = 1.0

  nsamp = 10000
  burn  = 2000

  ## Diabetes
  load("~/RPackage/BayesBridge/Code/C/diabetes.RData")
  cov.name = colnames(diabetes$x);
  y = diabetes$y;
  X = diabetes$x;

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }

  ## Synthetic 1
  n = 200
  p = 20
  A <- matrix(0.99, nrow=p, ncol=p);
  diag(A) = 1.0
  U = chol(A);
  z = matrix(rnorm(n*p), nrow=p, ncol=n);
  X = t( t(U) %*% z );

  tau.syn = 1.0
  sig2.syn = 1.0
  alpha.syn = 0.95
  beta.syn = rep(1.0, p)
  ## beta.syn = rpgnorm(p, alpha.syn, 0.0, sig.for.pg(tau.syn, alpha.syn))
  y = X %*% beta.syn + sig2.syn * rnorm(n);
  
  ## Synthetic 2
  n = 200
  p = 100
  A = matrix(0.9, ncol=p, nrow=p); diag(A) = 1.0
  U = chol(A);
  z = matrix(rnorm(n*p), nrow=n, ncol=p);
  X = z %*% U;

  tau.syn = 1e-11
  sig2.syn = 1.0
  alpha.syn = 0.1
  ## You can't generate data this way--marginal isn't stable, it is polynomial tilted stable.
  ## for(i in 1:p) lambda.syn[i] = 2 * retstable(0.5 * alpha.syn, 1.0, 0.0, method="LD");
  ## beta.syn = rnorm(p, 0.0, tau.syn / sqrt(lambda.syn));
  beta.syn = rpgnorm(p, alpha.syn, 0.0, sig.for.pg(tau.syn, alpha.syn))
  y = X %*% beta.syn + sig2.syn * rnorm(n);

  ## Synthetic 3
  n = 1000
  p = 500
  alpha.syn = 0.5

  ## Estimate
  LS = solve(t(X) %*% X, t(X) %*% y);

  nsamp = 10000
  burn  = 2000

  alpha = 0.9
  tau = 0.0
    
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

  alpha.a = 0.5
  alpha.b = 10
  
  alpha = 0.85
  tau = 0.1
  sig2 = 0.0
  
  out.C.tri = bridge.reg.tri(y, X, nsamp=nsamp*10, alpha=alpha, sig2.shape=sig2.shape, sig2.scale=sig2.scale,
    nu.shape=nu.shape, nu.rate=nu.rate, alpha.a=alpha.a, alpha.b=alpha.b,
    sig2.true=0.0, tau.true=tau, burn=burn)

  out.C.stb = bridge.reg.stb(y, X, nsamp*10, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate, alpha.a=alpha.a, alpha.b=alpha.b,
    sig2.true=0.0, tau.true=tau, burn=burn)

  P = ncol(X);
  bk = 100;

  out.C = out.C.stb
  for (i in 1:P) {
    hist(out.C$beta[,i], breaks=bk);
    readline("<ENTER>");
  }

  P = ncol(X);
  bk = 100;

  out = out.C.tri
  par(mfrow=c(2,1))
  for (i in 1:P) {
    beta.1 = out$beta[out$shape[,i]==1,i]
    beta.2 = out$beta[out$shape[,i]==2,i]
    h1 = hist(beta.1, breaks=bk, plot=FALSE)
    h2 = hist(beta.2, breaks=bk, plot=FALSE)
    hall = hist(out$beta[,i], breaks=bk, plot=FALSE);
    ymax=max(c(h1$counts, h2$counts));
    xmin = min(c(h1$breaks, h2$breaks));
    xmax = max(c(h1$breaks, h2$breaks))
    plot(h1, col="#FF000088", ylim=c(0,ymax), main=paste("hist", i))
    plot(h2, col="#0000FF66", add=TRUE)
    plot(hall);
    readline("<ENTER>");
  }
  par(mfrow=c(1,1));

  ################################################################################

  out = out.C.tri
  for (i in 1:P) {
    png(paste("Images/beta-", i,".png", sep=""), width=900, height=400);
    par(mfrow=c(1,2))
    beta.1 = out$beta[out$shape[,i]==1,i]
    beta.2 = out$beta[out$shape[,i]==2,i]
    h1 = hist(beta.1, breaks=bk, plot=FALSE)
    h2 = hist(beta.2, breaks=bk, plot=FALSE)
    hall = hist(out$beta[,i], breaks=bk, plot=FALSE);
    ymax=max(c(h1$counts, h2$counts));
    xmin = min(c(h1$breaks, h2$breaks));
    xmax = max(c(h1$breaks, h2$breaks))
    plot(hall, col="#22222244", border="#22222200", xlab=expression(beta),
         main=paste("Posterior draws of beta", i));
    plot(h1, col="#10101020", border="#10101000", ylim=c(0,ymax), xlab=expression(beta),
         main="Draws stratified by mixture component.")
    plot(h2, col="#20202020", border="#10101010", add=TRUE)
    ## readline("<ENTER>");
    dev.off();
  }
  par(mfrow=c(1,1));
 

  ##------------------------------------------------------------------------------

  sig2.shape = 0.0
  sig2.scale = 0.0
  nu.shape = 2.0
  nu.rate  = 2.0

  nsamp = 10000
  burn  = 2000

  alpha = 0.85
  tau = 0.1
  sig2 = 0.0
  
  out = bridge.reg.tri(y, X, nsamp=nsamp, alpha=alpha, sig2.shape=sig2.shape, sig2.scale=sig2.scale,
    nu.shape=nu.shape, nu.rate=nu.rate,
    sig2.true=0.0, tau.true=tau, burn=burn)
  
  png("beta-ex.png", width=900, height=200);
  par(mfrow=c(1,4), mar=c(4, 5, 4, 3))
  #### FONT SIZE ####
  cex = 2.0
  par(mfrow=c(1,4))
  for (i in c(2,3)) {
    beta.1 = out$beta[out$shape[,i]==1,i]
    beta.2 = out$beta[out$shape[,i]==2,i]
    h1 = hist(beta.1, breaks=bk, plot=FALSE)
    h2 = hist(beta.2, breaks=bk, plot=FALSE)
    hall = hist(out$beta[,i], breaks=bk, plot=FALSE);
    ymax=max(c(h1$counts, h2$counts));
    xmin = min(c(h1$breaks, h2$breaks));
    xmax = max(c(h1$breaks, h2$breaks))
    plot(hall, col="#22222244", border="#22222200", xlab="", ylab="", main="");
    title(main=paste("Posterior draws of beta", i), xlab=expression(beta), ylab="Freq.", cex.main=cex, cex.lab=cex);
    plot(h1, col="#10101020", border="#10101000", ylim=c(0,ymax), xlab="", ylab="", main="");
    title(main=paste("Stratified"), xlab=expression(beta), ylab="Freq.", cex.main=cex, cex.lab=cex);
    plot(h2, col="#20202020", border="#10101010", add=TRUE)
  }
  dev.off()

  ################################################################################

  alpha.a = 1.0
  alpha.b = 1.0
  
  alpha = 0
  tau = 0
  sig2 = 0.0

  nsamp = 100000
  
  ## out.C.tri = bridge.reg.tri(y, X, nsamp=nsamp*10, alpha=alpha, sig2.shape=sig2.shape, sig2.scale=sig2.scale,
  ##   nu.shape=nu.shape, nu.rate=nu.rate, alpha.a=alpha.a, alpha.b=alpha.b,
  ##   sig2.true=0.0, tau.true=tau, burn=burn)

  out.C.stb = bridge.reg.stb(y, X, nsamp, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate, alpha.a=alpha.a, alpha.b=alpha.b,
    sig2.true=0.0, tau.true=tau, burn=burn)

  alpha.k = 0.9
  out.k = bridge.reg.stb(y, X, nsamp, alpha.k, sig2.shape, sig2.scale, nu.shape, nu.rate, alpha.a=alpha.a, alpha.b=alpha.b,
    sig2.true=0.0, tau.true=tau, burn=burn)

  P = ncol(X);
  bk = 100;

  ## png("Images/alpha-bhi.png", width=600, height=300)
  par(mfrow=c(1,2))
  print(dim(X))
  hist(out.C.stb$alpha, main="", xlab="", ylab="")
  title(main="Posterior of alpha", xlab=expression(alpha), ylab="Freq.")
  acf(out.C.stb$alpha, main="", xlab="", ylab="")
  title(main="Autocorrelation of alpha", xlab="index", ylab="Autocorrelation")
  ## dev.off()

  P = ncol(X);
  bk = 100;
  
  for (i in 1:P) {
    png(paste("Images/alpha-known-unknown-", i, ".png", sep=""), width=600, height=300)
    par(mfrow=c(1,2))
    hist(out.C.stb$beta[,i], breaks=bk, main=paste(cnames[i], "(alpha unknown)"), xlab=expression(beta));
    hist(out.k$beta[,i], breaks=bk, main=paste(cnames[i], "(alpha=0.9)"), xlab=expression(beta));
    ## readline("<ENTER>");
    dev.off();
  }

  quantile(x, probs = seq(0, 1, 0.01))
           

  apply(out.C.stb$beta, 2, mean)
  apply(out.k$beta, 2, mean)

  bridge.EM.R(y, X, alpha=mean(out.C.stb$alpha), ratio=mean(out.C.stb$tau)/mean(sqrt(out.C.stb$sig)), max.iter=100);
  bridge.EM.R(y, X, alpha=alpha.k, ratio=mean(out.C.stb$tau)/mean(sqrt(out.C.stb$sig)), max.iter=100);
  
  
}

################################################################################
                                   ## TEST ##
################################################################################

if (FALSE) {

  ## source("BridgeByStable.R")
  ## source("bridge2-MCMC.R")
  ## source("Bridge.R")
  ## library("copula")
  library("coda")
  library("BayesBridge")
  library("chemometrics")
  
  ## data(diabetes)
  load("~/RPackage/BayesBridge/Code/C/diabetes.RData");
  lm1 = lm(y ~ x, data=diabetes);
  X = model.matrix(y ~ x, data=diabetes);
  y = diabetes$y

  data("NIR", package="chemometrics")
  y = NIR$yGlcEtOH[,1]
  X = cbind(1.0, NIR$xNIR);
  X = as.matrix(X);

  data("BostonHousing", package="mlbench")
  y = BostonHousing$medv
  X = model.matrix(medv ~ ., data=BostonHousing);

  nsamp = 4000;
  burn = 0

  ## sig2 = mean(lm1$res^2);
  sig2 = mean((y-PPM(X) %*% y)^2);
  tau = 100.0
  
  bridge.tri = bridge.reg.know.sig2(y, X, nsamp=nsamp, burn=burn, alpha=0.5, sig2=sig2, tau=tau);
  bridge.stb = bridge.by.nmix(y, X, nsamp=nsamp, burn=burn, alpha=0.5, sig2, tau=tau, verbose=1000);

  ## bridge.tri.R <- bridge.reg.know.sig2.R(y, X, nsamp=nsamp, alpha=0.5, burn=burn,
  ##                                        sig2=sig2, tau=tau, verbose=1000);

  bridge.tri.R <- bridge(y, X, nsamp=nsamp, alpha=0.5, burn=burn,
                         sig2=sig2, tau=tau, verbose=1000);

  bridge.em = bridge.EM(y, X, alpha=0.5, ratio=tau/sqrt(sig2));
  
  k = 3;
  sstat = matrix(nrow=ncol(X), ncol=k*3);
  
  sstat[,1] = apply(bridge.tri$beta, 1, mean);
  sstat[,4] = apply(bridge.stb$beta, 2, mean);
  sstat[,7] = apply(bridge.tri.R$beta, 2, mean);
  
  sstat[,2] = apply(bridge.tri$beta, 1, sd);
  sstat[,5] = apply(bridge.stb$beta, 2, sd);
  sstat[,8] = apply(bridge.tri.R$beta, 2, sd);

  sstat[,3] = effectiveSize(t(bridge.tri$beta))
  sstat[,6] = effectiveSize(bridge.stb$beta)
  sstat[,9] = effectiveSize(bridge.tri.R$beta)

  sstat

  for (i in 1:11) {
    par(mfrow=c(1,2))
    hist(bridge.tri$beta[i,], prob=TRUE, breaks=100, main=colnames(X)[i])
    abline(v=bridge.em[i], col=2)
    hist(bridge.stb$beta[,i], prob=TRUE, breaks=100)
    abline(v=bridge.em[i], col=2)
    readline("<ENTER>")
  }
  
}
