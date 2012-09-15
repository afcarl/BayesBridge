
bridge.stable <- function(y, x, niter, alpha, burn=100, verbose=500, sig2, tau=1.0)
{
  X <- as.matrix(x)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ## ixx <- chol2inv(chol(xx))

  p <- ncol(X)
  n <- length(y)
  gamma = alpha

  ##  w <- rgamma(p,1)
  ##  for ( i in 1:p )
  ##    if ( runif(1)>alpha )
  ##      w[i] <- rgamma(1,2)

  ## bhat <- drop(ixx%*%xy)
  beta <- rep(1,p)
  lambda = rep(1, p);
  ## tau <- 2

  ## beta = bhat
  ## sig2 = mean((y-X%*%bhat)^2)

  output <- list(lambda = matrix(nrow=niter, ncol=length(beta)),
                 beta   = matrix(nrow=niter, ncol=length(beta)),
                 sig2   = rep(0, niter)
                 )

  colnames(output$beta) = colnames(X);

  ## GIBBS
  for( i in 1:(niter+burn))
    {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")
      
      ## lambda
      for(j in 1:p)
        lambda[j] = 2 * retstable(0.5 * gamma, 1.0, beta[j]^2 / tau^2, method="LD");
      
      ## sig2 - before beta to set.
      rss = sum( (as.matrix(y)-x%*%as.matrix(beta))^2 )
      prec = rgamma(1, 2+n/2, 2+rss/2)
      sig2 = 1 / prec
      
      ## beta
      VInv = xx + diag(lambda * sig2 / tau^2, p);
      V = solve(VInv);
      U = chol(V);
      m = V %*% xy;

      beta = drop( m + sqrt(sig2) * t(U) %*% rnorm(p))

      if(niter > burn)
        {
          output$beta[i-burn,]   = beta
          output$lambda[i-burn,] = lambda
          output$sig2[i-burn]   = sig2
        }
    }

  output;
}

################################################################################
                                   ## TEST ##
################################################################################

if (FALSE) {

  ## source("BridgeByStable.R")
  ## source("bridge2-MCMC.R")
  ## source("Bridge.R")
  library("copula")
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
  bridge.stb = bridge.stable(y, X, niter=nsamp, burn=burn, alpha=0.5, sig2, tau=tau, verbose=1000);

  ## bridge.tri.R <- bridge.reg.know.sig2.R(y, X, nsamp=nsamp, alpha=0.5, burn=burn,
  ##                                        sig2=sig2, tau=tau, verbose=1000);

  bridge.tri.R <- bridge(y, X, niter=nsamp, alpha=0.5, burn=burn,
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
