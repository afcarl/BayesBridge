## Bridge regression using the mixture of triangles formulation.

## library(msm)
## library(mvtnorm)
## library(truncnorm)

draw.u <- function(tau, beta, w, alpha)
{
  m = 1-{abs(beta/tau)*w^(-1/alpha)}
  runif(length(beta),max=m)
}

draw.w <- function(tau, beta, u, alpha)
{
  p = length(beta)
  a = (abs(beta/tau)/(1-u))^alpha
  pr = alpha/(1-alpha*a)
  p1 = 0.5*(1+alpha)
  pr = p1/(1-p1*a)
  shape = (runif(p) < pr) + 1;
  w = rgamma(p, shape, 1);
  w+a
}

draw.tau <- function(beta, alpha, c, d)
{
  p = length(beta)
  nu = rgamma(1, c + p/alpha, rate=d + sum(abs(beta)^alpha))
  tau = nu^(-1/alpha)
  return(tau);
}

draw.sig2 <- function(beta, x, y, sig2.shape=0.0, sig2.scale=0.0)
{
  n = length(y)
  rss = sum( (as.matrix(y)-x%*%as.matrix(beta))^2 )
  prec = rgamma(1, sig2.shape+n/2, rate=sig2.scale+rss/2)
  return(1/prec)
}

## Geweke style Gibbs sampling
draw.beta.1 <- function(beta, bhat, xx, sig2, tau, u, w, alpha)
{
  p = length(bhat)
  b = (1-u)*{w^(1/alpha)}*tau
  for ( i in 1:p )
  {
    m = bhat[i] - crossprod(xx[i,-i],beta[-i]-bhat[-i])/xx[i,i]
    v = sig2/xx[i,i]
    beta[i] = rtnorm(1,m,sqrt(v),-b[i],b[i])
  }
  beta
}

## Rodriguez-Yam style Gibbs sampling - using SVD
draw.beta.2 <- function(beta, a, tV, d, sig2, tau, u, w, alpha)
{
  P = length(beta);
  
  b = (1-u)*{w^(1/alpha)}*tau
  z = tV %*% beta;

  for (i in 1:P) {
    lmax = -Inf
    rmin =  Inf

    for (j in 1:P) {
      vji = tV[i,j];
      vj  = tV[ ,j];
      ## rji = vj %*% z - vji * z[i];
      rji = vj[-i] %*% z[-i];
      Dif = b[j] - rji
      Sum = b[j] + rji
      left  = ifelse(vji > 0, -Sum, -Dif) / abs(vji);
      right = ifelse(vji > 0,  Dif,  Sum) / abs(vji);
      lmax = max(c(lmax,left))
      rmin = min(c(rmin,right))
    }

    if (d[i]!=0) {
      m = a[i] / (d[i]^2);
      s = sqrt(sig2) / d[i];
      z[i] = rtnorm(1, m, s, lmax, rmin);
    } else {
      cat("d =", d[i], "\n");
      z[i] = runif(1, lmax, rmin);
    }

  }

  beta = t(tV) %*% z;
  beta
}

## Alternate Rodriguez-Yam style Gibbs sampling.  Not as good.
draw.beta.3 <- function(beta, bhat, xx, sig2, tau, u, w, alpha)
{
  P = length(beta)
  evd = eigen(xx)
  rt  = evd$vectors %*% diag(sqrt(evd$values), P) %*% t(evd$vectors);  ## sqrt XX
  irt = evd$vectors %*% diag(sqrt(1/evd$values), P) %*% t(evd$vectors);
  b = (1-u)*{w^(1/alpha)}*tau
  
  z = rt %*% beta
  m = rt %*% bhat

  for (i in 1:P) {
  
    left  = rep(0, P)
    right = rep(0, P)
    
    for (j in 1:P) {
      rji = irt[j,-i] %*% z[-i]
      Dif = b[j] - rji;
      Sum = b[j] + rji;
      left[j]  = ifelse(irt[j,i] > 0, -Sum, -Dif) / abs(irt[j,i]);
      right[j] = ifelse(irt[j,i] > 0,  Dif,  Sum) / abs(irt[j,i]);
    }
    
    lmax = max(left)
    rmin = min(right)

    z[i] = rtnorm(1, m[i], sqrt(sig2), lmax, rmin);
    
  }
  
  beta = irt %*% z;
}

################################################################################

bridge.tmix.R <- function(y, X, nsamp, alpha=0.5, sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,
                         burn=100, sig2=0.0, tau=0.0, verbose=500)
{
  X <- as.matrix(X)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ixx <- chol2inv(chol(xx))

  known.sig2 = sig2 > 0
  known.tau  = tau > 0

  jsvd = svd(X);
  tV = t(jsvd$v);
  d  = jsvd$d
  A  = jsvd$u %*% diag(d);
  a  = t(A) %*% y;
  
  p <- ncol(X)

  bhat <- drop(ixx%*%xy)
  beta <- rep(1,p)
  beta <- bhat
  ## tau <- 2
  w = rep(2,p)
  u = rep(0,p)

  #beta = bhat + rnorm(1,0,1)
  #w = 2*abs(beta)
  #tau=1

  if (sig2 <= 0) sig2 = (1/length(y))*sum((y-X%*%bhat)^2)
  if (tau  <= 0) tau  = 1;
  
  output <- list(u = matrix(nrow=nsamp, ncol=length(beta)),
                 w = matrix(nrow=nsamp, ncol=length(beta)),
                 beta = matrix(nrow=nsamp, ncol=length(beta)),
                 sig2 = rep(0, nsamp),
                 tau = rep(0, nsamp)
                 )

  start.time = proc.time();
  
  for( i in 1:(nsamp+burn))
  {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")
      if( i==(burn+1) ) ess.time = proc.time();

      if (!known.tau) tau = draw.tau(beta, alpha, nu.shape, nu.rate)

      if (!known.sig2) sig2 = draw.sig2(beta, X, y, sig2.shape, sig2.scale)

      w <- draw.w(tau, beta, u, alpha)
      u <- draw.u(tau, beta, w, alpha)

      ##beta = draw.beta.1(beta, bhat, xx, sig2, tau, u, w, alpha)
      beta = draw.beta.2(beta, a, tV, d, sig2, tau, u, w, alpha)
      ## beta = draw.beta.3(beta, bhat, xx, sig2, tau, u, w, alpha)

      if(i > burn)
      {
          output$u[i-burn,]    = u
          output$w[i-burn,]    = w
          output$beta[i-burn,] = beta
          output$sig2[i-burn]  = sig2
          output$tau[i-burn]   = tau
      }
  }

  end.time = proc.time();
  output$runtime = (end.time - start.time)[1];

  output
}

bridge.EM.R = function(y, X, alpha, ratio=1.0, lambda.max=1e9*ratio, tol=1e-9, max.iter=30){
    X <- as.matrix(X)
    xx <- t(X)%*%X
    xy <- t(X)%*%y
    ixx <- chol2inv(chol(xx))

    p <- ncol(X)

    bhat <- drop(ixx%*%xy)
    Beta = bhat

    diff = 1

    # tau = (Nu)^{-1/alpha}
    sigma = 1;
    tau = ratio;

    iter = 0

    while(diff > tol && iter < max.iter)
    {
        # YHat = X %*% Beta
        # sigma = sqrt(sum( (Y-YHat)^2 )/(n-p))
        # sigma=1
        # EXPECTATION STEP
        Lambda = pmin( alpha*(tau^(2-alpha)) * abs(Beta)^(alpha-2), lambda.max)
        #OmegaInv = as.numeric((d+1)/(d*sigma^2+(Y-YHat)^2))
        # H = solve((1/tau^2)*diag(as.numeric(Lambda))+t(X) %*% X) %*% t(X)
        # BetaNew = H %*% y
        BetaNew = solve((1/tau^2)*diag(as.numeric(Lambda))+t(X) %*% X, xy);
        # S = X %*% H
        diff = sum(abs(Beta - BetaNew))
        Beta = BetaNew
        # print(Beta);
        #Nu = (b.nu + sum(abs(Beta)/sigma))/(p + a.nu - 1)
        iter = iter + 1;
    }

    Beta = drop(Beta)
}
