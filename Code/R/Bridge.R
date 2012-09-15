# library(msm)
# library(mvtnorm)
# library(truncnorm)

draw.u <- function(tau, beta, w, alpha)
{
  m = 1-{abs(beta/tau)*w^(-1/alpha)}
  runif(length(beta),max=m)
}

draw.w <- function(tau, beta, u, alpha)
{
  p = length(beta)
  a = (abs(beta/tau)/(1-u))^alpha
  w = rgamma(p,1,1)
  pr = alpha/(1-alpha*a)
  #p1 = 0.5*(1+alpha)
  #pr = p1/(1-p1*a)
  for ( i in 1:p )
    if ( runif(1) < pr[i] )
      w[i] = rgamma(1,2,1)
  w+a
}

draw.tau <- function(beta, alpha, c, d)
{
	p = length(beta)
	nu = rgamma(1, c + p/alpha, rate=d + sum(abs(beta)^alpha))
	tau = nu^(-1/alpha)
	return(tau);
}

draw.sig2 <- function(beta, x, y)
{
	n = length(y)
	rss = sum( (as.matrix(y)-x%*%as.matrix(beta))^2 )
	prec = rgamma(1, 2+n/2, rate=2+rss/2)
	return(1/prec)
}

draw.beta <- function(beta, xx, bhat, sig2, tau, u, w, alpha)
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


bridge.reg.R <- function(y, X, nsamp, alpha=0.5, nu.shape=0.5, nu.rate=0.5,
                         burn=100, verbose=500)
{
  X <- as.matrix(X)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ixx <- chol2inv(chol(xx))

  p <- ncol(X)

  c = nu.shape;
  d = nu.rate;

  bhat <- drop(ixx%*%xy)
  beta <- rep(1,p)
  tau <- 2
  w = rep(2,p)

  #beta = bhat + rnorm(1,0,1)
  #w = 2*abs(beta)
  #tau=1

  sig2 = (1/length(y))*sum((y-X%*%bhat)^2)

  output <- list(u = matrix(nrow=nsamp, ncol=length(beta)),
                 w = matrix(nrow=nsamp, ncol=length(beta)),
                 beta = matrix(nrow=nsamp, ncol=length(beta)),
                 sig2 = rep(0, nsamp),
                 tau = rep(0, nsamp)
                 )

  for( i in 1:(nsamp+burn))
  {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")

      u <- draw.u(tau, beta, w, alpha)
      w <- draw.w(tau, beta, u, alpha)
      tau = draw.tau(beta, alpha, c, d)
      # tau=50

      beta = draw.beta(beta, xx, bhat, sig2, tau, u, w, alpha)

      sig2 = draw.sig2(beta, X, y)
      # sig2 = 2500

      if(i > burn)
      {
          output$u[i-burn,] <- u
          output$w[i-burn,] <- w
          output$beta[i-burn,] <- beta
          output$sig2[i-burn] = sig2
          output$tau[i-burn] = tau
      }
  }

  output
}

bridge.reg.know.sig2.R <- function(y, X,  nsamp, alpha=0.5, sig2=var(lm(y~X)$residuals),
                                   tau=1.0, burn=100, verbose=500)
{
  X <- as.matrix(X)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ixx <- chol2inv(chol(xx))

  p <- ncol(X)

  bhat <- drop(ixx%*%xy)
  beta <- rep(1,p)
  tau <- 2
  w = rep(2,p)

  #beta = bhat + rnorm(1,0,1)
  #w = 2*abs(beta)
  #tau=1

  # sig2 = (1/length(y))*sum((y-X%*%bhat)^2)

  output <- list(u = matrix(nrow=nsamp, ncol=length(beta)),
                 w = matrix(nrow=nsamp, ncol=length(beta)),
                 beta = matrix(nrow=nsamp, ncol=length(beta))
                 # sig2 = rep(0, nsamp),
                 # tau = rep(0, nsamp)
                 )

  for( i in 1:(nsamp+burn))
  {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")

      u <- draw.u(tau, beta, w, alpha)
      w <- draw.w(tau, beta, u, alpha)
      # tau = draw.tau(beta, alpha, c, d)

      beta = draw.beta(beta, xx, bhat, sig2, tau, u, w, alpha)
      
      # sig2 = draw.sig2(beta, X, y)
      

      if(i > burn)
      {
          output$u[i-burn,] <- u
          output$w[i-burn,] <- w
          output$beta[i-burn,] <- beta
          # output$sig2[i-burn] = sig2
          # output$tau[i-burn] = tau
      }
  }

  output
}

bridge.reg.know.tau.R <- function(y, X, nsamp, alpha=0.5, tau=1.0, burn=100, verbose=500)
{
  X <- as.matrix(X)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ixx <- chol2inv(chol(xx))

  p <- ncol(X)

  bhat <- drop(ixx%*%xy)
  beta <- rep(1,p)
  ## tau <- 2
  w = rep(2,p)
  u = 0

  #beta = bhat + rnorm(1,0,1)
  #w = 2*abs(beta)
  #tau=1

  sig2 = (1/length(y))*sum((y-X%*%bhat)^2)

  output <- list(u = matrix(nrow=nsamp, ncol=length(beta)),
                 w = matrix(nrow=nsamp, ncol=length(beta)),
                 beta = matrix(nrow=nsamp, ncol=length(beta)),
                 sig2 = rep(0, nsamp)
                 # tau = rep(0, nsamp)
                 )

  svdl = svd(X)
  U = svdl$u
  V = svdl$v
  d = svdl$d
  A = U %*% diag(d, p)
  tV = t(V);
  a = t(A) %*% y;

  ## u = 0
  ## lm1 = lm(y~X);
  ## sig2 = var(lm1$res)
  
  for( i in 1:(nsamp+burn))
  {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")

      w <- draw.w(tau, beta, u, alpha)
      sig2 = draw.sig2(beta, X, y)
      ## sig2 = 2500
      # tau = draw.tau(beta, alpha, c, d)
      # tau=50

      ## beta = draw.beta(beta, xx, bhat, sig2, tau, u, w, alpha)
      beta = draw.beta.2(beta, a, tV, d, sig2, tau, u, w, alpha)
      ## beta = draw.beta.3(beta, bhat, xx, sig2, tau, u, w, alpha)

      u <- draw.u(tau, beta, w, alpha)
      
      if(i > burn)
      {
          output$u[i-burn,] <- u
          output$w[i-burn,] <- w
          output$beta[i-burn,] <- beta
          output$sig2[i-burn] = sig2
          # output$tau[i-burn] = tau
      }
  }

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
