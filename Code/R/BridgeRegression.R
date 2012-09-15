# library(msm)
# library(mvtnorm)
library(truncnorm)

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
	prec = rgamma(1, 2+n/2, 2+rss/2)
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

bridge.reg.R <- function(y, x, niter, alpha, c=1/2, d=1/2, burn=100, verbose=500)
{
  X <- as.matrix(x)
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

  sig2 = (1/length(y))*sum((y-X%*%bhat)^2)

  output <- list(u = matrix(nrow=niter, ncol=length(beta)),
                 w = matrix(nrow=niter, ncol=length(beta)),
                 beta = matrix(nrow=niter, ncol=length(beta)),
                 sig2 = rep(0, niter),
                 tau = rep(0, niter)
                 )

  for( i in 1:(niter+burn))
  {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")

      u <- draw.u(tau, beta, w, alpha)
      w <- draw.w(tau, beta, u, alpha)
      tau = draw.tau(beta, alpha, c, d)
      # tau=50

      beta = draw.beta(beta, xx, bhat, sig2, tau, u, w, alpha)

      sig2 = draw.sig2(beta, X, y)
      # sig2 = 2500

      if(niter > burn)
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

bridge.reg.know.sig2.R <- function(y, x, niter, alpha, sig2, tau, burn=100, verbose=500)
{
  X <- as.matrix(x)
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

  output <- list(u = matrix(nrow=niter, ncol=length(beta)),
                 w = matrix(nrow=niter, ncol=length(beta)),
                 beta = matrix(nrow=niter, ncol=length(beta))
                 # sig2 = rep(0, niter),
                 # tau = rep(0, niter)
                 )

  for( i in 1:(niter+burn))
  {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")

      u <- draw.u(tau, beta, w, alpha)
      w <- draw.w(tau, beta, u, alpha)
      # tau = draw.tau(beta, alpha, c, d)

      beta = draw.beta(beta, xx, bhat, sig2, tau, u, w, alpha)

      # sig2 = draw.sig2(beta, X, y)

      if(niter > burn)
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

bridge.reg.know.tau.R <- function(y, x, niter, alpha, tau, burn=100, verbose=500)
{
  X <- as.matrix(x)
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

  sig2 = (1/length(y))*sum((y-X%*%bhat)^2)

  output <- list(u = matrix(nrow=niter, ncol=length(beta)),
                 w = matrix(nrow=niter, ncol=length(beta)),
                 beta = matrix(nrow=niter, ncol=length(beta)),
                 sig2 = rep(0, niter)
                 # tau = rep(0, niter)
                 )

  for( i in 1:(niter+burn))
  {
      if( i%%verbose==0 ) cat("iteration ", i, "\n")

      u <- draw.u(tau, beta, w, alpha)
      w <- draw.w(tau, beta, u, alpha)
      # tau = draw.tau(beta, alpha, c, d)
      # tau=50

      beta = draw.beta(beta, xx, bhat, sig2, tau, u, w, alpha)

      sig2 = draw.sig2(beta, X, y)
      # sig2 = 2500

      if(niter > burn)
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
