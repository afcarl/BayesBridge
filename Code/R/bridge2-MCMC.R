library(msm)
library(mvtnorm)
## library(truncnorm)

PPM <- function(X)
{
  ## Psuedoinverse.
  p <- ncol(X)
  xx = t(X) %*% X;
  evd = eigen(xx);
  idc = abs(evd$values / evd$values[1]) < 1e-16
  inv.values = 1.0 / evd$values;
  inv.values[idc] = 0;
  pixx = evd$vectors %*% diag(inv.values, p) %*% t(evd$vectors);
  ppm = X %*% pixx %*% t(X);
  ppm
}

bridge <- function(y, x,niter,alpha,c=1/2,d=1/2,burn=100, verbose=500, sig2, tau=1.0)
{
  X <- as.matrix(x)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ixx <- chol2inv(chol(xx))

  p <- ncol(X)

  ppm = PPM(X);
  
  ## Print t(X) %*% X;
  ## print(xx);

  ##  if(mygrid == 0)
  ##  {
  ##  	mylm = lm(y~x-1)
  ##  	myr = range(c(mylm$coefficients - 2*summary(mylm)[[4]][,2], mylm$coefficients + 2*summary(mylm)[[4]][,2]))
  ##  	mygrid = seq(myr[1], myr[2], length=250)
  ##  }
  
  ##  w <- rgamma(p,1)
  ##  for ( i in 1:p )
  ##    if ( runif(1)>alpha )
  ##      w[i] <- rgamma(1,2)

  ## YOU MUST SET beta = 1 for the MCMC to start correctly.
  ## You also might want to do this in the event that XX is singular.
  beta <- rep(1,p)
  ## beta = bhat
  ## tau <- 2
  w = rep(2,p)

  #beta = bhat + rnorm(1,0,1)
  #w = 2*abs(beta)
  #tau=1

  bhat <- drop(ixx%*%xy)
  ## sig2 = mean((y-X%*%bhat)^2)
  sig2 = mean((y-ppm%*%y)^2)

  output <- list(u = matrix(nrow=niter, ncol=length(beta)),
                 w = matrix(nrow=niter, ncol=length(beta)),
                 beta = matrix(nrow=niter, ncol=length(beta)),
                 sig2 = rep(0, niter),
                 tau = rep(0, niter)
                 )
  #betadens = array(0, c(length(mygrid), length(beta), niter))
  # betadens = matrix(0, nrow=length(mygrid), ncol=length(beta))
  for( i in 1:(niter+burn))
  {
    if( i%%verbose==0 ) cat("iteration ", i, "\n")

    u <- draw.u(tau,beta,w,alpha)
    w <- draw.w(tau,beta,u,alpha)
    ## tau = draw.tau(beta,alpha,c,d)

    sig2 = draw.sig2(beta,X,y)
    ## sig2 = 2500 ## SET SET SET
    
    myout = draw.beta(beta,xx,xy,sig2,tau,u,w,alpha, bhat)
    beta <- myout$beta
    
	if(niter > burn)
	{
          output$u[i-burn,] <- u
          output$w[i-burn,] <- w
          output$beta[i-burn,] <- beta
          ##betadens[,,i-burn] = myout$dens
          output$sig2[i-burn] = sig2
          output$tau[i-burn] = tau
          ## if(length(which(myout$dens > 1e15)) == 0)
          ## {
          ## 	betadens = betadens + (1/niter)*myout$dens
          ## }
 	}
  }

  ## output$betadens = apply(betadens,c(1,2),mean)
  output;
}


draw.u <- function(tau,beta,w,alpha)
{
  m = 1-{abs(beta/tau)*w^(-1/alpha)}
  ## m = 1-{0.5*abs(beta/tau)*w^(-1/alpha)}
  runif(length(beta),max=m)
}

draw.w <- function(tau,beta,u,alpha)
{
  p = length(beta)
  a = (abs(beta/tau)/(1-u))^alpha
  ## a = (0.5*abs(beta/tau)/(1-u))^alpha
  w = rgamma(p,1,1)
  pr = alpha/(1+alpha*a)
  #p1 = 0.5*(1+alpha)
  #pr = p1/(1-p1*a)
  for ( i in 1:p )
    if ( runif(1) < pr[i] )
      w[i] = rgamma(1,2,1)
  w+a
}

draw.tau <- function(beta, alpha,c,d)
{
	p = length(beta)
	nu = rgamma(1, c + p/alpha, rate=d + sum(abs(beta)^alpha))
	tau = nu^(-1/alpha)
	return(tau);
}

draw.sig2 <- function(beta,x,y)
{
	n = length(y)
	rss = sum( (as.matrix(y)-x%*%as.matrix(beta))^2 )
	prec = rgamma(1, 2+n/2, 2+rss/2)
	return(1/prec)
}

draw.beta <- function(beta,xx,xy,sig2,tau,u,w,alpha,bhat)
{
  p = length(beta)
  b = (1-u)*{w^(1/alpha)}*tau
  ## b = (1-u)*2*{w^(1/alpha)}*tau

  ## ## GIBBS in beta
  ## for ( i in 1:p )
  ## {
  ##   ## m = bhat[i] - crossprod(xx[i,-i],beta[-i]-bhat[-i])/xx[i,i]
  ##   ## m = bhat[i] - xx[i,-i] %*% (beta[-i] - bhat[-i]) / xx[i,i];
  ##   m = (xy[i] - xx[i,-i] %*% beta[-i]) / xx[i,i];
  ##   v = sig2/xx[i,i]
    
  ##   beta[i] = rtnorm(1,m,sqrt(v),-b[i],b[i])
  ##   ## myind = which(dens[,i] < Inf)
  ##   #if(length(myind) < length(mygrid))
  ##   #{ cat("Density evaluates to infinity\n")
  ##   #	cat(i, m,sqrt(v),-b[i],b[i], "\n", sep=" ")
  ##   #}
  ## }

  ## ## Joint Draw.
  ## This appears to take too long.
  ## require("tmvtnorm")
  ## print(b);
  ## beta = rtmvnorm(1, bhat, sig2*solve(xx), -1*b, b);

  beta = tnorm.new.gibbs(beta, xx, sig2, b, bhat);
  
  list(beta=beta)
}


tnorm.new.gibbs <- function(beta, xx, sig2, b, bhat) {

  p = ncol(xx);
  
  evd = eigen(xx);
  A    = evd$vectors %*% diag(sqrt(evd$values)) %*% t(evd$vectors);
  AInv = evd$vectors %*% diag(1.0 / sqrt(evd$values)) %*% t(evd$vectors);
  m = A %*% bhat
  v = sig2;

  z = A %*% beta;
  
  for (i in 1:p) {
    left  = rep(0, p);
    right = rep(0, p);
    left.top = rep(0, p);
    right.top = rep(0, p);
    AInvZ = rep(0, p);
    for (j in 1:p) {
      AInvZ[j] = AInv[j,-i] %*% z[-i];
      ## cat("i,j: ", i, j, "\n");
      ## cat("A[j,i]", AInv[j,i], "b[j]", b[j], "AInvZ", AInvZ[j],
      ##     "diff", b[j]-AInvZ[j], ( b[j]-AInvZ[j]) / AInv[j,i], "\n");
      ## Problems the sign.  Division be small pos/neg number can sometimes be
      ## wrong.  There is numerical instability with that operation.  There also
      ## appears to be numerical instability in the intervals.  This is with
      ## Efron's diabetes data.  The columans are nearly colinear I think.     
      if (AInv[j,i]==0) {
        left[j] = -1e16;
        right[j] = 1e16;
      }
      else if (AInv[j,i] > 0) {
        ## cat("greater\n");
        left.top[j] = ( -1*b[j] - AInvZ[j] )
        right.top[j] = (   b[j] - AInvZ[j] ) 
        left[j]  = ( -1*b[j] - AInvZ[j] ) / abs(AInv[j,i]);
        right[j] = (    b[j] - AInvZ[j] ) / abs(AInv[j,i]);
      }
      else {
        ## cat("less\n");
        left.top[j] = -1 * (b[j] - AInvZ[j])
        right.top[j] = (b[j] + AInvZ[j] )
        left[j]  = -1 * (b[j] - AInvZ[j]) / abs(AInv[j,i]);
        right[j] = (b[j] + AInvZ[j] ) / abs(AInv[j,i]);
      }
    }
    lmax = max(left);
    rmin = min(right);
    z[i] = rtnorm(1, m[i], sqrt(v), lmax, rmin);
    ## cat("m[i]: ", m[i], "z[i]: ", z[i], "\n");
    ## cat("left: ", left, "\n");
    ## cat("right: ", right, "\n");
    ## cat("left.top:", left.top, "\n");
    ## cat("right.top:", right.top, "\n");
    ## cat("left.mult:", left * abs(AInv[,i]) + AInvZ, "\n");
    ## cat("right.mult:", right * abs(AInv[,i]) + AInvZ, "\n");
    ## cat("azji:", abs(AInv[i,]) * z[i], "\n");
    ## beta = AInv %*% z;
    ## cat("beta: ", beta, "\n");
    ## cat("b: ", b, "\n");
  }

  beta = AInv %*% z;

  beta
}

if (FALSE) {

  source("bridge2-MCMC.R")
  
  M = matrix(c(1,0.0, 0.0, 1), nrow=2);
  bhat = c(1.0, 0.0);
  beta = rep(0, 2);
  sig2 = 1.0;
  b = c(2, 1);
  MInv = solve(M)
  
  nsamp = 10000;
  output = list(beta=matrix(nrow=nsamp, ncol=2));

  for (i in 1:nsamp) {
    beta = tnorm.new.gibbs(beta, MInv, sig2, b, bhat);
    output$beta[i,] = beta
  }

  par(mfrow=c(1,2))
  hist(output$beta[,1])
  hist(output$beta[,2])

  apply(output$beta, 2, mean)
  apply(output$beta, 2, sd)

  
  
}
