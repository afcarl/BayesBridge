library(msm)
library(mvtnorm)
library(truncnorm)

bridge <- function(x,y,niter,alpha,c=1/2,d=1/2,mygrid,burn=100, verbose=500)
{
  X <- as.matrix(x)
  xx <- t(X)%*%X
  xy <- t(X)%*%y
  ixx <- chol2inv(chol(xx))

#  if(mygrid == 0)
#  {
#  	mylm = lm(y~x-1)
#  	myr = range(c(mylm$coefficients - 2*summary(mylm)[[4]][,2], mylm$coefficients + 2*summary(mylm)[[4]][,2]))
#  	mygrid = seq(myr[1], myr[2], length=250)
#  }

  p <- ncol(X)


#  w <- rgamma(p,1)
#  for ( i in 1:p )
#    if ( runif(1)>alpha )
#      w[i] <- rgamma(1,2)



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
                 tau = rep(0, niter),
                 betadens = matrix(0, nrow=length(mygrid), ncol=length(beta))
                 )
  #betadens = array(0, c(length(mygrid), length(beta), niter))
  betadens = matrix(0, nrow=length(mygrid), ncol=length(beta))
  for( i in 1:(niter+burn))
  {
    if( i%%verbose==0 ) cat("iteration ", i, "\n")

    u <- draw.u(tau,beta,w,alpha)
    w <- draw.w(tau,beta,u,alpha)
    tau = draw.tau(beta,alpha,c,d)
    # tau=50


    myout = draw.beta(beta,xx,bhat,sig2,tau,u,w,alpha, mygrid)
    beta <- myout$beta
    sig2 = draw.sig2(beta,X,y)
    sig2 = 2500

	if(niter > burn)
	{
    	output$u[i-burn,] <- u
    	output$w[i-burn,] <- w
    	output$beta[i-burn,] <- beta
 		#betadens[,,i-burn] = myout$dens
 		output$sig2[i-burn] = sig2
 		output$tau[i-burn] = tau
 		if(length(which(myout$dens > 1e15)) == 0)
 		{
 			betadens = betadens + (1/niter)*myout$dens
 		}
 	}
  }

  output$betadens = apply(betadens,c(1,2),mean)
  output
}


draw.u <- function(tau,beta,w,alpha)
{
  m = 1-{abs(beta/tau)*w^(-1/alpha)}
  runif(length(beta),max=m)
}

draw.w <- function(tau,beta,u,alpha)
{
  p = length(beta)
  a = (abs(beta/tau)/(1-u))^alpha
  w = rgamma(p,1,1)
  pr = alpha/(1+alpha*a)
  #p1 = 0.5*(1+alpha)
  #pr = p1/(1-p1*a)
  for ( i in 1:p )
    if ( runif(1) < pr[i] )
      w[i] = rgamma(1,2,1)
  w+a
}

draw.tau = function(beta, alpha,c,d)
{
	p = length(beta)
	nu = rgamma(1, c + p/alpha, rate=d + sum(abs(beta)^alpha))
	tau = nu^(-1/alpha)
	return(tau);
}

draw.sig2 = function(beta,x,y)
{
	n = length(y)
	rss = sum( (as.matrix(y)-x%*%as.matrix(beta))^2 )
	prec = rgamma(1, 2+n/2, 2+rss/2)
	return(1/prec)
}

draw.beta <- function(beta,xx,bhat,sig2,tau,u,w,alpha,mygrid)
{
  p = length(bhat)
  b = (1-u)*{w^(1/alpha)}*tau
  dens = matrix(0, nrow=length(mygrid), ncol=length(beta))

  for ( i in 1:p )
  {
    ## Factoring joint.
    m = bhat[i] - crossprod(xx[i,-i],beta[-i]-bhat[-i])/xx[i,i]
    v = sig2/xx[i,i]
    beta[i] = rtnorm(1,m,sqrt(v),-b[i],b[i])
    dens[,i] = dtnorm(mygrid,m,sqrt(v),-b[i],b[i])
    myind = which(dens[,i] < Inf)
    #if(length(myind) < length(mygrid))
    #{ cat("Density evaluates to infinity\n")
    #	cat(i, m,sqrt(v),-b[i],b[i], "\n", sep=" ")
    #}

  }
  list(beta=beta,dens=dens)
}
