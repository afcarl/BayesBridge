## Load C++ code.
## dyn.load("libbridge.so");
## dyn.load("Bridge.so");
## dyn.load("Bridge.so", PACKAGE=BayesBridge);


################################################################################
                             ## HELPER FUNCTIONS ##
################################################################################

# Check if param >= val.  "name" is the name of the param.
is.above <- function(param, val, name){
    above = TRUE;
    if (param < val){
        alert = paste("Error: ", name, "<", val, sep="");
        print(alert);
        above = FALSE;
    }
    # While we are at it, check that we are working with a number.
    if (!is.numeric(param)) {
        alert = paste("Error:", name, "is not numeric.");
        print(alert);
        above = FALSE;
    }
    above;
}

# Check that the parameters are valid.
check.parameters <- function(N, R, M, sig2, tau, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate){
    ok = TRUE;
    if (N!=R)    { print("Error: y and X do not conform."); ok=FALSE; }
    ok = ok *
         is.above(M         , 1, "niter") *
         is.above(sig2      , 0, "sig2") *
         is.above(tau       , 0, "tau") *
         is.above(alpha     , 0, "alpha") *
         is.above(sig2.shape, 0, "sig2.shape") *
         is.above(sig2.scale, 0, "sig2.scale") *
         is.above(nu.shape  , 0, "nu.shape") *
         is.above(nu.rate   , 0, "nu.rate");

    ## if (M > 1000){
    ##     ans = readline(paste("niter =", M, "> 1000.  Do you really want to proceed? [n] "));
    ##     ans = substr(ans, 1, 1);
    ##     ok = ok * (ans=="y" || ans=="Y");
    ## }

    ok
}

# Check the extra parameters in expectation maximization.
check.EM <- function(lambda.max, tol, max.iter)
{
    ok = TRUE;
    ok = ok *
        is.above(lambda.max, 0.0, "lambda.max") *
        is.above(tol       , 0.0, "tolerance") *
        is.above(max.iter  , 1.0, "max.iter");

    ok
}

################################################################################
                         ## EXPECTATION MAXIMIZATION ##
################################################################################

bridge.EM <- function(y, X,
                      alpha=0.5,
                      ratio=1.0,
                      lambda.max=1e9*ratio, tol=1e-9, max.iter=30,
                      use.cg=FALSE, ret.solves=FALSE)
{
    N = length(y);
    R = dim(X)[1];
    P = dim(X)[2];

    sig2 = 1.0;
    tau = ratio;

    ok = check.parameters(N, R, 1, sig2, tau, alpha, 1.0, 1.0, 1.0, 1.0) *
         check.EM(lambda.max, tol, max.iter);

    if (!ok) { break; }

    beta = array(0, dim=c(P));

    OUT = .C("bridge_EM",
             beta,
             as.double(y), as.double(X),
             ratio, alpha,
             as.integer(P), as.integer(N),
             lambda.max, tol, as.integer(max.iter),
             as.integer(use.cg),
             PACKAGE="Bridge");

    output = OUT[[1]]; # beta
    rownames(output) = colnames(X);

    if(ret.solves) {
        output = list("beta"=output, "num.solves"=OUT[[10]]);
    }

    output
}

################################################################################
                            ## BRIDGE REGRESSION ##
################################################################################

bridge.reg.know.sig2 <- function(y, X,
                                 nsamp,
                                 alpha=0.5,
                                 sig2=var(lm(y~X)$residuals),
                                 tau=1.0,
                                 burn=500){
    N = length(y);
    R = dim(X)[1];
    P = dim(X)[2];
    M = nsamp;
    rt = 0;

    ok = check.parameters(N, R, M, sig2, tau, alpha, 1.0, 1.0, 1.0, 1.0);
    if (!ok) { break; }

    beta  = array(0, dim=c(P, M));
    u     = array(0, dim=c(P, M));
    omega = array(0, dim=c(P, M));

    OUT = .C("bridge_reg_know_sig2",
             beta, u, omega,
             as.double(y), as.double(X),
             sig2, tau, alpha,
             as.integer(P), as.integer(N), as.integer(M), as.integer(burn), rt,
             PACKAGE="Bridge");

    beta = OUT[[1]];
    rownames(beta) = colnames(X);

    output = list("beta"=beta, "u"=OUT[[2]], "w"=OUT[[3]], "runtime"=OUT[[13]]);
}

bridge.reg.know.tau <- function(y, X,
                                nsamp,
                                alpha=0.5,
                                tau=1.0,
                                sig2.shape=0.0, sig2.scale=0.0,
                                burn=500, beta.iter=1){
    N = length(y);
    R = dim(X)[1];
    P = dim(X)[2];
    M = nsamp;
    rt = 0;

    ok = check.parameters(N, R, M, 1.0, tau, alpha, sig2.shape, sig2.scale, 1.0, 1.0);
    if (!ok) { break; }

    beta  = array(0, dim=c(P, M));
    u     = array(0, dim=c(P, M));
    omega = array(0, dim=c(P, M));
    sig2  = array(0, dim=c(M));

    OUT = .C("bridge_reg_know_tau",
             beta, u, omega, sig2,
             as.double(y), as.double(X),
             tau, alpha,
             sig2.shape, sig2.scale,
             as.integer(P), as.integer(N), as.integer(M), as.integer(burn), rt, as.integer(beta.iter),
             PACKAGE="Bridge");

    output = list("beta"=OUT[[1]], "u"=OUT[[2]], "w"=OUT[[3]], "sig2"=OUT[[4]], "runtime"=OUT[[15]])
    rownames(output$beta) = colnames(X);

    output
}

bridge.reg <- function(y, X,
                       nsamp,
                       alpha=0.5,
                       sig2.shape=0.0, sig2.scale=0.0,
                       nu.shape=0.5, nu.rate=0.5,
                       burn=500){
    N = length(y);
    R = dim(X)[1];
    P = dim(X)[2];
    M = nsamp;
    rt = 0;

    ok = check.parameters(N, R, M, 1.0, 1.0, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate);
    if (!ok) { break; }

    beta  = array(0, dim=c(P, M));
    u     = array(0, dim=c(P, M));
    omega = array(0, dim=c(P, M));
    sig2  = array(0, dim=c(M));
    tau   = array(0, dim=c(M));

    OUT <- .C("bridge_regression_all",
              beta, u, omega, sig2, tau,
              as.double(y), as.double(X),
              alpha,
              sig2.shape, sig2.scale,
              nu.shape, nu.rate,
              true.sig2=0.0, true.tau=0.0,
              as.integer(P), as.integer(N), as.integer(M), as.integer(burn), rt,
              PACKAGE="Bridge");

    output <- list("beta"=OUT[[1]], "u"=OUT[[2]], "w"=OUT[[3]], "sig2"=OUT[[4]], "tau"=OUT[[5]],
                   "runtime"=OUT[[19]]);

    rownames(output$beta) = colnames(X);

    output
}

## Draw truncated normal
##------------------------------------------------------------------------------
rtnorm.left <- function(num=1, left=0.0, mu=0.0, sig=1.0)
{
    ## Check Parameters.
    if (sum(sig<=0)!=0) {
        print("sig must be greater than zero.");
        return(NA);
    }
    if (! (num>0) ) {
      print("num must be greater than zero.");
      return(NA);
    }

    x = rep(0, num);
    
    if (length(mu)  != num) { mu  = array(mu,  num); }
    if (length(sig) != num) { sig = array(sig, num); }

    if (length(left)  != num) { left  = array(left,  num); }
    
    OUT = .C("rtnorm_left", x, left, mu, sig, as.integer(num), PACKAGE="Bridge");

    OUT[[1]]
}

rtnorm.both <- function(num=1, left=-1.0, right=1.0, mu=0.0, sig=1.0)
{
  LGER = left>=right;
    ## Check Parameters.
    if (sum(sig<=0)!=0) {
        print("sig must be greater than zero.");
        return(NA);
    }
    if (sum(LGER)!=0) {
      print("left must be less than right.");
      return(NA);
    }
    if (! (num>0) ) {
      print("num must be greater than zero.");
      return(NA);
    }

    x = rep(0, num);

    if (length(mu)  != num) { mu  = array(mu,  num); }
    if (length(sig) != num) { sig = array(sig, num); }

    if (length(left)  != num) { left  = array(left,  num); }
    if (length(right) != num) { right = array(right, num); }

    OUT = .C("rtnorm_both", x, left, right, mu, sig, as.integer(num), PACKAGE="Bridge");

    OUT[[1]]
}

rtnorm <- function(num=1, mu=0.0, sig=1.0, left=-Inf, right=Inf)
{
  LGER = left>=right;
  ## Check Parameters.
  if (sum(sig<=0)!=0) {
    print("sig must be greater than zero.");
    return(NA);
  }
  if (sum(LGER)!=0) {
    print("left must be less than right.");
    cat("left:", left[LGER], "\n");
    cat("rght:", right[LGER], "\n");
    return(NA);
  }

  if (length(mu)  != num) { mu  = array(mu,  num); }
  if (length(sig) != num) { sig = array(sig, num); }
  
  if (length(left)  != num) { left  = array(left,  num); }
  if (length(right) != num) { right = array(right, num); }
  
  x = rep(0, num);
  
  u  = (left ==-Inf) & (right == Inf);
  l  = (left !=-Inf) & (right == Inf);
  r  = (left ==-Inf) & (right != Inf);
  b  = (left !=-Inf) & (right != Inf);

  n.u = sum(u);
  n.l = sum(l);
  n.r = sum(r);
  n.b = sum(b);

  if (n.b > 0) {
    x[b] = .C("rtnorm_both", x[b], left[b], right[b], mu[b], sig[b], as.integer(n.b), PACKAGE="Bridge")[[1]];
  }
  if (n.r > 0) {
    x[r] = -1*.C("rtnorm_left", x[r], -1*right[r], -1*mu[r], sig[r], as.integer(n.r), PACKAGE="Bridge")[[1]];
  }
  if (n.l > 0) {
    x[l] = .C("rtnorm_left", x[l], left[l], mu[l], sig[l], as.integer(n.l), PACKAGE="Bridge")[[1]]
  }
  if (n.u > 0) {
    x[u] = rnorm(n.u, mu[u], sig[u]);
  }

  x
}

##------------------------------------------------------------------------------

bridge.reg.know.tau.stable <- function(y, X,
                                       nsamp,
                                       alpha=0.5,
                                       tau=1.0,
                                       sig2.shape=0.0, sig2.scale=0.0,
                                       burn=500){
    N = length(y);
    R = dim(X)[1];
    P = dim(X)[2];
    M = nsamp;
    rt = 0;
    
    ok = check.parameters(N, R, M, 1.0, tau, alpha, sig2.shape, sig2.scale, 1.0, 1.0);
    if (!ok) { break; }

    beta   = array(0, dim=c(P, M));
    lambda = array(0, dim=c(P, M));
    sig2  = array(0, dim=c(M));

    OUT = .C("bridge_reg_know_tau_stable",
             beta, lambda, sig2,
             as.double(y), as.double(X),
             tau, alpha,
             sig2.shape, sig2.scale,
             as.integer(P), as.integer(N), as.integer(M), as.integer(burn), rt,
             PACKAGE="Bridge");

    output = list("beta"=OUT[[1]], "lambda"=OUT[[2]], "sig2"=OUT[[3]], "runtime"=OUT[[14]])
    rownames(output$beta) = colnames(X);

    output
}

bridge.reg.stable <- function(y, X,
                              nsamp,
                              alpha=0.5,
                              sig2.shape=0.0, sig2.scale=0.0,
                              nu.shape=2.0, nu.rate=1/2,
                              burn=500){
    N = length(y);
    R = dim(X)[1];
    P = dim(X)[2];
    M = nsamp;
    rt = 0;
    
    ok = check.parameters(N, R, M, 1.0, 1.0, alpha, sig2.shape, sig2.scale, nu.shape, nu.rate);
    if (!ok) { break; }

    beta   = array(0, dim=c(P, M));
    lambda = array(0, dim=c(P, M));
    sig2  = array(0, dim=c(M));
    tau   = array(0, dim=c(M));

    OUT <- .C("bridge_reg_stable",
              beta, lambda, sig2, tau,
              as.double(y), as.double(X),
              alpha,
              sig2.shape, sig2.scale,
              nu.shape, nu.rate,
              as.integer(P), as.integer(N), as.integer(M), as.integer(burn), rt,
              PACKAGE="Bridge");

    output = list("beta"=OUT[[1]], "lambda"=OUT[[2]], "sig2"=OUT[[3]], "tau"=OUT[[4]], "runtime"=OUT[[16]])
    rownames(output$beta) = colnames(X);

    output
}
