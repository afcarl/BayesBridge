library("coda")

dyn.unload("Bridge.so")
dyn.load("Bridge.so");

# Load the functions.
source("BridgeWrapper.R");

# Load the data.
## data(diabetes, package="lars");
load("diabetes.RData")
cov.name = colnames(diabetes$x);
y = diabetes$y;
X = diabetes$x;

y = y - mean(y);
mX = colMeans(X);
for(i in 1:442){ X[i,] = X[i,] - mX; }

LS = solve(t(X) %*% X, t(X) %*% y);

# Well, wrapping C in R doesn't necessarily speed things up that much.
# http://fluff.info/blog/arch/00000172.htm
# It can though... See Appendix!

nsamp = 10000

gb1 = bridge.reg.know.sig2(y, X, nsamp, 0.5, 2500.0, 100.0, 500);
gb3 = bridge.reg(y, X, nsamp, 0.5, 0.0, 0.0, 0.5, 0.5, 500);

gb2 = bridge.reg.know.tau(y, X, nsamp, 0.5, 100.0, 0.0, 0.0, 500, 1);
gb4 = bridge.reg.know.tau.stable(y, X, nsamp, 0.5, 100.0, 0.0, 0.0, 500);

for(i in 1:10){
    par(mfrow=c(1,2));
    hist(gb2$beta[i,], breaks=100, prob=TRUE);
    hist(gb4$beta[i,], breaks=100, prob=TRUE);
    ## acf(ans$beta[,i]);
    ## acf(gb3$beta[i,]);
    print(paste(mean(gb2$beta[i,]), mean(gb4$beta[i,])));
    print(paste(sd(gb2$beta[i,]), sd(gb4$beta[i,])));
    print(paste(effectiveSize(gb2$beta[i,]), effectiveSize(gb4$beta[i,])));
    readline("Press Enter...");
}

sstat = matrix(nrow=10, ncol=6);

sstat[,1] = apply(gb2$beta, 1, mean);
sstat[,2] = apply(gb4$beta, 1, mean);
sstat[,3] = apply(gb2$beta, 1, sd);
sstat[,4] = apply(gb4$beta, 1, sd);
sstat[,5] = apply(gb2$beta, 1, effectiveSize);
sstat[,6] = apply(gb4$beta, 1, effectiveSize);

# .C("donothing", y, X, 10.0, 0.5, 0.0, 0.0, 1000, 500);

# .C("passnothing");

# y, X, sig2, tau, alpha, lambda.max, tol, max.iter, use.cg
L = length(NuGrid);
P = dim(X)[2];
alpha = 0.5;
the.beta = array(0, dim=c(P, L));
for(i in 1:L){
    tau = NuGrid[i]^{-1/alpha};
    the.beta[,i] = bridge.EM(y, X, 1.0, tau, alpha, 10e8, 10e-9, 30, use.cg=FALSE);
}

gb = gb1;

bridge(X, y, 10000, 0.5, 0.5, 0.5, 500) -> mcmc.R

par(mfrow=c(1, 2));
for(i in 1:10){
    hist(gb$beta[i,], breaks=100, main=paste(cov.name[i], "um"), prob=TRUE);
    m.gbbeta = mean(gb$beta[i,]);
    s.gbbeta = sd(gb$beta[i,]);
    # hist(gb3$beta[i,], breaks=100, main=paste(cov.name[i], 3), prob=TRUE);
    hist(mcmc.R$beta[,i], breaks=100, main="R code", prob=TRUE);
    m.Rbeta = mean(mcmc.R$beta[,i]);
    s.Rbeta = sd(mcmc.R$beta[,i]);
    print(paste(m.gbbeta, s.gbbeta));
    print(paste(m.Rbeta, s.Rbeta));
    readline("Continue...");
}

bridge(X, y, 10000, 0.5, 0.5, 0.5, 500) -> mcmc.R

names(mcmc.R) = c("uuu", "www", "beta", "sig2", "tau");
beta = array(scan("beta.post"), dim=c(10, 10000));
uuu = array(scan("u.post"), dim=c(10, 10000));
www = array(scan("omega.post"), dim=c(10, 10000));

par(mfrow=c(1, 2));
for(i in 1:10){
    hist(beta[i,], breaks=100, main=paste(cov.name[i], "um"), prob=TRUE);
    m.gbbeta = mean(beta[i,]);
    s.gbbeta = sd(beta[i,]);
    # hist(gb3$beta[i,], breaks=100, main=paste(cov.name[i], 3), prob=TRUE);
    hist(mcmc.R$beta[,i], breaks=100, main="R code", prob=TRUE);
    m.Rbeta = mean(mcmc.R$beta[,i]);
    s.Rbeta = sd(mcmc.R$beta[,i]);
    print(paste(m.gbbeta, s.gbbeta));
    print(paste(m.Rbeta, s.Rbeta));
    readline("Continue...");
}

################################################################################
                         ## Expectation Maximization ##
################################################################################

GridSize = 10;
#NuGrid = c(seq(0.001,0.02*sum(abs(bhat)),length=GridSize/2), seq(0.02*sum(abs(bhat)),2*sum(abs(bhat)),length=GridSize/2))
NuGrid = sort(c(seq(-1,4.5,length=GridSize)))
NuGrid = rev(10^NuGrid)

data(diabetes, package="lars")
X = diabetes$x
Y = diabetes$y
Y = (Y - mean(Y))
p = ncol(X)
n = length(Y)
for(j in 1:p)
{
	X[,j] = (X[,j] - mean(X[,j]))
}

lse = lm(Y~X-1)
bhat = lse$coefficients

alpha = 0.5
tau = (NuGrid[4])^(-1/alpha);

Beta = bhat
Beta = bhat + rnorm(p,0,abs(bhat)/10)

diff = 1;

while(diff > 1e-9)
{
    YHat = X %*% Beta
    # sigma = sqrt(sum( (Y-YHat)^2 )/(n-p))
    sigma=1
    # EXPECTATION STEP
    LambdaInv = pmin( alpha*(tau^(2-alpha)) * abs(Beta)^(alpha-2), tau*1e7)
    #OmegaInv = as.numeric((d+1)/(d*sigma^2+(Y-YHat)^2))
    OmegaInv = rep(1, n)
    H = solve( (1/tau^2) * diag(as.numeric(LambdaInv)) + t(X) %*% diag(OmegaInv) %*% X) %*%
        t(X) %*% diag(OmegaInv)
    BetaNew = H %*% Y
    S = X %*% H
    diff = sum(abs(Beta - BetaNew))
    Beta = BetaNew
    # print(Beta);
    #Nu = (b.nu + sum(abs(Beta)/sigma))/(p + a.nu - 1)
}

######################################################################
                            ## APPENDIX ##
######################################################################

# When I run with num.samples = 1000, burn.in=500 I get a different
# timing ratio than when I run in with num.samples = 10000, burn.in =
# 500.  In turn I get a different timing ratio when using just C!

# Using R
# Burn-In Time Samples Time
#   500   6.67  10000  133.3
#   500   .658   2000  2.604
#   500   .271   1000  .5403

# Using C
# Burn-In Time Samples Time
#   500  .0106  10000  0.200

# To try and figure out what R is doing I created several different
# functions, each of which did the same thing, each of which used the
# data data passed to it in different ways.

# I found that when R knows that the data it passes is NOT MODIFIED,
# then the code will run (essentially) as fast as if I ran it from the
# command line directly.  When I do something that could be used to
# alter the R memory directly, then the code slows down A LOT.

# For instance, I thought I was being clever by creating MatrixFrames.
# One could then pass, for instance, the design matrix X directly into
# a MatrixFrame and go from there.  Clearly, this is not the case as R
# slows down considerably.  Presumably, this is a secruity issue.  If
# R didn't check what I was doing with the memory it had allocated for
# itself, then I might be able to overflow the memory and insert some
# malliscious code.
