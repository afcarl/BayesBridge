# Load the data.
data(diabetes, package="lars");
cov.name = colnames(diabetes$x);
y = diabetes$y;
X = diabetes$x;
id = colnames(X);

# Center.
y = y - mean(y);
mX = colMeans(X);
for(i in 1:442){ X[i,] = X[i,] - mX; }

# The least squares solution.
LS = solve(t(X) %*% X, t(X) %*% y);

# Set parameters.
alpha = 0.5;
sig2  = 2500;
tau   = sqrt(sig2)*1e-5;
N = 2000;

# Expectation maximization
beta.EM = bridge.EM(y, X, alpha, tau/sqrt(sig2), 10e8, 10e-9, 30, use.cg=FALSE);

# Bridge Regression using R routine.  This will take a while.
readline("");
g.R = bridge.reg.R(y, X, N, alpha, 0.5, 0.5, 500);

# Bridge Regression using C routine.  Same number of samples.
readline("");
g.C = bridge.reg(y, X, N, alpha, 0.0, 0.0, 0.5, 0.5, 500);

readline("Histograms: Press <Return> to continue.");

par(mfrow=c(1,2));
par(mai=c(1.02, 0.82, 0.82, 0.42));
for(i in 1:10){
    # Summary statistics.
    print("Mean and SD using R and C.");
    print(paste(mean(g.R$beta[,i]), mean(g.C$beta[i,])));
    print(paste(sd(g.R$beta[,i]), sd(g.C$beta[i,])));
    # Marginal Posterior.
    hist(g.R$beta[,i], breaks=100, prob=TRUE, ylab="dens", xlab=id[i],
         main=paste(id[i], "using R"));
    hist(g.C$beta[i,], breaks=100, prob=TRUE, ylab="dens", xlab=id[i],
         main=paste(id[i], "using C"));
    # Pause.
    # readline("Press <Return> to continue.");
}

# Now do this with lots of samples.
N = 50000;
# Bridge Regression using C routine.  Lots of samples
readline("");
g.C = bridge.reg(y, X, N, alpha, 0.0, 0.0, 0.5, 0.5, 500);

readline("Histograms: Press <Return> to continue.");
par(mfrow=c(5,2));
par(mai=c(0.4,0.4,0.1,0.1));
for(i in 1:10){
    # Marginal Posterior.
    hist(g.C$beta[i,], breaks=100, prob=TRUE, ylab="dens", xlab=id[i],
         main=id[i]);
    # Joint Modes.
    abline(v=LS[i], col=2);
    abline(v=beta.EM[i], col=3);
    legend("topright", legend=c("LS", "EM"), col=c(2,3), lty=c(1,1));
}
