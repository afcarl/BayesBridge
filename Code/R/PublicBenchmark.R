## Nicholas Poslon, James Scott, and Jesse Windle, 12-10-27.

## This script benchmarks the mixture of triangles method against the mixture of
## normals method.  The function compare.it returns an object info, which is a
## list with the following information:

## tri.info - summary statistics for mixture of triangles
## stb.info - summary statistics for mixture of normals
## tri.list - list of statistics for each run.
## stb.list - list of statistics for each run.
## tri.stat - array of statistics for each run -- good for using apply.
## stb.stat - array of statistics for each run -- good for using apply.
## tri.gb   - a single MC from the run of simluations.
## stb.gb   - a single MC from the run of simulations.
## n        - the number of response observations.

library("coda")

################################################################################
                                  ## Setup ##
################################################################################

## Set to TRUE to run benchmark.

run <- list("EFRON" = TRUE, # Efron's diabetes data
            "DBI"   = TRUE, # Efron's diabetes data with interactions and squared terms.
            "BH"    = TRUE, # Boston Housing
            "BHI"   = TRUE) # Boston Housing with interactions and squared terms.

## Ortogonalizing matrices by QR.
oth <- list("EFRON" = TRUE,
            "DBI"   = TRUE,
            "BH"    = TRUE,
            "BHI"   = TRUE)

## RUN INFO
nsamp = 10000   ## Number of samples
burn = 10000     ## Burn-in.
alpha = 0.5
ntrials = 2     ## Number of simulations.
tau = 0          ## Set to <= 0 for unknown tau.
betaburn = 0     ## You can burn the beta sampler.
## use.hmc = FALSE ## DEPRECATED ## Set to TRUE if you want to use HMC.
extras = FALSE   ## To take place of use.hmc for testing purposes.
inflate = 1      ## Just leave at 1 for now.

save.it  = FALSE  ## Write output to file.
plot.it  = FALSE  ## Plot histograms.
print.it = TRUE   ## Print summary.

################################################################################
                           ## Load .so or package ##
################################################################################

library("BayesBridge");

################################################################################
                     ## Make sure everything is working ##
################################################################################

if (FALSE) {

  data("diabetes", package="BayesBridge");
  cov.name = colnames(diabetes$x);
  y = diabetes$y;
  X = diabetes$x;

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }

  LS = solve(t(X) %*% X, t(X) %*% y);

  nsamp = 100000
  burn  = 1000

  gb1 = bridge.reg.tri(y, X, nsamp, 0.5, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, sig2.true=0, tau.true=0, burn=burn, ortho=FALSE, betaburn=betaburn);
  gb2 = bridge.reg.stb(y, X, nsamp, 0.5, 0.0, 0.0, 2.0, 2.0, 1.0, 1.0, sig2.true=0, tau.true=0, burn=burn, ortho=FALSE);

  sstat.1 = sum.stat(gb1);
  sstat.2 = sum.stat(gb2);

  sstat.1
  sstat.2

  ## Make sure things look right.
  for(i in 1:10){
    par(mfrow=c(1,2));
    hist(gb1$beta[,i], breaks=100, prob=TRUE);
    ## hist(gb2$beta[,i], breaks=100, prob=TRUE);
    hist(gb2$beta[,i], breaks=100, prob=TRUE);
    readline("Press Enter...");
  }

}

################################################################################

################################################################################

## For printing out header.
print.header <- function(the.title, char="#", n=70)
{
    bar = paste(rep(char, n), collapse="")
    cat("\n", bar, "\n", sep="")
    cat(char, the.title, "\n")
    cat(bar, "\n\n", sep="")
}

## Generates summary statistics.
sum.stat <- function(gbs, thin=1)
{
  n = nrow(gbs$beta);
  p = ncol(gbs$beta);
  gbs$beta = gbs$beta[seq(1, n, thin),];
  if (p==1) {gbs$beta = as.matrix(gbs$beta);}
  sstat = matrix(nrow=p, ncol=7);

  sstat[,1] = apply(gbs$beta, 2, mean);
  sstat[,2] = apply(gbs$beta, 2, sd);
  sstat[,3] = apply(gbs$beta, 2, effectiveSize);
  sstat[,4] = sstat[,3] / gbs$runtime;
  ## sstat[,5] = apply(gbs$beta, 2, ESS);
  ## sstat[,6] = sstat[,5] / gbs$runtime;
  sstat[,5] = sstat[,1] / sstat[,2];
  sstat[,6] = apply(gbs$beta, 2, function(x){quantile(x, 0.1)});
  sstat[,7] = apply(gbs$beta, 2, function(x){quantile(x, 0.9)});

  ## colnames(sstat) = c("mean", "sd", "ESS", "ESS.sec", "myESS", "myESS.sec", "t", "Q.1", "Q.9");
  colnames(sstat) = c("mean", "sd", "ESS", "ESS.sec", "t", "Q.1", "Q.9");

  sstat
}

################################################################################
                             ## Compare Routines ##
################################################################################

compare.it <- function(y, X, nsamp=10000,
                       alpha=0.5,
                       sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,
                       tau=0.0, 
                       burn=0, ntrials=1, ortho=FALSE, betaburn=0, extras=FALSE, inflate=1.0)
{
  ## Returns a list of summary statistics, the last Gibbs MCMC, and a
  ## table with the average statistics over several runs.
  
  tri.info = matrix(nrow=ntrials, ncol=5);
  colnames(tri.info) = c("ave.ESS", "sd.ESS", "ave.ESS.sec", "sd.ESS.sec", "runtime");
  stb.info = tri.info;

  tri.stat = list();
  stb.stat = list();

  for (i in 1:ntrials){
    
    gb.tri = bridge.reg.tri(y, X, round(nsamp*inflate), alpha=alpha,
      sig2.shape=sig2.shape, sig2.scale=sig2.scale, nu.shape=nu.shape, nu.rate=nu.rate,
      sig2=0.0, tau=tau, burn=round(burn*inflate), ortho=ortho, betaburn=betaburn, extras=extras);
    gb.stb = bridge.reg.stb(y, X, nsamp, alpha=alpha,
      sig2.shape=sig2.shape, sig2.scale=sig2.scale, nu.shape=nu.shape, nu.rate=nu.rate,
      sig2=0.0, tau=tau, burn=burn, ortho=ortho)
    
    sstat.tri = sum.stat(gb.tri);
    sstat.stb = sum.stat(gb.stb);

    tri.stat[[i]] = sstat.tri;
    stb.stat[[i]] = sstat.stb;

    tri.info[i,"ave.ESS"] = mean(sstat.tri[,"ESS"]);
    tri.info[i,"sd.ESS" ] = sd(sstat.tri[,"ESS"]);
    tri.info[i,"ave.ESS.sec"] = mean(sstat.tri[,"ESS.sec"]);
    tri.info[i,"sd.ESS.sec"]  = sd(sstat.tri[,"ESS.sec"]);
    tri.info[i,"runtime"] = gb.tri$runtime

    stb.info[i,"ave.ESS"] = mean(sstat.stb[,"ESS"]);
    stb.info[i,"sd.ESS" ] = sd(sstat.stb[,"ESS"]);
    stb.info[i,"ave.ESS.sec"] = mean(sstat.stb[,"ESS.sec"]);
    stb.info[i,"sd.ESS.sec"]  = sd(sstat.stb[,"ESS.sec"]);
    stb.info[i,"runtime"] = gb.stb$runtime
  }

  OUT <- list("tri.info"=tri.info, "stb.info"=stb.info,
              "tri.list"=tri.stat, "stb.list"=stb.stat,
              "tri.stat"=simplify2array(tri.stat), "stb.stat"=simplify2array(stb.stat),
              "tri.gb"=gb.tri, "stb.gb"=gb.stb);

  OUT
}

##------------------------------------------------------------------------------

plot.info <- function(info, cnames, P=ncol(info$tri.gb$beta))
{
  ## sig2
  par(mfrow=c(1,2));
  hist(info$tri.gb$sig2, breaks=100, prob=TRUE, main=paste("tri: sig2"));
  hist(info$stb.gb$sig2, breaks=100, prob=TRUE, main=paste("stb: sig2"));
  cat("sig2: \n");
  cat("tri:", mean(info$tri.gb$sig2), sd(info$tri.gb$sig2), "\n");
  cat("stb:", mean(info$stb.gb$sig2), sd(info$stb.gb$sig2), "\n");
  readline("Press Enter...");

  ## tau
  if (!is.null(info$tri.gb$tau)) {
    par(mfrow=c(1,2));
    hist(info$tri.gb$tau, breaks=100, prob=TRUE, main=paste("tri: tau"));
    hist(info$stb.gb$tau, breaks=100, prob=TRUE, main=paste("stb: tau"));
    cat("tau: \n");
    cat("tri:", mean(info$tri.gb$tau), sd(info$tri.gb$tau), "\n");
    cat("stb:", mean(info$stb.gb$tau), sd(info$stb.gb$tau), "\n");
    readline("Press Enter...");
  }

  ## beta
  for(i in 1:P){
    par(mfrow=c(1,2));
    hist(info$tri.gb$beta[,i], breaks=100, prob=TRUE, main=paste("tri:", cnames[i]));
    hist(info$stb.gb$beta[,i], breaks=100, prob=TRUE, main=paste("stb:", cnames[i]));
    cat("Variable", i, "\n");
    cat("tri:", mean(info$tri.gb$beta[,i]), sd(info$tri.gb$beta[,i]), "\n");
    cat("stb:", mean(info$stb.gb$beta[,i]), sd(info$stb.gb$beta[,i]), "\n");
    readline("Press Enter...");
  }

}

hist.info <- function(info, cnames, P=ncol(info$tri.gb$beta), breaks=40)
{
  par(mfrow=c(1,1));
  for (i in 1:P) {
    beta.1 = info$tri.gb$beta[info$tri.gb$shape[,i]==1,i]
    beta.2 = info$tri.gb$beta[info$tri.gb$shape[,i]==2,i]
    h1 = hist(beta.1, breaks=breaks)
    h2 = hist(beta.2, breaks=breaks)
    ymax=max(c(h1$counts, h2$counts));
    xmin = min(c(h1$breaks, h2$breaks));
    xmax = max(c(h1$breaks, h2$breaks))
    plot(h1, col="#FF000088", ylim=c(0,ymax), main=cnames[i])
    plot(h2, col="#0000FF66", add=TRUE)
    readline("<ENTER>");
  }
}

##------------------------------------------------------------------------------

run.it <- function(y, X, nsamp=1000,  burn=100, 
                   alpha=0.5, sig2.shape=0.0, sig2.scale=0.0, nu.shape=2.0, nu.rate=2.0,
                   ntrials=1, tau=0.0,
                   save.it=FALSE, print.it=FALSE, plot.it=FALSE, name="somerun", ortho=FALSE, betaburn=0, extras=FALSE, inflate=1.0)
{
  ## Runs the comparison and plots/prints/saves the resulting data.

  info = compare.it(y, X, nsamp=nsamp,
    alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, nu.shape, nu.rate, tau=tau,
    burn=burn, ntrials=ntrials, ortho=ortho, betaburn=betaburn, extras=extras, inflate=inflate)

  ## Make sure things look right.
  if (plot.it) {
    cnames = colnames(X);
    plot.info(info, cnames);
  }

  if (print.it) {
    print("tri.stat:")
    print(info$tri.list[[ntrials]])
    print("stb.stat:")
    print(info$stb.list[[ntrials]])
  }

  cat("Info size:", object.size(info) / 2^20, "Mb\n")
  print("tri:")
  print(info$tri.info)
  print("stb:")
  print(info$stb.info)

  info$n = length(y);

  filename = paste("info", name, "RData", sep=".");
  if (save.it) save(info, file=filename, compress=TRUE);

  info
}

##------------------------------------------------------------------------------

table.info <- function(info, colnum=3)
{
  n = info$n
  p = ncol(info$tri.gb$beta)
  
  tab = matrix(nrow=2, ncol=7);
  colnames(tab) = c("n", "p", "Time", "Min", "Med", "Max", "SD");

  tab[,"n"]=n
  tab[,"p"]=p
  
  tri.med = apply(info$tri.stat[,colnum,], 1, median);
  stb.med = apply(info$stb.stat[,colnum,], 1, median);

  tab[1,"Min"] = min(tri.med);
  tab[1,"Max"] = max(tri.med);
  tab[1,"Med"] = median(tri.med);
  tab[1,"SD" ] = sd(tri.med);
  tab[1,"Time"] = median(info$tri.info[,"runtime"]);

  tab[2,"Min"] = min(stb.med);
  tab[2,"Max"] = max(stb.med);
  tab[2,"Med"] = median(stb.med);
  tab[2,"SD" ] = sd(stb.med);
  tab[2,"Time"] = median(info$stb.info[,"runtime"]);

  tab
}

################################################################################
                              ## SUITE OF TESTS ##
################################################################################

##------------------------------------------------------------------------------
## Efron's diabetes data.

if (run$EFRON) {

  print.header("Efron's diabetes data.")

  if (tau>0) tau = 41
  
  ## Load data.
  data("diabetes", package="BayesBridge");
  cov.name = colnames(diabetes$x);
  y = diabetes$y;
  X = diabetes$x;
  n = length(y);
  cnames = cov.name
  
  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:n){ X[i,] = X[i,] - mX; }

  info <- run.it(y, X, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=print.it, plot.it=plot.it, name="efron",
                 ortho=FALSE, betaburn=betaburn, extras=extras, inflate=inflate)

}

##------------------------------------------------------------------------------
## Efron's diabetes data - with interactions and squared terms.

if (run$DBI) {

  print.header("Efron's diabetes data - with interactions and squared terms.")

  if (tau>0) tau = 41
  
  ## Load data.
  ## load("~/RPackage/BayesBridge/Code/C/diabetes.RData")
  data("diabetes", package="BayesBridge");
  cov.name = colnames(diabetes$x2);
  y = diabetes$y;
  X = diabetes$x2;
  n = length(y);
  cnames = cov.name
  
  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:n){ X[i,] = X[i,] - mX; }

  info <- run.it(y, X, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=print.it, plot.it=plot.it, name="DBI",
                 ortho=FALSE, betaburn=betaburn, extras=extras, inflate=inflate)

}

##------------------------------------------------------------------------------
## Boston Housing

if (run$BH) {

  print.header("Boston Housing")
    
  if (tau>0) tau   = 0.15
  
  data("BostonHousing", package="mlbench")

  ## Setup
  y = BostonHousing$medv
  n = length(y)

  ## No interactions or squared terms.
  X = model.matrix(medv ~ ., data=BostonHousing);

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:n){ X[i,] = X[i,] - mX; }
  X = X[,-1]
  
  cnames = colnames(X);

  info <- run.it(y, X, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=print.it, plot.it=plot.it, name="BH",
                 ortho=FALSE, betaburn=betaburn, extras=extras, inflate=inflate)  

}

##------------------------------------------------------------------------------
## Boston Housing with interactions

if (run$BHI) {

  print.header("Boston Housing with interactions")

  if (tau>0) tau = 0.15
  
  data("BostonHousing", package="mlbench")

  ## Setup
  y = BostonHousing$medv
  n = length(y)

  ## No interactions or squared terms.
  X = model.matrix(medv ~ ., data=BostonHousing);
  cnames = colnames(X);

  ## Include interations and squared terms.
  idc  = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13);
  nidc = length(idc)

  X.base  = model.matrix(medv ~ . * ., data=BostonHousing);

  X = X.base
  cnames = colnames(X.base);
  for (i in 1:nidc){
    nm = paste(colnames(BostonHousing)[idc[i]], 2, sep="");
    cnames = c(cnames, nm);
    X = cbind(X, BostonHousing[,idc[i]] * BostonHousing[,idc[i]]);
  }

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:n){ X[i,] = X[i,] - mX; }
  X = X[,-1]

  cnames = colnames(X);

  info <- run.it(y, X, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=FALSE, plot.it=plot.it, name="BHI",
                 ortho=FALSE, betaburn=betaburn, extras=extras, inflate=inflate)
  
}

################################################################################
                         ## ORTHOGONAL DESIGN MATRIX ##
################################################################################

##------------------------------------------------------------------------------
## Efron - orthogonalized

if (oth$EFRON) {

  print.header("Efron - orthogonalized")

  if (tau>0) tau = 41;
  
  ## Load data.
  ## load("~/RPackage/BayesBridge/Code/C/diabetes.RData")
  data("diabetes", package="BayesBridge");
  y = diabetes$y;
  X = diabetes$x;

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }

  ## Ortongonalize
  Q = qr.Q(qr(X))
  colnames(Q) = c(1:ncol(Q))

  info <- run.it(y, Q, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=TRUE, plot.it=plot.it, name="efron-ortho", ortho=TRUE)
  
}

##------------------------------------------------------------------------------
## Efron - orthogonalized with interactions

if (oth$DBI) {

  print.header("Efron - orthogonalized with interactions")

  if (tau>0) tau = 41;
  
  ## Load data.
  ## load("~/RPackage/BayesBridge/Code/C/diabetes.RData")
  data("diabetes", package="BayesBridge");
  y = diabetes$y;
  X = diabetes$x2;

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }

  ## Ortongonalize
  Q = qr.Q(qr(X))
  colnames(Q) = c(1:ncol(Q))

  info <- run.it(y, Q, nsamp=nsamp, burn=burn,
                 alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                 save.it=save.it, print.it=TRUE, plot.it=plot.it, name="DBI-ortho", ortho=TRUE)
  
}

##------------------------------------------------------------------------------
## Boston Housing - orthogonalized

if (oth$BH) {

  print.header("Boston Housing - orthogonalized")

  if (tau>0) tau = 0.15
  
  data("BostonHousing", package="mlbench")

  ## Setup
  y = BostonHousing$medv

  ## No interactions or squared terms.
  X = model.matrix(medv ~ ., data=BostonHousing);

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }
  X = X[,-1]

  ## Ortongonalize
  Q = qr.Q(qr(X))
  colnames(Q) = c(1:ncol(Q))

  info <- run.it(y, Q, nsamp=nsamp, burn=burn,
                        alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                        save.it=save.it, print.it=TRUE, plot.it=plot.it, name="BH-ortho", ortho=TRUE)  

}

##------------------------------------------------------------------------------
## Boston Housing with interactions

if (oth$BHI) {

  print.header("Boston Housing with interactions - orthogonalized")

  if (tau>0) tau=0.15
  
  data("BostonHousing", package="mlbench")

  ## Setup
  y = BostonHousing$medv

  ## No interactions or squared terms.
  X = model.matrix(medv ~ ., data=BostonHousing);
  cnames = colnames(X);

  ## Include interations and squared terms.
  idc  = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13);
  nidc = length(idc)

  X.base  = model.matrix(medv ~ . * ., data=BostonHousing);

  X = X.base
  cnames = colnames(X.base);
  for (i in 1:nidc){
    nm = paste(colnames(BostonHousing)[idc[i]], 2, sep="");
    cnames = c(cnames, nm);
    X = cbind(X, BostonHousing[,idc[i]] * BostonHousing[,idc[i]]);
  }

  ## Center things.
  y = y - mean(y);
  mX = colMeans(X);
  for(i in 1:442){ X[i,] = X[i,] - mX; }
  X = X[,-1]

  ## Ortongonalize
  Q = qr.Q(qr(X))
  colnames(Q) = c(1:ncol(Q))

  info <- run.it(y, Q, nsamp=nsamp, burn=burn,
                        alpha=alpha, sig2.shape=0.0, sig2.scale=0.0, ntrials=ntrials, tau=tau, 
                        save.it=save.it, print.it=FALSE, plot.it=plot.it, name="BHI-ortho", ortho=TRUE)
  
}

################################################################################
                                 ## APPENDIX ##
################################################################################
