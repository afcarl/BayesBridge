######################################################################

## Copyright 2012 Nicholas G. Polson, James G. Scott, Jesse Windle
## Contact info: <jwindle@ices.utexas.edu>.

## This file is part of BayesBridge.

## BayesBridge is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
  
## BayesBridge is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
  
## You should have received a copy of the GNU General Public License
## along with BayesBridge.  If not, see <http:##www.gnu.org/licenses/>.
			      
######################################################################

log.tau.grid = seq(-10, 0, 0.5);

trace.beta <- function(y, X, alpha=0.5, ratio.grid=exp(seq(-20,20,0.1)),
                       tol=1e-9, max.iter=30, use.cg=FALSE, plot.it=FALSE)
{
  X = as.matrix(X);
  N = dim(X)[1];
  P = dim(X)[2];
  L = length(ratio.grid);
  beta = array(0, dim=c(L, P));

  colnames(beta) = colnames(X)
  
  for (i in 1:L) {
    beta[i,] = bridge.EM(y, X, alpha, ratio=ratio.grid[i],
                         lambda.max=ratio.grid[i]/tol, tol, max.iter, use.cg);
  }

  log.grid = log(ratio.grid);
  width = log.grid[L] - log.grid[1];
  ymin = min(beta);
  ymax = max(beta);
  
  plot(log.grid, beta[,1], col=1, type="l",
       ylim=c(ymin, ymax), xlim=c(log.grid[1]-0.1*width, log.grid[L]),
       ylab="Coefficient", xlab="Log Ratio",
       main="Coefficients vs. log(ratio)");
  if (P > 1) {
    for (i in 2:P) {
      lines(log.grid, beta[,i], col=i, lty=i/8+1);
    }
  }
  
  legend("bottomleft", legend=colnames(beta), col=seq(1:P), lty=seq(1:P)/8+1);
  
  list("beta"=beta, "grid"=ratio.grid, "log.grid"=log.grid)
}


