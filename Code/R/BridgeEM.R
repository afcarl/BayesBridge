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

bridge.EM.R = function(y, X, alpha, ratio=1.0, lambda.max=1e9*ratio, tol=1e-9, max.iter=30, init=NULL){
    X <- as.matrix(X)
    xx <- t(X)%*%X
    xy <- t(X)%*%y
    ixx <- chol2inv(chol(xx))

    p <- ncol(X)

    bhat <- drop(ixx%*%xy)
    Beta = bhat
    if (!is.null(init)) Beta = init

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

    out = list("beta"=drop(Beta), "iter"=iter, "diff"=diff)
    out
}
