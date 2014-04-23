

## dyn.unload("BayesBridge.so"); dyn.load("BayesBridge.so"); source("../C/BridgeWrapper.R")

test.weird <- function()
{
    rtexp(1, 0, Inf, 1.0)
    rtexp(1, Inf, Inf, 1.0)
    rtexp(1, Inf, 1.0, 1.0)
    rtexp(1, NaN, Inf, 1.0)
    rtexp(1, NA, 1.0, 1.0)
    rtexp(1, 0, 1, NaN)
    rtexp(1, 0, 1, NA)

    rtnorm(num=1, mu=0.0, sig=1.0, left=-1.0, right=1.0)
    rtnorm(num=1, mu=0.0, sig=1.0, left=-Inf, right=0.0)
    rtnorm(num=1, mu=0.0, sig=1.0, left=0.0, right=Inf)
    rtnorm(num=1, mu=0.0, sig=1.0, left=Inf, right=0.0)
    rtnorm(num=1, mu=0.0, sig=1.0, left=0.0, right=-Inf)
    rtnorm(num=1, mu=0.0, sig=1.0, left=-Inf, right=Inf)

}

texpon.left.test <- function()
{
    n    = 1e3

    left = 1
    rate = 2

    set.seed(100)
    x1 = rtexp.left(n, as.double(left),as.double(rate))

    set.seed(100)
    x2 = rexp(n, rate=rate) + left

    cbind(head(x1), head(x2))
    
}

texpon.test <- function()
{
    n     = 1e6

    left  = 1
    right = 2
    rate  = 1.0
    x1    = rtexp(n, left, right, rate)

    x2    = rexp(10*n, rate=rate) + left
    x2    = x2[x2<right]

    summary(x1)
    summary(x2)

}

tnorm.test <- function()
{
    require("truncnorm")
    n = 1e7
    N = 10
    left.grid = seq(-3, 3, length.out=N)

    cat("Summary statistics are for difference between trunacted\nnorm samplers using BayesBridge and truncnorm package.\n")
    
    for (i in 1:N) {
        out1 = rtnorm.left(n, left=left.grid[i], 0.0, 1.0)
        out2 = rtruncnorm(n, a=left.grid[i], b=Inf, 0.0, 1.0)
        print(summary(out1) - summary(out2))
        cat("........................................\n")
    }

    for (i in 1:(N-1)) {
        for (j in (i+1):N) {
            out1 = rtnorm.both(n, left=left.grid[i], right=left.grid[j], 0.0, 1.0)
            out2 = rtruncnorm(n, a=left.grid[i], b=left.grid[j], 0.0, 1.0)
            print(summary(out1) - summary(out2))
            cat("........................................\n")
        }
    }
    
}
