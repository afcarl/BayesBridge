source("bridge-MCMC.R")

library(lars)

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

mylabs = c("Age", "Female", "Body Mass Index", "Blood Pressure", "TC", "LDL", "HDL", "TCH", "LTG", "Glucose")

#write.table(X, "X.csv", sep=",", eol="\n", row.names=FALSE, col.names=TRUE)
#write.table(Y, file="Y.csv", eol="\n", row.names=FALSE, col.names=FALSE)

### Now the MAP

alpha = 1/2

# GridSize = 200
GridSize = 10;
#NuGrid = c(seq(0.001,0.02*sum(abs(bhat)),length=GridSize/2), seq(0.02*sum(abs(bhat)),2*sum(abs(bhat)),length=GridSize/2))
NuGrid = sort(c(seq(-1,4.5,length=GridSize)))
NuGrid = rev(10^NuGrid)
GCV = rep(0, GridSize)
BetaMatrix = matrix(ncol=p, nrow=GridSize)

for(i in 1:GridSize)
{
	if(i == 1)
	{
		Beta = bhat
	}
	else
	{
		## Warm start at the previous solution
		## But add some jitter to the initial guess avoid screwy solution paths
		Beta = BetaMatrix[i-1,] + rnorm(p,0,abs(bhat)/10)
		#Beta = rep(1,p)
	}
	diff = 1
	Nu = NuGrid[i]
	tau = (Nu)^{-1/alpha}
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
	GCV[i] = sum( ((Y - YHat)/(1-sum(diag(S))/n))^2 )
	BetaMatrix[i,] = Beta
}

mybest = which.min(GCV)
bhat.bridge = round(BetaMatrix[mybest,],3)
nuhat = NuGrid[mybest]


BetaLa = abs(BetaMatrix)^alpha
BetaLa = apply(BetaLa, 1, sum)
plot(BetaLa/sum(abs(bhat)^alpha), GCV)


plot(0, xlim=c(0,1), ylim=1.1*range(BetaMatrix), bty='n', type='n', xlab="|Beta|^a / |Beta_OLS|^a", ylab="Coefficients", main="Solution path: Bridge (a=0.5) with Student-t error")
for(i in 1:p)
{
	lines(BetaLa/sum(abs(bhat)^alpha), BetaMatrix[,i])
}
abline(v = BetaLa[mybest]/sum(abs(bhat)^alpha), col='red')



plot(log10(NuGrid), BetaLa/sum(abs(bhat)^alpha), xlim=rev(range(log10(NuGrid))), bty='n', type='l', xlab="Log10(nu)", ylab="|Beta|^a / |Beta_OLS|^a", main="Alpha norm versus nu")


plot(0, xlim=rev(range(log10(NuGrid))), ylim=1.1*range(BetaMatrix), bty='n', type='n', xlab="Log10(nu)", ylab="Coefficients", main="Solution path: Bridge (a=0.5)")
for(i in 1:p)
{
	lines(log10(NuGrid), BetaMatrix[,i])
}
abline(v=log10(NuGrid[mybest]), col='red')



#### Now the MCMC

mygrid=sort(union(seq(-800,800,length=1001), seq(-10,10,length=201)))


ans = bridge(X,Y,20000,0.5,mygrid=mygrid,verbose=500)


pdf("diabetes-betas.pdf", height=10, width=7.5)
par(mfrow=c(5,2), mar=c(4,4,2,1), las=1)
bhat.bayesmode = rep(0,p)
for(j in 1:10)
{
	hist(ans$beta[,j],100, prob=TRUE, col='white',border='white',
		ylim=c(0, max(ans$betadens[,j])),
		#xlim=c(max(min(mygrid), min(ans$beta[,j])), min(max(mygrid), max(ans$beta[,j]))),
		xlim=c(min(mygrid[ans$betadens[,j] > 1e-6]), max(mygrid[ans$betadens[,j] > 1e-6])),
		main = mylabs[j], xlab="", ylab="")
	lines(mygrid,ans$betadens[,j])
	bhat.bayesmode[j] = mygrid[which.max(ans$betadens[,j])]
	abline(v=mean(ans$beta[,j]), col='blue', lty='dashed')
	abline(v=mygrid[which.max(ans$betadens[,j])], col='black', lty="dotted")
		abline(v=bhat.bridge[j], col='red')
}
dev.off()





bhat.bayes = round(apply(ans$beta,2,mean),3)

#res = 200
#zz1 = sort(unique(c(0, seq(-650,50,length=res), seq(-490,-450,length=res),seq(-190,-150,length=res))))
#zz2 = sort(unique(c(0, seq(-50,250,length=res), seq(260,300,length=res), seq(-20,20,length=res))))



res = 500
zz1= sort(unique(c(0, seq(-750,150,length=res))))
zz2= sort(unique(c(0, seq(-150,350,length=res))))

ProfLike1 = matrix(0, nrow=length(zz1),ncol=length(zz2))
ProfLike2 = ProfLike1
ProfLike3 = ProfLike1
for(i in 1:length(zz1))
{
	for(j in 1:length(zz2))
	{
		betatry1 = bhat.bridge
		betatry1[c(5,6)] = c(zz1[i], zz2[j])
		betatry2 = bhat.bayes
		betatry2[c(5,6)] = c(zz1[i], zz2[j])
		betatry3 = bhat.bayesmode
		betatry3[c(5,6)] = c(zz1[i], zz2[j])
		YHat1 = X %*% betatry1
		YHat2 = X %*% betatry2
		YHat3 = X %*% betatry3
		ProfLike1[i,j] = -0.5*sum((Y-YHat1)^2) - nuhat*sum(abs(betatry1)^alpha)
		ProfLike2[i,j] = -0.5*sum((Y-YHat2)^2) - nuhat*sum(abs(betatry2)^alpha)
		ProfLike3[i,j] = -0.5*sum((Y-YHat3)^2) - nuhat*sum(abs(betatry3)^alpha)
	}
}


par(mfrow=c(2,2))

contour(zz1,zz2,ProfLike1, nlevels=250)
points(bhat.bridge[5], bhat.bridge[6], col='red', pch=19)
points(bhat.bayes[5], bhat.bayes[6], col='blue', pch=19)
points(bhat.bayesmode[5], bhat.bayes[6], col='darkgreen', pch=19)

contour(zz1,zz2,ProfLike2, nlevels=250)
points(bhat.bridge[5], bhat.bridge[6], col='red', pch=19)
points(bhat.bayes[5], bhat.bayes[6], col='blue', pch=19)
points(bhat.bayesmode[5], bhat.bayes[6], col='darkgreen', pch=19)

contour(zz1,zz2,ProfLike3, nlevels=250)
points(bhat.bridge[5], bhat.bridge[6], col='red', pch=19)
points(bhat.bayes[5], bhat.bayes[6], col='blue', pch=19)
points(bhat.bayesmode[5], bhat.bayes[6], col='darkgreen', pch=19)
