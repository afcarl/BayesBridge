
source("~/RV-Project/Code/C_Examples/MyLib/Matrix/MatrixIO.R");

# Read the files.
beta.small = drop(read.mat("beta.small"));
beta.sig2 = drop(read.mat("beta.sig2"));
beta.full = drop(read.mat("beta.full"));

stats = array(0, dim=c(2, 3, 10))

for(i in 1:10){

#h.small = hist(beta.small[i,5001:20000], breaks=50);
#h.sig2 = hist(beta.sig2[i,5001:20000], breaks=50);
#h.full = hist(beta.full[i,5001:20000], breaks=50);

stats[1,1,i] = mean(beta.small[i, 5001:20000]);
stats[1,2,i] = mean(beta.sig2[i, 5001:20000]);
stats[1,3,i] = mean(beta.full[i, 5001:20000]);

stats[2,1,i] = sd(beta.small[i, 5001:20000]);
stats[2,2,i] = sd(beta.sig2[i, 5001:20000]);
stats[2,3,i] = sd(beta.full[i, 5001:20000]);

}
