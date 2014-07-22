# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load Craig's data
conv = read.csv(paste(dpath,"basis.csv",sep=""),header=FALSE)[,1]
hrf = read.csv(paste(dpath,"hrf.csv",sep=""),header=FALSE)[,1]
regions = c("region-frontalpole-left","region-IPS-left","region-IPS-right","region-primaryvisual","region-secondaryvisual-left","region-secondaryvisual-right")
nt = length(conv)
Nr = length(regions)
y = array(NA, c(125, nt+6, Nr))
for(i in 1:Nr) y[,,i] = as.matrix(read.csv(paste(dpath,regions[i],".csv",sep=""),header=F))

# Plot fMRI time series versus convolution and hrf
pdf(file=paste(gpath,"fmri-craig-data..pdf",sep=""))
par(mfrow=c(3,1),mar=c(5,6,4,2))
plot(1:nt,y[27,-(1:6),2],type="l",xlab="",ylab=expression(y[t]),cex.lab=2)
title("Voxel 27 in IPS-left",cex.main=2)
plot(1:nt,conv,type="l",xlab="",ylab=expression(conv[t]),cex.lab=2)
title(eval(bquote(expression(paste(beta," = (",.(paste(sim$true.params$theta[1:d],sep="",collapse=",")),")",", ",phi," = ",.(sim$true.params$theta[d+1]),", ",sigma[s]^2," = ",.(sim$true.params$theta[d+2]),", ",sigma[m]^2," = ",.(sim$true.params$theta[d+3]))))),cex.main=2)
plot(1:length(hrf),hrf,type="l",xlab="TR",ylab="hrf",cex.lab=2)
dev.off()

# Load ARMA comparison data
arma.comp = read.csv(paste(dpath,"arma.model.comparison.csv",sep=""))
avg.ar <- function(region, crit) mean(arma.comp[which(arma.comp$Region == region & arma.comp$Criterion == crit),5])
avg.ma <- function(region, crit) mean(arma.comp[which(arma.comp$Region == region & arma.comp$Criterion == crit),6])
mydata = expand.grid(region=c("fontalpole-left","IPS-left","IPS-right","primaryvisual","secondaryvisual-left","secondaryvisual-right"),crit=c("AICC","AIC","BIC"),stringsAsFactors=F)
require(plyr)
mean.ar = maply(mydata, avg.ar)
mean.ma = maply(mydata, avg.ma)