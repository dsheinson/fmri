# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data
load(paste(dpath,"fmri-ar-sim-M010.rdata",sep=""))

# Fit data using OLS, prewhitening, and REML
require(nlme)
N = 1000
np = 4
nb = 4
M = 3
pval = array(NA, c(N, M, np, nb))
for(i in 1:np)
{
  for(j in 1:N)
  {
    for(k in 1:nb)
    {
      ind = (i-1)*nb + k
      y = sims[[1]][[ind]][[j]]$y
      X = sims[[1]][[ind]][[j]]$true.params$U
      
      # Fit OLS - t-test
      fit.ols = lm(y ~ X - 1)
      df = length(y) - 2
      t = fit.ols$coef[2] / summary(fit.ols)$coef[2,2]
      pval[j,1,i,k] = pt(t, df, lower=FALSE)
      
      # Fit prewhitened data - t-test
      fit.res = arima(fit.ols$residuals, order = c(1,0,0), include.mean=FALSE)
      V = fit.res$coef^as.matrix(dist(0:(length(y)-1)))
      Sinv = solve(t(chol(as.matrix(V))))
      y.trans = (Sinv%*%y)[,1]
      X.trans = Sinv%*%X
      fit.pw = lm(y.trans ~ X.trans - 1)     
      df = length(y) - 2
      t = fit.pw$coef[2] / summary(fit.pw)$coef[2,2]
      pval[j,2,i,k] = pt(t, df, lower=FALSE)
      
      # Fit using REML - t-test with adjusted df
      fit.reml = gls(y ~ X - 1, correlation=corAR1())
      phi = corMatrix(fit.reml$modelStruct$corStruct)[1,2]
      t = fit.reml$coef[2] / sqrt(fit.reml$varBeta[2,2])
      Neff = length(y)*((1-phi)/(1+phi))
      df = Neff - 2
      pval[j,3,i,k] = pt(t, df, lower=FALSE)
      
      print(c(i,j,k))
    }
  }
}

# Compute FPRs for all phi's
alpha = seq(0,1,0.001)
fpr = array(NA, c(length(alpha),M,np))
neg.df = numeric(np)
for(i in 1:np)
{
  ind = which(is.nan(pval[,3,i,1]))
  neg.df[i] = length(ind) / N
  if(length(ind) > 0) pval.df = pval[-ind,,i,1] else pval.df = pval[,,i,1]
  for(j in 1:length(alpha)) fpr[j,,i] = apply(pval.df,2,function(x) mean(x < alpha[j]))
}

# Plot false positive rates versus alpha
phi = unique(sims[[2]]$phi)
ycut = which(alpha == 0.2)
pdf(file=paste(gpath,"simstudy-FPR.pdf",sep=""),width=10,height=10)
par(mfrow=c(2,2),mar=c(5,6,4,2)+0.1)
for(i in 1:np)
{
  plot(alpha,fpr[,1,i],type="l",xlim=c(0,.2),ylim=c(0,max(fpr[1:ycut,,])),axes=ifelse(i==1,TRUE,FALSE),xlab=ifelse(i==1,expression(alpha),""),ylab=ifelse(i==1,"FPR",""),main=eval(bquote(expression(paste(phi,"=",.(phi[i]))))),cex.lab=1.75,cex.main=1.75)
  if(i != 1) box()
  lines(alpha,fpr[,2,i],col=2)
  lines(alpha,fpr[,3,i],col=4)
  abline(0,1,lty=2)
  mtext(paste("Negative DF rate:",neg.df[i]),side=3,cex=0.85)
  if(i == 1) legend("bottomright",c("OLS","PW","REML"),lty=c(1,1,1),col=c(1,2,4),cex=1.5)
}
dev.off()

# Make FPR table
ind = which(alpha %in% c(0.001,0.01,0.05))
round(fpr[ind,,],3)

# Calculate true postive rates curves
tpr = array(NA, c(length(alpha),M,np,nb-1))
neg.df = matrix(NA,nr=np,nc=nb-1)
for(k in 2:nb)
{
  for(i in 1:np)
  {
    ind = which(is.nan(pval[,3,i,k]))
    neg.df[i,k-1] = length(ind) / N
    if(length(ind) > 0) pval.df = pval[-ind,,i,k] else pval.df = pval[,,i,k]
    for(j in 1:length(alpha)) tpr[j,,i,k-1] = apply(pval.df,2,function(x) mean(x < alpha[j]))
  }
}

# Plot ROC curve for beta = 5 and phi = 0.95
beta = unique(sims[[2]]$beta)[-1]
k = 1; i = 4
pdf(file=paste(gpath,"simstudy-ROC.pdf",sep=""))
par(mar=c(5,6,4,2)+0.1)
plot(fpr[,1,i],tpr[,1,i,k],type="l",xlim=c(0,1),ylim=c(0,1),xlab="FPR",ylab="TPR",main=eval(bquote(expression(paste(phi,"=",.(phi[i]),", ",beta[1],"=",.(beta[k]))))),cex.lab=1.75,cex.main=1.75)
lines(fpr[,2,i],tpr[,2,i,k],col=2)
lines(fpr[,3,i],tpr[,3,i,k],col=4)
abline(0,1,lty=2)
mtext(paste("Negative DF rate:",neg.df[i,k]),side=3,cex=0.85)
legend("bottomright",c("OLS","PW","REML"),lty=c(1,1,1),col=c(1,2,4),cex=1.5)
dev.off()