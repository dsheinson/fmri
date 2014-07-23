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
pval = array(NA, c(N, 4, np, nb))
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
      
      # Box.test on OLS residuals
      pval.res[j,1,i,k] = Box.test(fit.ols$residuals)$p.val
      
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
      pval[j,3,i,k] = pt(t, length(y)-2, lower=FALSE)
      Neff = length(y)*((1-phi)/(1+phi))
      df = Neff - 2
      pval[j,4,i,k] = pt(t, df, lower=FALSE)
      
      # Fit AR(1) and AR(2) using ML and get pvalue of Box.test
      fit.ar1 = arima(y,order=c(1,0,0),xreg=X[,2])
      pval.res[j,2,i,k] = Box.test(fit.ar1$residuals)$p.val
      fit.ar2 = arima(y,order=c(2,0,0),xreg=X[,2])
      pval.res[j,3,i,k] = Box.test(fit.ar2$residuals)$p.val
      
      print(c(i,j,k))
    }
  }
}

# Save data
mypvals = list(pval=pval,pval.res=pval.res)
save(mypvals,file=paste(dpath,"simstudy_AR1.rdata",sep=""))

# Compute FPRs for all phi's
alpha = seq(0,1,0.001)
fpr = array(NA, c(length(alpha),4,np))
neg.df = numeric(np)
for(i in 1:np)
{
  ind = which(is.nan(pval[,4,i,1]))
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
  lines(alpha,fpr[,1,i]-1.96*sqrt(fpr[,1,i]*(1-fpr[,1,i])/(N-neg.df[i]*N)),lty=2)
  lines(alpha,fpr[,1,i]+1.96*sqrt(fpr[,1,i]*(1-fpr[,1,i])/(N-neg.df[i]*N)),lty=2)
  lines(alpha,fpr[,2,i],col=2)
  lines(alpha,fpr[,2,i]-1.96*sqrt(fpr[,2,i]*(1-fpr[,2,i])/(N-neg.df[i]*N)),col=2,lty=2)
  lines(alpha,fpr[,2,i]+1.96*sqrt(fpr[,2,i]*(1-fpr[,2,i])/(N-neg.df[i]*N)),col=2,lty=2)
  lines(alpha,fpr[,3,i],col=3)
  lines(alpha,fpr[,3,i]-1.96*sqrt(fpr[,3,i]*(1-fpr[,3,i])/(N-neg.df[i]*N)),col=3,lty=2)
  lines(alpha,fpr[,3,i]+1.96*sqrt(fpr[,3,i]*(1-fpr[,3,i])/(N-neg.df[i]*N)),col=3,lty=2)
  lines(alpha,fpr[,4,i],col=4)
  lines(alpha,fpr[,4,i]-1.96*sqrt(fpr[,4,i]*(1-fpr[,4,i])/(N-neg.df[i]*N)),col=4,lty=2)
  lines(alpha,fpr[,4,i]+1.96*sqrt(fpr[,4,i]*(1-fpr[,4,i])/(N-neg.df[i]*N)),col=4,lty=2)
  abline(0,1,lwd=2,col='gray')
  mtext(paste("Negative DF rate:",neg.df[i]),side=3,cex=0.85)
  if(i == 1) legend("topleft",c("OLS","PW","REML","REMLc","95% CI",expression(alpha)),lty=c(1,1,1,1,2,1),col=c(1,2,3,4,1,'gray'),lwd=c(1,1,1,1,1,2),cex=1.4)
}
dev.off()

# Make FPR table
ind = which(alpha %in% c(0.001,0.01,0.05))
round(fpr[ind,,],3)

# Calculate true postive rates curves
tpr = array(NA, c(length(alpha),4,np,nb-1))
neg.df = matrix(NA,nr=np,nc=nb-1)
for(k in 2:nb)
{
  for(i in 1:np)
  {
    ind = which(is.nan(pval[,4,i,k]))
    neg.df[i,k-1] = length(ind) / N
    if(length(ind) > 0) pval.df = pval[-ind,,i,k] else pval.df = pval[,,i,k]
    for(j in 1:length(alpha)) tpr[j,,i,k-1] = apply(pval.df,2,function(x) mean(x < alpha[j]))
  }
}

# Plot ROC curves: 4 x 4 plots, beta = 3, 6, 9
beta = unique(sims[[2]]$beta)[-1]
for(k in 1:(nb-1))
{
  pdf(file=paste(gpath,"simstudy-ROC-",beta[k],".pdf",sep=""),width=10,height=10)
  par(mfrow=c(2,2),mar=c(5,6,4,2)+0.1)
  for(i in 1:np)
  {
    plot(fpr[,1,i],tpr[,1,i,k],type="l",xlim=c(0,1),ylim=c(0,1),axes=ifelse(i==1,TRUE,FALSE),xlab=ifelse(i==1,"FPR",""),ylab=ifelse(i==1,"TPR",""),main=eval(bquote(expression(paste(phi,"=",.(phi[i]))))),cex.lab=1.75,cex.main=1.75)
    if(i != 1) box()
    lines(fpr[,2,i],tpr[,2,i,k],col=2)
    lines(fpr[,3,i],tpr[,3,i,k],col=3)
    lines(fpr[,4,i],tpr[,4,i,k],col=4)
    abline(0,1,lty=2)
    mtext(paste("Negative DF rate:",neg.df[i,k]),side=3,cex=0.85)
    legend("bottomright",c("OLS","PW","REML","REMLc"),lty=c(1,1,1,1),col=c(1,2,3,4),cex=1.5)
  }
  dev.off()
}

# Calculate proportion of whitened residuals
pwn = array(NA, c(length(alpha),3,np,nb))
for(k in 1:nb)
{
  for(i in 1:np)
  {
    for(j in 1:length(alpha)) pwn[j,,i,k] = apply(pval.res[,,i,k],2,function(x) mean(x > alpha[j]))
  }
}

# PWN table
ind = which(alpha %in% c(0.001,0.01,0.05))
pwn[ind,,,1]
