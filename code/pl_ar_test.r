source("fmri_ar_functions.r")
source("pl_ar_functions.r")
source("pl.r")
source("pf_functions.r")
source("invgamma.r")
require(dlm)

# Set data and graphics paths
dpath = "../data/"
gpath = "../graphs/"

# Load basis function
conv = read.csv(paste(dpath,"basis.csv",sep=""),header=FALSE)[,1]
U = cbind(1,conv)
nt = length(conv)

# Simulate observations from M101, M011, or M111
beta = c(750, 15)
phi = 0.9
sigma2s = 16
sigma2m = 9
FF = conv
sim = dlm.sim(nt, beta, phi, sigma2s, sigma2m, U, FF)
windows()
par(mfrow=c(3,1),mar=c(5,6,4,2))
plot(1:nt,sim$y,type="l",xlab=expression(t),ylab=expression(y[t]),cex.lab=2)
title(eval(bquote(expression(paste("Simulated from ",M['101'],sep="")))),cex.main=2)
plot(1:nt,conv,type="l",xlab="",ylab=expression(conv[t]),main=substitute(paste(beta[0],"=",a0,",",beta[1],"=",a1,",",phi,"=",a2,",",sigma[s]^2,"=",a3,",",sigma[m]^2,"=",a4),list(a0=sim$true.params$theta[1],a1=sim$true.params$theta[2],a2=sim$true.params$theta[3],a3=sim$true.params$theta[4],a4=sim$true.params$theta[5])),cex.main=2,cex.lab=2)
plot(0:nt,sim$x,type="l",xlab="",ylab=expression(x[t]),main="",cex.lab=2)

# Fit data to all models using particle learning
n = 100
abm = ig.mom(sigma2m, 500)
abs = ig.mom(sigma2s, 500)
prior = list(b0=c(beta[1], beta[0]),B0=matrix(c(1000,0,0,225),nr=2), am0=abm[1], bm0=abm[2], phi0=phi, Phi0=1, as0=abs[1], bs0=abs[2])
rprior <- function(j) rprior.pl(prior)
F1 = F2 = F3 = function(q, s, w) sum(w*pnorm(q, s[1,], sqrt(s[2,])))
F4 <- function(q, s, w)
{
  np = length(w)
  b = numeric(np)
  for(j in 1:np) b[j] = .5*s[4,j] + .5*(prior$phi0^2)/prior$Phi0 - .5*(s[1,j]^2)/s[2,j] + prior$bm0
  return(sum(w*pinvgamma(q, s[3,], b)))
}
F5 <- function(q, s, w)
{
  np = length(w)
  b = numeric(np)
  for(j in 1:np) b[j] = .5*s[8,j] + .5*t(prior$b0)%*%solve(prior$B0)%*%prior$b0 - .5*t(s[1:2,j])%*%solve(matrix(s[3:6,j],nr=2))%*%s[1:2,j] + prior$bm0
  return(sum(w*pinvgamma(q, s[7,], b)))
}
Fq=list(F1,F2,F3,F4,F5) # To compute quantiles from mixture distribution
int = list(c(-1000,1000),c(-1000,1000),c(-1000,1000),c(0,1000),c(0,1000))

# M101
FF = conv
fit = dlmMLE(sim$y, c(phi, log(sigma2s), log(sigma2m)), function(par) build.ar1(par, U, FF))
s = dlmSmooth(dlmFilter(sim$y, build.ar1(fit$par, U, FF)))
abm = ig.mom(exp(fit$par[3]), 500)
abs = ig.mom(exp(fit$par[2]), 500)
prior$b0 = s$s[nt+1,1:2]; prior$phi0 = fit$par[1]; prior$am0 = abm[1]; prior$bm0 = abm[2]; prior$as0 = abs[1]; prior$bs0 = abs[2]
rprior <- function(j) rprior.pl(prior)
dlpred <- function(y, x, suff.x, theta) dlpred.pl(y, x, theta, U, FF)
revo <- function(y, x, suff.x, theta) revo.pl(y, x, theta, U, FF)
rmove <- function(suff.theta) rmove.pl(suff.theta, prior)
smap.theta <- function(suff.theta, y, x) smap.theta.pl(suff.theta, y, x, U, FF)
smap.state <- function(suff.x, y, theta) smap.state.pl(suff.x, y, theta)
out = pl(sim$y, dlpred, revo, rprior, rmove, smap.theta, smap.state, n, TRUE, method="stratified", threshold=0.8, nonuniformity="ess", log=F)
out.s = pl.smooth(out$state, out$suff.theta, out$weight[,nt+1], dlevo.pl, rmove, method="stratified", nonuniformity="ess", threshold=0.8, log=F)
state.quant = pf.quantile(out$state, out.s$weight, function(x, j) x, c(.025,.975))
tq = list(out$suff.theta[c(1,3),,], out$suff.theta[c(2,6),,], out$suff.theta[9:10,,], out$suff.theta[9:12,,], out$suff.theta[1:8,,])
theta.quant = pf.mix.quantile(tq, out$weight, Fq, int=int, print.iter=T)

# M011
FF = rep(1,nt)
fit = dlmMLE(sim$y, c(phi, log(sigma2s), log(sigma2m)), function(par) build.ar1(par, U, FF))
s = dlmSmooth(dlmFilter(sim$y, build.ar1(fit$par, U, FF)))
abm = ig.mom(exp(fit$par[3]), 500)
abs = ig.mom(exp(fit$par[2]), 500)
prior$b0 = s$s[nt+1,1:2]; prior$phi0 = fit$par[1]; prior$am0 = abm[1]; prior$bm0 = abm[2]; prior$as0 = abs[1]; prior$bs0 = abs[2]
rprior <- function(j) rprior.pl(prior)
dlpred <- function(y, x, suff.x, theta) dlpred.pl(y, x, theta, U, FF)
revo <- function(y, x, suff.x, theta) revo.pl(y, x, theta, U, FF)
rmove <- function(suff.theta) rmove.pl(suff.theta, prior)
smap.theta <- function(suff.theta, y, x) smap.theta.pl(suff.theta, y, x, U, FF)
smap.state <- function(suff.x, y, theta) smap.state.pl(suff.x, y, theta)
out2 = pl(sim$y, dlpred, revo, rprior, rmove, smap.theta, smap.state, n, TRUE, method="stratified", threshold=0.8, nonuniformity="ess", log=F)
out2.s = pl.smooth(out2$state, out2$suff.theta, out2$weight[,nt+1], dlevo.pl, rmove, method="stratified", nonuniformity="ess", threshold=0.8, log=F)
state.quant2 = pf.quantile(out2$state, out2.s$weight, function(x, j) x, c(.025,.975))
tq = list(out2$suff.theta[c(1,3),,], out2$suff.theta[c(2,6),,], out2$suff.theta[9:10,,], out2$suff.theta[9:12,,], out2$suff.theta[1:8,,])
theta.quant2 = pf.mix.quantile(tq, out2$weight, Fq, int=int, print.iter=T)

# M111
FF = conv + 1
fit = dlmMLE(sim$y, c(phi, log(sigma2s), log(sigma2m)), function(par) build.ar1(par, U, FF))
s = dlmSmooth(dlmFilter(sim$y, build.ar1(fit$par, U, FF)))
abm = ig.mom(exp(fit$par[3]), 500)
abs = ig.mom(exp(fit$par[2]), 500)
prior$b0 = s$s[nt+1,1:2]; prior$phi0 = fit$par[1]; prior$am0 = abm[1]; prior$bm0 = abm[2]; prior$as0 = abs[1]; prior$bs0 = abs[2]
rprior <- function(j) rprior.pl(prior)
dlpred <- function(y, x, suff.x, theta) dlpred.pl(y, x, theta, U, FF)
revo <- function(y, x, suff.x, theta) revo.pl(y, x, theta, U, FF)
rmove <- function(suff.theta) rmove.pl(suff.theta, prior)
smap.theta <- function(suff.theta, y, x) smap.theta.pl(suff.theta, y, x, U, FF)
smap.state <- function(suff.x, y, theta) smap.state.pl(suff.x, y, theta)
out3 = pl(sim$y, dlpred, revo, rprior, rmove, smap.theta, smap.state, n, TRUE, method="stratified", threshold=0.8, nonuniformity="ess", log=F)
out3.s = pl.smooth(out3$state, out3$suff.theta, out3$weight[,nt+1], dlevo.pl, rmove, method="stratified", nonuniformity="ess", threshold=0.8, log=F)
state.quant3 = pf.quantile(out3$state, out3.s$weight, function(x, j) x, c(.025,.975))
tq = list(out3$suff.theta[c(1,3),,], out3$suff.theta[c(2,6),,], out3$suff.theta[9:10,,], out3$suff.theta[9:12,,], out3$suff.theta[1:8,,])
theta.quant3 = pf.mix.quantile(tq, out3$weight, Fq, int=int, print.iter=T)

# Plot sequential credible intervals
windows(width=15,height=10)
par(mfrow=c(2,3), mar=c(5,6,4,2)+0.1)
burn = 10
ymin = apply(theta.quant[-(1:burn),,1], 2, min)
ymax = apply(theta.quant[-(1:burn),,2], 2, max)
ymin2 = apply(theta.quant2[-(1:burn),,1], 2, min)
ymax2 = apply(theta.quant2[-(1:burn),,2], 2, max)
ymin3 = apply(theta.quant3[-(1:burn),,1], 2, min)
ymax3 = apply(theta.quant3[-(1:burn),,2], 2, max)

plot(0:nt, rep(sim$true.params$theta[1], nt+1), type="l", col='gray', lwd=2, ylim=c(min(ymin[1],ymin2[1],ymin3[1]),max(ymax[1],ymax2[1],ymax3[1])), xlab=expression(t), ylab=expression(beta[0]), cex.lab=2)
#title(substitute(paste(B0['1,1'],"=",a0,",",B0['2,2'],"=",a1,",",Phi[0],"=",a2),list(a0=round(prior$B0[1,1],2),a1=round(prior$B0[2,2],2),a2=round(prior$Phi0,2))))
mtext(substitute(paste(M['101'],": ",a0,sep=""),list(a0=round(pf.lmarglik(out),3))),side=3,cex=.85)
lines(0:nt, theta.quant[,1,1])
lines(0:nt, theta.quant[,1,2])
lines(0:nt, theta.quant2[,1,1],col=2)
lines(0:nt, theta.quant2[,1,2],col=2)
lines(0:nt, theta.quant3[,1,1],col=4)
lines(0:nt, theta.quant3[,1,2],col=4)
plot(0:nt, rep(sim$true.params$theta[2], nt+1), type="l", col='gray', lwd=2, ylim=c(min(ymin[2],ymin2[2],ymin3[2]),max(ymax[2],ymax2[2],ymax3[2])), xlab="", ylab=expression(beta[1]), main = expression(paste("Simulated from ",M['101'],sep="")), col.main=1, cex.lab=2, cex.main=2)
mtext(substitute(paste(M['011'],": ",a0,sep=""),list(a0=round(pf.lmarglik(out2),3))),side=3,col=2,cex=0.85)
lines(0:nt, theta.quant[,2,1])
lines(0:nt, theta.quant[,2,2])
lines(0:nt, theta.quant2[,2,1],col=2)
lines(0:nt, theta.quant2[,2,2],col=2)
lines(0:nt, theta.quant3[,2,1],col=4)
lines(0:nt, theta.quant3[,2,2],col=4)
plot(0:nt, rep(sim$true.params$theta[3], nt+1), type="l", col='gray', lwd=2, ylim=c(min(ymin[3],ymin2[3],ymin3[3]),max(ymax[3],ymax2[3],ymax3[3])), xlab="", ylab=expression(phi), cex.lab=2)
#title(substitute(paste(a['m0'],"=",a3,",",b['m0'],"=",a4,",",a['s0'],"=",a5,",",b['s0'],"=",a6),list(a3=round(prior$am0,2),a4=round(prior$bm0,2),a5=round(prior$as0,2),a6=round(prior$bs0,2))))
mtext(substitute(paste(M['111'],": ",a0,sep=""),list(a0=round(pf.lmarglik(out3),3))),side=3,col=4,cex=0.85)
lines(0:nt, theta.quant[,3,1])
lines(0:nt, theta.quant[,3,2])
lines(0:nt, theta.quant2[,3,1],col=2)
lines(0:nt, theta.quant2[,3,2],col=2)
lines(0:nt, theta.quant3[,3,1],col=4)
lines(0:nt, theta.quant3[,3,2],col=4)
plot(0:nt, rep(sim$true.params$theta[4], nt+1), type="l", col='gray', lwd=2, ylim=c(min(ymin[4],ymin2[4],ymin3[4]),max(ymax[4],ymax2[4],ymax3[4])), xlab="", ylab=expression(sigma[s]^2), main = "", cex.lab=2)
mtext(substitute(paste(beta[0],"=",a0,",",beta[1],"=",a1,",",phi,"=",a2,",",sigma[s]^2,"=",a3,",",sigma[m]^2,"=",a4),list(a0=sim$true.params$theta[1],a1=sim$true.params$theta[2],a2=sim$true.params$theta[3],a3=sim$true.params$theta[4],a4=sim$true.params$theta[5])),side=3,cex=0.85)
lines(0:nt, theta.quant[,4,1])
lines(0:nt, theta.quant[,4,2])
lines(0:nt, theta.quant2[,4,1],col=2)
lines(0:nt, theta.quant2[,4,2],col=2)
lines(0:nt, theta.quant3[,4,1],col=4)
lines(0:nt, theta.quant3[,4,2],col=4)
plot(0:nt, rep(sim$true.params$theta[5], nt+1), type="l", col='gray', lwd=2, ylim=c(min(ymin[5],ymin2[5],ymin3[5]),max(ymax[5],ymax2[5],ymax3[5])), xlab="", ylab=expression(sigma[m]^2), main = "", cex.lab=2)
mtext(substitute(paste(B0['1,1'],"=",a2,",",B0['2,2'],"=",a3,",",Phi[0],"=",a5),list(a2=prior$B0[1,1],a3=prior$B0[2,2],a5=prior$Phi0)),side=3,cex=0.85)
lines(0:nt, theta.quant[,5,1])
lines(0:nt, theta.quant[,5,2])
lines(0:nt, theta.quant2[,5,1],col=2)
lines(0:nt, theta.quant2[,5,2],col=2)
lines(0:nt, theta.quant3[,5,1],col=4)
lines(0:nt, theta.quant3[,5,2],col=4)
plot(0:nt, sim$x, type="l", col='gray', lwd=2, ylim=c(min(sim$x,state.quant[-(1:burn),1,1],state.quant2[-(1:burn),1,1],state.quant3[-(1:burn),1,1]),max(sim$x,state.quant[-(1:burn),1,2],state.quant2[-(1:burn),1,2],state.quant3[-(1:burn),1,2])), xlab="", ylab=expression(x[t]), main = "", cex.lab=2)
mtext(substitute(paste(a['m0'],"=",a0,",",b['m0'],"=",a1,",",a['s0'],"=",a2,",",b['s0'],"=",a3),list(a0=round(prior$am0,2),a1=round(prior$bm0,2),a2=round(prior$as0,2),a3=round(prior$bs0,2))),side=3,cex=0.85)
lines(0:nt, state.quant[,1,1])
lines(0:nt, state.quant[,1,2])
lines(0:nt, state.quant2[,1,1],col=2)
lines(0:nt, state.quant2[,1,2],col=2)
lines(0:nt, state.quant3[,1,1],col=4)
lines(0:nt, state.quant3[,1,2],col=4)
