source("fmri_ar_functions.r")
source("pl_ar_functions.r")
source("pl.r")
source("pf_functions.r")
source("invgamma.r")
require(dlm)

# Set data path
dpath = "../data/"

fmri_pl <- function(mod.sim, mod.est, n.par, n.sim, n.run, np, dyn, sd.fac=1, mle=TRUE, smooth=FALSE, alpha = 0.05, mix=FALSE, save.all=FALSE, progress = FALSE)
{
  rnorm(1)
  
  # Load simulated data
  load(paste(dpath,"fmri-ar-sim-",mod.sim,".rdata",sep=""))
  sim = sims[[1]][[n.par]][[n.sim]]
  rm(sims)
  nt = length(sim$y)
  d = length(sim$true.params$theta) - 3
  
  # Create FF depending on which model to fit
  if(mod.est == "M101")
  {
    FF = sim$true.params$U[,d]
  } else if(mod.est == "M011") {
    FF = rep(1,nt)
  } else {
    stop("mod.est must be 'M101' or 'M011'")
  }
  
  if(missing(dyn)) dyn = ifelse(mod.est=="M101", d, 1)
  
  # Set prior
  abm = ig.mom(sim$true.params$theta[d+3], (sd.fac^2)*500)
  abs = ig.mom(sim$true.params$theta[d+2], (sd.fac^2)*500)
  prior = list(b0=sim$true.params$theta[1:d],B0=matrix(c((sd.fac^2)*1000,0,0,(sd.fac^2)*225),nr=2), am0=abm[1], bm0=abm[2], phi0=sim$true.params$theta[d+1], Phi0=(sd.fac^2)*0.25, as0=abs[1], bs0=abs[2])

  # Center priors at MLEs?
  fit0 = lm(sim$y ~ sim$true.params$U - 1)
  phi.init = cor(fit0$residuals[1:(nt-1)], fit0$residuals[2:nt])
  sigma2.init = summary(fit0)$sigma^2
  if(mle)
  {
    try(fit <- dlmMLE(sim$y, c(phi.init, log(sigma2.init/2), log(sigma2.init/2)), function(par) build.ar1(par, sim$true.params$U, FF), method="Nelder-Mead"), silent=T)
    if(!(class(fit) == "try-error"))
    {
      print(fit$convergence)
      s = dlmSmooth(dlmFilter(sim$y, build.ar1(fit$par, sim$true.params$U, FF)))
      abm = ig.mom(exp(fit$par[3]), (sd.fac^2)*500)
      abs = ig.mom(exp(fit$par[2]), (sd.fac^2)*500)
      prior$b0 = s$s[nt+1,1:d]; prior$phi0 = fit$par[1]
      prior$am0 = abm[1]; prior$bm0 = abm[2]
      prior$as0 = abs[1]; prior$bs0 = abs[2]
    } else {
      print("priors not centered at MLE")
      fit = c(fit0$coef,phi.init,sigma2.init)
    }
  } else {
    fit = c(fit0$coef,phi.init,sigma2.init)
  }
  
  # Run particle filter
  rprior <- function(j) rprior.pl(prior)
  dlpred <- function(y, x, suff.x, theta) dlpred.pl(y, x, theta, sim$true.params$U, FF)
  revo <- function(y, x, suff.x, theta) revo.pl(y, x, theta, sim$true.params$U, FF)
  rmove <- function(suff.theta) rmove.pl(suff.theta, prior)
  smap.theta <- function(suff.theta, y, x) smap.theta.pl(suff.theta, y, x, sim$true.params$U, FF, prior)
  smap.state <- function(suff.x, y, theta) smap.state.pl(suff.x, y, theta)
  time = system.time(out <- pl(sim$y, dlpred, revo, rprior, rmove, smap.theta, smap.state, np, progress, method="stratified", threshold=0.8, nonuniformity="ess", log=F))
  print(paste(mod.sim,mod.est,n.par,n.sim,n.run,np,sd.fac,mle,smooth,mix,round(time[3],2),sep="-"))
  
  # Smooth states?
  if(smooth)
  {
    time = system.time(out.s <- pl.smooth(out$state, out$suff.theta, out$weight[,nt+1], dlevo.pl, rmove, progress, method="stratified", nonuniformity="ess", threshold=0.8, log=F))
    print(paste(mod.sim,mod.est,n.par,n.sim,n.run,np,sd.fac,mle,smooth,mix,"smooth",round(time[3],2),sep="-"))
  } else {
    out.s = out
  }
  
  # Compute state quantiles
  state.quant.filt = pf.quantile(out$theta[dyn,,,drop=FALSE]+out$state[1,,,drop=FALSE], out$weight, function(x, j) x, c(alpha/2,1-alpha/2))
  dimnames(state.quant.filt)[[3]] = c(alpha/2,1-alpha/2)
  if(smooth)
  {
    tmp = resample(out$weight[,nt+1], method="stratified")
    r.theta = out$theta[dyn,tmp$indices,nt+1]
    state.quant.smooth = pf.quantile((r.theta + out.s$state)[1,,,drop=FALSE], matrix(1/np,nr=np,nc=nt+1), function(x, j) x, c(alpha/2,1-alpha/2))
    dimnames(state.quant.smooth)[[3]] = c(alpha/2,1-alpha/2)
  } else {
    state.quant.smooth = "No smoothing done"
  }
  
  # Compute quantiles for theta using sufficient statistics or particles?
  if(mix)
  {
    F1 = function(q, s, w) sum(w*pnorm(q, s[1,], sqrt(s[2,])))
    F2 = F3 = function(q, s, w) sum(w*pinvgamma(q,s[1,],s[2,]))
    Fq=c(rep(list(F1),d+1),list(F2,F3)) # To compute quantiles from mixture distribution
    int = list(c(-10000,10000),c(-10000,10000),c(-10000,10000),c(0,10000),c(0,10000))
    tq = list(out$suff.theta[c(1,3),,], out$suff.theta[c(2,6),,], out$suff.theta[9:10,,], out$suff.theta[c(11,14),,], out$suff.theta[c(7,13),,])
    theta.quant = pf.mix.quantile(tq, out$weight, Fq, c(alpha/2,1-alpha/2), int)
    print(paste(mod.sim,mod.est,n.par,n.sim,n.run,np,sd.fac,mle,smooth,mix,"mix quant",sep="-"))
  } else {
    theta.quant = pf.quantile(out$theta, out$weight, function(x, j) x, c(alpha/2,1-alpha/2))
  }
  dimnames(theta.quant)[[3]] = c(alpha/2,1-alpha/2)

  # Compute log-marginal likelihood
  lmarglik = pf.lmarglik(out)
  
  # Save particles?
  if(!save.all) out = "didn't save particles"
  
  pf.out = list(out=out,fit=fit,lmarglik=lmarglik,state.quant.filt=state.quant.filt,state.quant.smooth=state.quant.smooth,theta.quant=theta.quant,prior=prior)
  save(pf.out, file=paste(dpath,paste(mod.sim,mod.est,n.par,n.sim,n.run,np,sd.fac,mle,smooth,mix,sep="-"),".rdata",sep=""))
}

require(plyr)
require(doMC)
registerDoMC()
data1 = expand.grid(mod.sim="M101",mod.est="M101",n.run=1:5,n.par=560,n.sim=1,np=c(500,1000,5000,10000,20000),sd.fac=1,mle=FALSE)
data2 = expand.grid(mod.sim="M011",mod.est="M011",n.run=1:5,n.par=560,n.sim=1,np=c(500,1000,5000,10000,20000),sd.fac=1,mle=FALSE)
mydata = rbind(data1,data2)
m_ply(mydata, fmri_pl, .parallel=T)
