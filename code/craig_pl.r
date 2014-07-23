source("fmri_ar_functions.r")
source("pl_ar_functions.r")
source("pl.r")
source("pf_functions.r")
source("invgamma.r")
require(dlm)

# Set data path
dpath = "../data/"

# Load basis function
conv = read.csv(paste(dpath,"basis.csv",sep=""),header=FALSE)[,1]
U = cbind(1,conv)
regions = c("region-frontalpole-left","region-IPS-left","region-IPS-right","region-primaryvisual","region-secondaryvisual-left","region-secondaryvisual-right")
mods = c("M011","M101","M001")
nt = length(conv)
Nr = length(regions)
Nm = length(mods)
Nv = 125

# Load data: 5 by 5 by 5 voxel clusters from 6 different brain regions
voxels = array(NA, c(Nv, nt+6, Nr))
for(i in 1:Nr) voxels[,,i] = as.matrix(read.csv(paste(dpath,regions[i],".csv",sep=""),header=F))

# Calculate marginal averages of mles
load(paste(dpath,"craig_mle.rdata",sep=""))
mle.converge = array(NA, dim(out$mle))
for(j in 1:dim(mle.converge)[3])
{
  for(k in 1:dim(mle.converge)[4])
  {
    ind = which(out$converge[,j,k] == 0 & out$error[,j,k] == 0)
    mle.converge[ind,,j,k] = out$mle[ind,,j,k]
  }
}
avg.mle = apply(mle.converge,2:4,function(x) mean(x,na.rm=TRUE))

# Divide voxels into 2 clusters for secondary-visual cortex, based on MLEs
cluster = array(NA, c(Nv, Nr, Nm))
centers = list(); length(centers) = length(mods)
for(k in 1:Nm)
{
  centers[[k]] = list()
  length(centers[[k]]) = Nr
  for(j in 1:Nr)
  {
    if(j < 5)
    {
      cluster[,j,k] = 1
      centers[[k]][[j]] = matrix(avg.mle[,j,k],nr=1)
    } else {
      if(mods[k] == "M001") ind = c(1,2,5) else ind = 1:5
      to.clust = apply(mle.converge[,ind,j,k], 1, function(x) all(!is.na(x)))
      km = kmeans(mle.converge[to.clust,ind,j,k], centers=2)
      cluster[to.clust,j,k] = km$cluster
      centers[[k]][[j]] = km$centers
    }
  }
}

craig_pl <- function(n.vox, n.clust, mod.est, n.run, np, dyn, sd.fac=1, mle=FALSE, smooth=FALSE, alpha = 0.05, mix=FALSE, save.all=FALSE, progress = FALSE)
{
  rnorm(1)
  
  # Load data
  y = voxels[n.vox,-(1:6),n.clust]
  nt = length(y)
  d = dim(U)[2]
  
  # Create FF depending on which model to fit
  if(mod.est == "M101")
  {
    FF = U[,d]
    m = 2
  } else if(mod.est == "M011") {
    FF = rep(1,nt)
    m = 1
  } else {
    stop("mod.est must be 'M101' or 'M011'")
  }
  
  if(missing(dyn)) dyn = ifelse(mod.est=="M101", d, 1)
  
  # Set prior
  if(is.na(cluster[n.vox,n.clust,m])){
    abm = ig.mom(avg.mle[5,n.clust,m], (sd.fac^2)*500)
    abs = ig.mom(avg.mle[4,n.clust,m], (sd.fac^2)*500)
    prior = list(b0=avg.mle[1:2,n.clust,m],B0=matrix(c((sd.fac^2)*1000,0,0,(sd.fac^2)*225),nr=2), am0=abm[1], bm0=abm[2], phi0=avg.mle[3,n.clust,m], Phi0=(sd.fac^2)*0.25, as0=abs[1], bs0=abs[2])
  } else {
    abm = ig.mom(centers[[m]][[n.clust]][cluster[n.vox,n.clust,m],5], (sd.fac^2)*500)
    abs = ig.mom(centers[[m]][[n.clust]][cluster[n.vox,n.clust,m],4], (sd.fac^2)*500)
    prior = list(b0=centers[[m]][[n.clust]][cluster[n.vox,n.clust,m],1:2],B0=matrix(c((sd.fac^2)*1000,0,0,(sd.fac^2)*225),nr=2), am0=abm[1], bm0=abm[2], phi0=centers[[m]][[n.clust]][cluster[n.vox,n.clust,m],3], Phi0=(sd.fac^2)*0.25, as0=abs[1], bs0=abs[2])    
  }
  
  # Center priors at MLEs?
  fit0 = lm(y ~ U - 1)
  phi.init = cor(fit0$residuals[1:(nt-1)], fit0$residuals[2:nt])
  sigma2.init = summary(fit0)$sigma^2
  if(mle)
  {
    try(fit <- dlmMLE(y, c(phi.init, log(sigma2.init/2), log(sigma2.init/2)), function(par) build.ar1(par, U, FF), method="Nelder-Mead"), silent=T)
    if(!(class(fit) == "try-error"))
    {
      print(fit$convergence)
      s = dlmSmooth(dlmFilter(y, build.ar1(fit$par, U, FF)))
      abm = ig.mom(exp(fit$par[3]), (sd.fac^2)*500)
      abs = ig.mom(exp(fit$par[2]), (sd.fac^2)*500)
      prior$b0 = s$s[nt+1,1:d]; prior$phi0 = fit$par[1]
      prior$am0 = abm[1]; prior$bm0 = abm[2]
      prior$as0 = abs[1]; prior$bs0 = abs[2]
    } else {
      fit = fit0
      print("priors not centered at MLE")
    }
  } else {
    fit = c(fit0$coef,phi.init,sigma2.init)
  }
  
  # Run particle filter
  rprior <- function(j) rprior.pl(prior)
  dlpred <- function(y, x, suff.x, theta) dlpred.pl(y, x, theta, U, FF)
  revo <- function(y, x, suff.x, theta) revo.pl(y, x, theta, U, FF)
  rmove <- function(suff.theta) rmove.pl(suff.theta, prior)
  smap.theta <- function(suff.theta, y, x) smap.theta.pl(suff.theta, y, x, U, FF, prior)
  smap.state <- function(suff.x, y, theta) smap.state.pl(suff.x, y, theta)
  time = system.time(out <- pl(y, dlpred, revo, rprior, rmove, smap.theta, smap.state, np, progress, method="stratified", threshold=0.8, nonuniformity="ess", log=F))
  print(paste(n.vox,n.clust,mod.est,n.run,np,sd.fac,mle,smooth,mix,round(time[3],2),sep="-"))
    
  # Smooth states?
  if(smooth)
  {
    time = system.time(out.s <- pl.smooth(out$state, out$suff.theta, out$weight[,nt+1], dlevo.pl, rmove, progress, method="stratified", nonuniformity="ess", threshold=0.8, log=F))
    print(paste(n.vox,n.clust,mod.est,n.run,np,sd.fac,mle,smooth,mix,"smooth",round(time[3],2),sep="-"))
  } else {
    out.s = out
  }
  
  # Compute state quantiles
  state.quant.filt = pf.quantile(out$theta[dyn,,,drop=FALSE]+out$state[1,,,drop=FALSE], out$weight, function(x, j) x, c(alpha/2,.5,1-alpha/2))
  dimnames(state.quant.filt)[[3]] = c(alpha/2,.5,1-alpha/2)
  if(smooth)
  {
    tmp = resample(out$weight[,nt+1], method="stratified")
    r.theta = out$theta[dyn,tmp$indices,nt+1]
    state.quant.smooth = pf.quantile((r.theta + out.s$state)[1,,,drop=FALSE], matrix(1/np,nr=np,nc=nt+1), function(x, j) x, c(alpha/2,.5,1-alpha/2))
    dimnames(state.quant.smooth)[[3]] = c(alpha/2,.5,1-alpha/2)
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
    theta.quant = pf.mix.quantile(tq, out$weight, Fq, c(alpha/2,.5,1-alpha/2), int)
    print(paste(n.vox,n.clust,mod.est,n.run,np,sd.fac,mle,smooth,mix,"mix quant",sep="-"))
  } else {
    theta.quant = pf.quantile(out$theta, out$weight, function(x, j) x, c(alpha/2,.5,1-alpha/2))
  }
  dimnames(theta.quant)[[3]] = c(alpha/2,.5,1-alpha/2)
  
  # Compute log-marginal likelihood
  lmarglik = pf.lmarglik(out)
  
  # Save particles?
  if(!save.all) out = "didn't save particles"
  
  pf.out = list(out=out,fit=fit,lmarglik=lmarglik,state.quant.filt=state.quant.filt,state.quant.smooth=state.quant.smooth,theta.quant=theta.quant,prior=prior)
  save(pf.out, file=paste(dpath,paste(n.vox,n.clust,mod.est,n.run,np,sd.fac,mle,smooth,mix,sep="-"),".rdata",sep=""))
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(n.vox=1:125,n.clust=1:6,mod.est=c("M101"),n.run=1,np=5000)
m_ply(mydata, craig_pl, .parallel=T)
