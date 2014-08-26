source("pf_functions.r")
source("fmri_ar_functions.r")
require(coda)
require(plyr)

# Set data and graphics paths
dpath = "../data/"
gpath = "../graphs/"

# Traceplots of mcmc compared with MLEs
fmri_mcmc_trace <- function(mod.sim, mod.est, n.par, n.sim, dyn, mcmc)
{
  # Load simulated data
  load(paste(dpath,"fmri-ar-sim-",mod.sim,".rdata",sep=""))
  sim = sims[[1]][[n.par]][[n.sim]]
  rm(sims)
  nt = length(sim$y)
  
  # Find length of beta
  load(paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$n.chains[1],mcmc$nsims[1],mcmc$nburn[1],mcmc$nthin[1],mcmc$sd.fac,mcmc$diffuse,sep="-"),".rdata",sep=""))
  d = length(out.est$prior$b0)
  expr=c()
  for(i in 1:d) expr = c(expr,eval(bquote(expression(beta[.(i-1)]))))
  expr[d+1] = expression(phi); expr[d+2] = expression(sigma[s]^2); expr[d+3] = expression(sigma[m]^2)
  
  # Find plot mins and maxs for fixed parameters
  gmin = rep(Inf, d+3)
  gmax = rep(-Inf, d+3)
  for(i in 1:length(mcmc$n.chains))
  {
    load(paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$n.chains[i],mcmc$nsims[i],mcmc$nburn[i],mcmc$nthin[i],mcmc$sd.fac,mcmc$diffuse,sep="-"),".rdata",sep=""))
    gmin = apply(rbind(gmin,out.est$out$theta), 2, min)
    gmax = apply(rbind(gmax,out.est$out$theta), 2, max)
  }
  
  # Adjust for dynamic component
  if(missing(dyn)) dyn = ifelse(mod.est == "M101", d, 1)
  expr = expr[-dyn]
  gmin = gmin[-dyn]
  gmax = gmax[-dyn]
  par.ind = (1:(d+3))[-dyn]
  
  # Construct plots for fixed parameters
  n.keep = (mcmc$nsims - mcmc$nburn) %/% mcmc$nthin
  iter = list()
  for(i in 1:length(mcmc$n.chains)) iter[[i]] = (1:n.keep[i])*mcmc$nthin[i]
  ncol = ceiling(sqrt(d+2))
  nrow = ceiling((d+2) / ncol)
  pdf(file=paste(gpath,"fmri_mcmc_trace-theta-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$sd.fac,mcmc$diffuse,sep="-"),".pdf",sep=""),width=5*ncol,height=5*nrow)
  par(mfrow=c(nrow,ncol), mar = c(5,6,6,2) + 0.1)
  for(k in 1:length(par.ind))
  {
    if(k == 1) 
    {
      plot(iter[[1]],rep(sim$true.params$theta[par.ind[k]],n.keep[1]),type="l",col='gray',lwd=3,xlim=c(0,max(sapply(iter,max))),ylim=c(gmin[k],gmax[k]),xlab="Iteration",ylab=expr[k],cex.lab=2,cex.axis=1.2)
    } else {
      plot(iter[[1]],rep(sim$true.params$theta[k],n.keep[1]),type="l",col='gray',lwd=3,xlim=c(0,max(sapply(iter,max))),ylim=c(gmin[k],gmax[k]),xlab="",ylab=expr[k],cex.lab=2,cex.axis=1.2)
    }
    mtext(sim$true.params$theta[par.ind[k]], at=sim$true.params$theta[par.ind[k]],side=4)
    ess = 0
    for(j in 1:length(mcmc$n.chains))
    {
      load(paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$n.chains[j],mcmc$nsims[j],mcmc$nburn[j],mcmc$nthin[j],mcmc$sd.fac,mcmc$diffuse,sep="-"),".rdata",sep=""))
      lines(iter[[j]],out.est$out$theta[,par.ind[k]],col=j+1)
      ess = ess + effectiveSize(mcmc(out.est$out$theta[,par.ind[k]]))
    }
    mtext(paste("ESS:",round(ess,2)),side=3,cex=2)
    abline(h=c(sim$true.params$theta[par.ind[k]],out.est$theta.mle[par.ind[k]]),col=c('gray',1),lwd=c(3,3))
    if(k == 1) legend("topright",c(mcmc$n.chains,"MLE","Truth"),lty=rep(1,length(mcmc$n.chains)+2),lwd=c(rep(1,length(mcmc$n.chains)),3,3),col=c(2:(length(mcmc$n.chains)+1),1,'gray'),bg="white")
  }
  dev.off()
  
  # Find plot mins and maxs for states - 9 selected states
  tt = dim(out.est$out$x)[2]
  npts = floor(seq(1,tt,len=9))
  gmin = rep(Inf, 9)
  gmax = rep(-Inf, 9)
  for(i in 1:length(mcmc$n.chains))
  {
    load(paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$n.chains[i],mcmc$nsims[i],mcmc$nburn[i],mcmc$nthin[i],mcmc$sd.fac,mcmc$diffuse,sep="-"),".rdata",sep=""))
    gmin = apply(rbind(gmin,out.est$out$x[,npts]+out.est$out$theta[,dyn]), 2, min)
    gmax = apply(rbind(gmax,out.est$out$x[,npts]+out.est$out$theta[,dyn]), 2, max)
  }
  
  # Construct plots for states - 9 selected states
  pdf(file=paste(gpath,"fmri_mcmc_trace-states-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$sd.fac,mcmc$diffuse,sep="-"),".pdf",sep=""),width=15,height=15)
  par(mfrow=c(3,3), mar = c(5,6,6,2) + 0.1)
  for(k in 1:9)
  {
    if(k == 1)
    {
      plot(iter[[1]],rep(sim$x[npts[k]] + sim$true.params$theta[dyn],n.keep[1]),type="l",col='gray',lwd=3,xlim=c(0,max(sapply(iter,max))),ylim=c(gmin[k],gmax[k]),xlab="Iteration",ylab=eval(bquote(expression(beta[.(dyn-1)]))),cex.lab=2.25,cex.axis=1.2)
    } else {
      plot(iter[[1]],rep(sim$x[npts[k]],n.keep[1]),type="l",col='gray',lwd=3,xlim=c(0,max(sapply(iter,max))),ylim=c(gmin[k],gmax[k]),xlab="",ylab=eval(bquote(expression(paste(beta[.(dyn-1)]," + ",x[.(npts[k]-1)])))),cex.lab=2.25,cex.axis=1.2)
    }
    mtext(round(sim$x[k]+sim$true.params$theta[dyn],2), at=sim$x[npts[k]]+sim$true.params$theta[dyn],side=4)
    ess = 0
    for(j in 1:length(mcmc$n.chains))
    {
      load(paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$n.chains[j],mcmc$nsims[j],mcmc$nburn[j],mcmc$nthin[j],mcmc$sd.fac,mcmc$diffuse,sep="-"),".rdata",sep=""))
      lines(iter[[j]],out.est$out$x[,npts[k]]+out.est$out$theta[,dyn],col=j+1)
      ess = ess + effectiveSize(mcmc(out.est$out$x[,npts[k]]+out.est$out$theta[,dyn]))
    }
    mtext(paste("ESS:",round(ess,2)),side=3,cex=2)
    abline(h=c(sim$x[npts[k]]+sim$true.params$theta[dyn],out.est$x.mle[npts[k]]+out.est$theta.mle[dyn]),col=c('gray',1),lwd=c(3,3))
    if(k == 1) legend("topright",c(mcmc$n.chains,"MLE","Truth"),lty=rep(1,length(mcmc$n.chains)+2),lwd=c(rep(1,length(mcmc$n.chains)),3,3),col=c(2:(length(mcmc$n.chains)+1),1,'gray'),cex=1.75,bg="white")
  }
  dev.off()
}

mydata = data.frame(mod.sim=c("M101","M011"),mod.est=c("M101","M011"),n.par=560,n.sim=rep(1,2))
m_ply(mydata, function(mod.sim,mod.est,n.par,n.sim) fmri_mcmc_trace(mod.sim,mod.est,n.par,n.sim,mcmc=list(n.chains=1:3,nsims=rep(25000,3),nburn=rep(5000,3),nthin=rep(20,3),sd.fac=1,diffuse=TRUE)))
m_ply(mydata, function(mod.sim,mod.est,n.par,n.sim) fmri_mcmc_trace(mod.sim,mod.est,n.par,n.sim,mcmc=list(n.chains=1:3,nsims=rep(25000,3),nburn=rep(5000,3),nthin=rep(20,3),sd.fac=1,diffuse=FALSE)))

# Histograms of mcmc compared with MLEs
fmri_mcmc_hist <- function(mod.sim, mod.est, n.par, n.sim, dyn, mcmc)
{
  # Load simulated data
  load(paste(dpath,"fmri-ar-sim-",mod.sim,".rdata",sep=""))
  sim = sims[[1]][[n.par]][[n.sim]]
  rm(sims)
  nt = length(sim$y)
  
  # Find length of beta
  load(paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$n.chains[1],mcmc$nsims[1],mcmc$nburn[1],mcmc$nthin[1],mcmc$sd.fac,mcmc$diffuse,sep="-"),".rdata",sep=""))
  d = length(out.est$prior$b0)
  expr=c()
  for(i in 1:d) expr = c(expr,eval(bquote(expression(beta[.(i-1)]))))
  expr[d+1] = expression(phi); expr[d+2] = expression(sigma[s]^2); expr[d+3] = expression(sigma[m]^2)

  # Load mcmc samples
  mcmc.theta = matrix(NA, nr=0, nc=d+3)
  mcmc.x = matrix(NA, nr=0, nc=nt+1)
  for(i in 1:length(mcmc$n.chains))
  {
    load(paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$n.chains[i],mcmc$nsims[i],mcmc$nburn[i],mcmc$nthin[i],mcmc$sd.fac,mcmc$diffuse,sep="-"),".rdata",sep=""))
    mcmc.theta = rbind(mcmc.theta, out.est$out$theta)
    mcmc.x = rbind(mcmc.x, out.est$out$x)
  }
  
  # Adjust for dynamic component
  if(missing(dyn)) dyn = ifelse(mod.est == "M101", d, 1)
  expr = expr[-dyn]
  par.ind = (1:(d+3))[-dyn]
  
  # Construct plots for fixed parameters
  ncol = ceiling(sqrt(d+2))
  nrow = ceiling((d+2) / ncol)
  pdf(file=paste(gpath,"fmri_mcmc_hist-theta-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$sd.fac,mcmc$diffuse,sep="-"),".pdf",sep=""),width=5*ncol,height=5*nrow)
  par(mfrow=c(nrow,ncol), mar = c(5,6,6,2) + 0.1)
  for(k in 1:length(par.ind))
  {
    if(k == 1) 
    {
      hist(mcmc.theta[,par.ind[k]], xlab=expr[k], main="", cex.lab=2, cex.axis=1.25)
      legend("topright",c("Truth","MLE"),lty=c(1,1),lwd=c(3,3),col=c("gray",2))
    } else {
      hist(mcmc.theta[,par.ind[k]], xlab=expr[k], ylab = "", main="", cex.lab=2, cex.axis=1.25)
    }
    abline(v=c(sim$true.params$theta[par.ind[k]], out.est$theta.mle[par.ind[k]]),lwd=c(3,3),col=c('gray',2))
  }
  dev.off()
  
  # Construct plots for states - 9 selected states
  tt = dim(out.est$out$x)[2]
  npts = floor(seq(1,tt,len=9))
  pdf(file=paste(gpath,"fmri_mcmc_hist-states-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$sd.fac,mcmc$diffuse,sep="-"),".pdf",sep=""),width=15,height=nrow*5)
  par(mfrow=c(3,3), mar = c(5,6,6,2) + 0.1)
  for(k in 1:9)
  {
    if(k == 1) 
    {
      hist(mcmc.x[,npts[k]] + mcmc.theta[,dyn], xlab=eval(bquote(expression(beta[.(dyn-1)]))), main="", cex.lab=2.25,cex.axis=1.2)
      legend("topright",c("Truth","MLE"),lty=c(1,1),lwd=c(3,3),col=c("gray",2))
    } else {
      hist(mcmc.x[,npts[k]] + mcmc.theta[,dyn], xlab=eval(bquote(expression(paste(beta[.(dyn-1)]," + ",x[.(npts[k]-1)],sep="")))), ylab="", main="", cex.lab=2.25,cex.axis=1.2)
    }
    abline(v=c(sim$x[npts[k]]+sim$true.params$theta[dyn], out.est$x.mle[npts[k]]+out.est$theta.mle[dyn]),lwd=c(3,3),col=c('gray',2))
  }
  dev.off()
}

mydata = data.frame(mod.sim=c("M101","M011"),mod.est=c("M101","M011"),n.par=560,n.sim=rep(1,2))
m_ply(mydata, function(mod.sim,mod.est,n.par,n.sim) fmri_mcmc_hist(mod.sim,mod.est,n.par,n.sim,mcmc=list(n.chains=1:3,nsims=rep(25000,3),nburn=rep(5000,3),nthin=rep(20,3),sd.fac=1,diffuse=TRUE)))

# Plot quantiles for multiple runs of a given sim and model with increasing number of particles
fmri_pl_quant <- function(mod.sim, mod.est, n.par, n.sim, nruns, np, dyn, sd.fac=1, mle=TRUE, smooth=FALSE, mix=FALSE, mcmc, burn = 25)
{
  # Load simulated data
  load(paste(dpath,"fmri-ar-sim-",mod.sim,".rdata",sep=""))
  sim = sims[[1]][[n.par]][[n.sim]]
  rm(sims)
  nt = length(sim$y)
  
  # Find length of beta and which quantiles computed
  load(paste(dpath,paste(mod.sim,mod.est,n.par,n.sim,nruns[1],np[1],sd.fac,mle,smooth,mix,sep="-"),".rdata",sep=""))
  d = length(pf.out$prior$b0)
  if(smooth) pf.out$state.quant = pf.out$state.quant.smooth else pf.out$state.quant = pf.out$state.quant.filt
  probs = as.numeric(dimnames(pf.out$state.quant)[[3]])
  expr=c()
  for(i in 1:d) expr = c(expr,eval(bquote(expression(beta[.(i-1)]))))
  expr[d+1] = expression(phi); expr[d+2] = expression(sigma[s]^2); expr[d+3] = expression(sigma[m]^2)

  # Adjust for dynamic component
  if(missing(dyn)) dyn = ifelse(mod.est == "M101", d, 1)
  expr = expr[-dyn]
  par.ind = (1:(d+3))[-dyn]
  
  # Load mcmc and compute quantiles
  if(!missing(mcmc))
  {
    mcmc.theta = matrix(NA, nr=0, nc=d+3)
    mcmc.x = matrix(NA, nr=0, nc=nt+1)
    for(i in 1:length(mcmc$n.chains))
    {
      load(paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,mcmc$n.chains[i],mcmc$nsims[i],mcmc$nburn[i],mcmc$nthin[i],mcmc$sd.fac,mcmc$diffuse,sep="-"),".rdata",sep=""))
      mcmc.theta = rbind(mcmc.theta, out.est$out$theta)
      mcmc.x = rbind(mcmc.x, out.est$out$x)
    }
    mcmc.theta.quant = t(apply(mcmc.theta[,-dyn,drop=FALSE], 2, function(x) quantile(x, probs)))
    mcmc.x.quant = t(apply(mcmc.x, 2, function(x) quantile(x+mcmc.theta[,dyn], probs)))
  } else {
    mcmc.theta.quant=cbind(rep(Inf,d+2),rep(-Inf,d+2))
    mcmc.x.quant=cbind(rep(Inf,nt+1),rep(-Inf,nt+1))
  }
  
  # Find plot mins and maxs
  gmin = rep(Inf,d+3)
  gmax = rep(-Inf,d+3)
  for(j in 1:length(np))
  {
    for(i in nruns)
    {
      load(paste(dpath,paste(mod.sim,mod.est,n.par,n.sim,i,np[j],sd.fac,mle,smooth,mix,sep="-"),".rdata",sep=""))
      if(smooth) pf.out$state.quant = pf.out$state.quant.smooth else pf.out$state.quant = pf.out$state.quant.filt
      if(burn > 0) tmp = pf.out$state.quant[-(1:burn),1,] else tmp = pf.out$state.quant[,1,]
      gmin[d+3] = min(gmin[d+3], tmp[,1], sim$x+sim$true.params$theta[dyn], mcmc.x.quant[,1])
      gmax[d+3] = max(gmax[d+3], tmp[,2], sim$x+sim$true.params$theta[dyn], mcmc.x.quant[,2])
      for(k in 1:length(par.ind))
      {
        if(burn > 0) tmp = pf.out$theta.quant[-(1:burn),par.ind[k],] else tmp = pf.out$theta.quant[,par.ind[k],]
        gmin[k] = min(gmin[k], tmp[,1], sim$true.params$theta[par.ind[k]], mcmc.theta.quant[k,1])
        gmax[k] = max(gmax[k], tmp[,2], sim$true.params$theta[par.ind[k]], mcmc.theta.quant[k,2])
      }
    }
  }
  
  # Plot
  if(mod.sim == "M101") mod.sim.name = "M011"
  if(mod.sim == "M011") mod.sim.name = "M101"
  if(mod.est == "M101") mod.est.name = "M011"
  if(mod.est == "M011") mod.est.name = "M101"
  ncol = ceiling(sqrt(d+3))
  nrow = ceiling((d+3) / ncol)
  pdf(file=paste(gpath,"fmri_pl_quant-",paste(mod.sim,mod.est,n.par,n.sim,sd.fac,mle,smooth,mix,sep="-"),".pdf",sep=""),width=5*ncol,height=5*nrow)
  par(mfrow=c(nrow,ncol), mar = c(7,9,9,2) + 0.1,mgp=c(4,1,0))
  plot(0:nt,sim$x+sim$true.params$theta[dyn],type="l",col='gray',lwd=3,ylim=c(gmin[d+3],gmax[d+3]),xlab=expression(t),ylab=eval(bquote(expression(paste(beta[.(dyn-1)]," + ",x[t])))),cex.lab=4,cex.axis=1.9)
  mtext(substitute(paste(B0['1,1'],"=",a2,",",B0['2,2'],"=",a3,",",Phi[0],"=",a5),list(a2=pf.out$prior$B0[1,1],a3=pf.out$prior$B0[2,2],a5=pf.out$prior$Phi0)),side=3,cex=1.5)
  for(j in 1:length(np))
  {
    for(i in nruns)
    {
      load(paste(dpath,paste(mod.sim,mod.est,n.par,n.sim,i,np[j],sd.fac,mle,smooth,mix,sep="-"),".rdata",sep=""))
      if(smooth) pf.out$state.quant = pf.out$state.quant.smooth else pf.out$state.quant = pf.out$state.quant.filt
      lines(0:nt, pf.out$state.quant[,1,1], col=j+1)
      lines(0:nt, pf.out$state.quant[,1,2], col=j+1)
    }
  }
  if(!missing(mcmc))
  {
    if(smooth)
    {
      lines(0:nt,mcmc.x.quant[,1],lwd=2)
      lines(0:nt,mcmc.x.quant[,2],lwd=2)
    } else {
      points(c(nt,nt),mcmc.x.quant[nt+1,],lwd=3, pch='x',cex=3.5)
    }
  }
  for(k in 1:length(par.ind))
  {
    plot(0:nt,rep(sim$true.params$theta[par.ind[k]],nt+1),type="l",col='gray',lwd=3,ylim=c(gmin[k],gmax[k]),xlab="",ylab=expr[k],cex.lab=4,cex.axis=1.9)
    if(k==1)
    {
      title(eval(bquote(expression(paste("Sim: ",M[.(paste(strsplit(mod.sim.name,"")[[1]][2:4],sep="",collapse=""))],", Fit: ", M[.(paste(strsplit(mod.est.name,"")[[1]][2:4],sep="",collapse=""))])))), cex.main=4,line=5)
      mtext(eval(bquote(expression(paste(beta," = (",.(paste(sim$true.params$theta[1:d],sep="",collapse=",")),")",", ",phi," = ",.(sim$true.params$theta[d+1]),", ",sigma[s]^2," = ",.(sim$true.params$theta[d+2]),", ",sigma[m]^2," = ",.(sim$true.params$theta[d+3]))))),side=3,cex=1.5)
    } else if(k==ncol-1) {
      mtext(substitute(paste(a['m0'],"=",a0,",",b['m0'],"=",a1,",",a['s0'],"=",a2,",",b['s0'],"=",a3),list(a0=round(pf.out$prior$am0,2),a1=round(pf.out$prior$bm0,2),a2=round(pf.out$prior$as0,2),a3=round(pf.out$prior$bs0,2))),side=3,cex=1.5)
    }
    for(j in 1:length(np))
    {
      for(i in nruns)
      {
        load(paste(dpath,paste(mod.sim,mod.est,n.par,n.sim,i,np[j],sd.fac,mle,smooth,mix,sep="-"),".rdata",sep=""))
        lines(0:nt, pf.out$theta.quant[,par.ind[k],1], col=j+1)
        lines(0:nt, pf.out$theta.quant[,par.ind[k],2], col=j+1)
      }
    }
    if(!missing(mcmc)) points(c(nt,nt),mcmc.theta.quant[k,],lwd=3, pch='x',cex=3.5) 
    if(k == ncol-1) legend("bottomright",c(paste(np,"particles"),"MCMC","Truth"),pch=c(rep(NA,length(np)),'x',NA),lty=c(rep(1,length(np)),NA,1),lwd=c(rep(1,length(np)),3,3),col=c(2:(length(np)+1),1,'gray'),cex=1.75)
  }
  dev.off()
}

fmri_pl_quant("M101","M101",560,1,1:5,c(500,1000,5000,10000,20000),mle=FALSE,mcmc=list(n.chains=1:3,nsims=rep(25000,3),nburn=rep(5000,3),nthin=rep(20,3),sd.fac=1,diffuse=FALSE))
fmri_pl_quant("M011","M011",560,1,1:5,c(500,1000,5000,10000,20000),mle=FALSE,mcmc=list(n.chains=1:3,nsims=rep(25000,3),nburn=rep(5000,3),nthin=rep(20,3),sd.fac=1,diffuse=FALSE))

# Plot log marginal likelihood of each model against parameter values for data simulated from M101
# npars is a matrix of parameter values
fmri_pl_loglik <- function(mod.sim, nsims, np, mods, npars, wpars, nruns, sd.fac=1, mle=TRUE, smooth=FALSE, mix=FALSE)
{
  # Track whether models are distinguishable
  if("phi" %in% wpars)
  {
    phi.snr = list()
    length(phi.snr) = length(nsims)
  }
  
  # Switch model labels
  M.ind = c(which(mods == "M101"),which(mods == "M011"))
  mod.names = mods
  mod.names[M.ind[1]] = mods[M.ind[2]]
  mod.names[M.ind[2]] = mods[M.ind[1]]
  if(mod.sim == "M101") mod.sim.name = "M011"
  if(mod.sim == "M011") mod.sim.name = "M101"
  if(mod.sim == "M001") mod.sim.name = "M001"
  
  # model labels
  mlabels = rep(NA,length(mods))
  for(i in 1:length(mods)) mlabels[i] = eval(bquote(expression(M[.(paste(strsplit(mod.names[i],"")[[1]][-1],sep="",collapse=""))])))
  
  
  for(s in 1:length(nsims))
  {
    # Load simulated data and create matrix of parameters
    n.sim = nsims[s]
    load(paste(dpath,"fmri-ar-sim-",mod.sim,".rdata",sep=""))
    pars = sims$data[npars,which(colnames(sims$data) %in% c("phi","sigma2s","sigma2m"))]
    xexpr = expression(phi,sigma[s]^2,sigma[m]^2)
    
    # Index over each parameter
    for(k in 1:length(wpars))
    { 
      # Find unique parameter values
      upar = list(); length(upar) = length(wpars[-k])
      kind = which(colnames(pars) == wpars[k])
      wind = which(colnames(pars) %in% wpars[-k])
      for(j in 1:length(wpars[-k])) upar[[j]] = unique(pars[,wind[j]])
      
      # Find where models distinguishable?
      if(wpars[k] == "phi")
      {
        phi.snr[[s]] = array(NA, c(prod(sapply(upar, length)),2,length(sd.fac)))
        tracker = 1
      }
      
      # Set up plot panels
      ncol = ceiling(sqrt(prod(sapply(upar,length))))
      nrow = ceiling(prod(sapply(upar, length)) / ncol)
      pdf(file=paste(gpath,"fmri_pl_loglik-",paste(mod.sim,n.sim,np,wpars[k],npars[1],npars[length(npars)],paste(sd.fac,sep="",collapse="-"),sep="-"),".pdf",sep=""),width=4*ncol,height=4*nrow)
      par(mfrow=c(nrow,ncol),mar=c(8,7.5,5,2)+0.1,mgp=c(5,1,0))
      
      # Index over each set of unique fixed parameter values
      uind = rep(1, length(upar))
      done = FALSE
      first = TRUE
      while(!done)
      {
        # Load log marginal likelihoods
        opar = numeric(length(uind))
        for(i in 1:length(uind)) opar[i] = upar[[i]][uind[i]]
        ind = which(apply(as.matrix(pars[,wind]), 1, function(x) all(x == opar)))
        if(length(ind) > 0)
        {
          par = pars[ind,kind]
          lmargliks = array(NA, c(length(nruns),length(par),length(mods)+length(sd.fac)-1))
          for(l in 1:dim(lmargliks)[3])
          {
            for(j in 1:length(par))
            {
              for(i in 1:length(nruns))
              {
                if(l <= length(mods))
                {
                  if(mods[l] != "M001")
                  {
                    loaded = try(load(paste(dpath,paste(mod.sim,mods[l],npars[ind[j]],n.sim,nruns[i],np,sd.fac[1],mle,smooth,mix,sep="-"),".rdata",sep="")),silent=TRUE)
                    if(!(class(loaded) == "try-error")) lmargliks[i,j,l] = pf.out$lmarglik
                  } else {
                    # Calculate marginal likelihood for M001
                    y = sims[[1]][[npars[ind[j]]]][[n.sim]]$y
                    U = sims[[1]][[npars[ind[j]]]][[n.sim]]$true.params$U
                    nt = length(y)
                    d = dim(U)[2]
                    if(mle)
                    {
                      fit = lm(y ~ U - 1)
                      prior = list()
                      prior$b0 = fit$coef
                      prior$B0 = matrix(c((sd.fac[1]^2)*1000,0,0,(sd.fac[1]^2)*225),nr=2)
                      prior$am0 = ig.mom(summary(fit)$sigma^2, (sd.fac[1]^2)*500)[1]
                      prior$bm0 = ig.mom(summary(fit)$sigma^2, (sd.fac[1]^2)*500)[2]
                    } else {
                      prior = list()
                      prior$b0 = rep(0,d)
                      prior$B0 = 1e6*diag(d)
                      prior$am0 = 1e-6
                      prior$bm0 = 1e-6
                    }
                    lmargliks[i,j,l] = br.lmarglik(y, U, prior$b0, solve(prior$B0), prior$am0, prior$bm0)
                  }
                } else {
                  sind = which(mods == mod.sim)
                  if(length(sind) > 0)
                  {
                    if(mods[sind] != "M001")
                    {
                      loaded = try(load(paste(dpath,paste(mod.sim,mods[sind],npars[ind[j]],n.sim,nruns[i],np,sd.fac[l-length(mods)+1],mle,smooth,mix,sep="-"),".rdata",sep="")),silent=TRUE)
                      if(!(class(loaded) == "try-error")) lmargliks[i,j,l] = pf.out$lmarglik
                    } else {
                      # Calculate marginal likelihood for M001
                      y = sims[[1]][[npars[ind[j]]]][[n.sim]]$y
                      U = sims[[1]][[npars[ind[j]]]][[n.sim]]$true.params$U
                      nt = length(y)
                      d = dim(U)[2]
                      if(mle)
                      {
                        fit = lm(y ~ U - 1)
                        prior = list()
                        prior$b0 = fit$coef
                        prior$B0 = matrix(c((sd.fac[l-length(mods)+1]^2)*1000,0,0,(sd.fac[l-length(mods)+1]^2)*225),nr=2)
                        prior$am0 = ig.mom(summary(fit)$sigma^2, (sd.fac[l-length(mods)+1]^2)*500)[1]
                        prior$bm0 = ig.mom(summary(fit)$sigma^2, (sd.fac[l-length(mods)+1]^2)*500)[2]
                      } else {
                        prior = list()
                        prior$b0 = rep(0,d)
                        prior$B0 = 1e6*diag(d)
                        prior$am0 = 1e-6
                        prior$bm0 = 1e-6
                      }
                      lmargliks[i,j,l] = br.lmarglik(y, U, prior$b0, solve(prior$B0), prior$am0, prior$bm0)
                    }
                  }
                }
              }
            }
          }
          
          # Where do models become distinguishable?
          if(wpars[k] == "phi")
          {
            for(i in 1:length(sd.fac))
            {
              phi.snr[[s]][tracker,1,i] = pars[ind[1],2] / pars[ind[1],3]
              s.ind = ifelse(i == 1, which(mods == mod.sim), length(mods)+i-1)
              smargliks = lmargliks[,,c(s.ind,which(mods != mod.sim)),drop=FALSE]
              max.ind = apply(smargliks, 1:2, function(x) which(x == max(x)))
              w.max = apply(max.ind, 2, function(y) all(y == 1))
              largestF = which(!w.max)
              if(length(largestF) > 0) par.ind = max(largestF) + 1 else par.ind = 1
              if(par.ind <= length(par)) phi.snr[[s]][tracker,2,i] = par[par.ind]
            }
            tracker = tracker + 1
          }
          
          # Construct plots
          if(first)
          {
            xlab = xexpr[kind]
            ylab = "Log marginal likelihood"
            main = eval(bquote(expression(paste("Simulated from ",M[.(paste(strsplit(mod.sim.name,"")[[1]][2:4],sep="",collapse=""))]))))
            first = FALSE
            leg = TRUE
          } else {
            xlab=""
            ylab=""
            main = ""
            leg = FALSE
          }
          inf.ind = which(lmargliks == -Inf | lmargliks == Inf)
          if(length(inf.ind) > 0)
          {
            ymin = min(lmargliks[-inf.ind],na.rm=T)
            ymax = max(lmargliks[-inf.ind],na.rm=T)
          } else {
            ymin = min(lmargliks,na.rm=T)
            ymax = max(lmargliks,na.rm=T)
          }
          plot(par, lmargliks[1,,1], type="b", ylim=c(ymin,ymax), xlab=xlab, ylab=ylab, main=main, cex.main=2.25, cex.lab=2.25, cex.axis=1.25)
          if(wpars[k] == "phi")
          {
            expr = substitute(paste(sigma[s]^2," = ",a1,", ",sigma[m]^2," = ",a2,sep=""),list(a1=pars[ind[1],2],a2=pars[ind[1],3]))
          } else if(wpars[k] == "sigma2s"){
            expr = substitute(paste(phi," = ",a1,", ",sigma[m]^2," = ",a2,sep=""),list(a1=pars[ind[1],1],a2=pars[ind[1],3]))
          } else {
            expr = substitute(paste(phi," = ",a1,", ",sigma[s]^2," = ",a2,sep=""),list(a1=pars[ind[1],1],a2=pars[ind[1],2]))
          }
          mtext(expr,side=3,cex=1.1)
          if(length(nruns) > 1)
          {
            for(i in 2:length(nruns)) lines(par,lmargliks[i,,1],type="b")
          }
          if(dim(lmargliks)[3] > 1)
          {
            cols=1
            ltys=1
            for(l in 2:dim(lmargliks)[3])
            {
              for(i in 1:length(nruns))
              {
                if(l <= length(mods))
                {
                  col = 2*(l-1)
                  lty = 1
                } else {
                  col = cols[which(mods == mod.sim)]
                  lty = 1 + l - length(mods)
                }
                lines(par,lmargliks[i,,l],type="b",col=col,lty=lty)
              }
              cols = c(cols,col)
              ltys = c(ltys,lty)
            }
          }
          if(leg)
          {
            if(length(sd.fac) > 1)
            {
              k.names = rep(NA, length(sd.fac)-1)
              for(i in 2:length(sd.fac)) k.names[i-1] = eval(bquote(expression(paste(kappa,"=",.(sd.fac[i])))))
              leg.names = c(mlabels,k.names)
            } else {
              leg.names = mlabels
            }
            legend(ifelse(mod.sim=="M001","topright","bottomleft"),leg.names,lty=ltys,pch=rep(1,dim(lmargliks)[3]),col=cols,cex=1.75,bg="white")
          }
        }
        
        # Move index over parameters
        for(i in 1:length(uind))
        {
          if(uind[i] < length(upar[[i]]))
          {
            uind[i] = uind[i] + 1
            break
          } else {
            uind[i] = 1
          }
        }
        if(all(uind == 1)) done = TRUE
      }
      dev.off()
    }
  }
  
  # Plot SNR versus phi?
  if("phi" %in% wpars)
  {
    ncol = ceiling(sqrt(length(sd.fac)))
    nrow = ceiling(length(sd.fac) / ncol)
    pdf(file=paste(gpath,"fmri_pl_loglik-phiSNR-",paste(mod.sim,n.sim,np,npars[1],npars[length(npars)],paste(sd.fac,sep="",collapse="-"),sep="-"),".pdf",sep=""),width=5*ncol,height=5*nrow)
    par(mfrow=c(nrow,ncol),mar=c(7,6,4,2)+0.1,mgp=c(4,1,0))
    xmin = min(sapply(phi.snr, function(x) min(x[,1,],na.rm=T)))
    xmax = max(sapply(phi.snr, function(x) max(x[,1,],na.rm=T)))
    ymin = min(sapply(phi.snr, function(x) min(x[,2,],na.rm=T)))
    ymax = max(sapply(phi.snr, function(x) max(x[,2,],na.rm=T)))
    for(i in 1:length(sd.fac))
    {
      if(i == 1)
      {
        xlab = expression(sigma[s]^2/sigma[m]^2)
        ylab = expression(phi)
      } else {
        xlab=""
        ylab=""
      }
      plot(phi.snr[[1]][,1,i],phi.snr[[1]][,2,i],xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,main=substitute(paste(kappa,"=",aa),list(aa=sd.fac[i])),cex.lab=2.25,cex.main=2.25,cex.axis=1.2)
      if(length(phi.snr) > 1)
      {
        for(j in 2:length(phi.snr)) points(phi.snr[[j]][,1,i],phi.snr[[j]][,2,i])
      }
      if(i == 1)
      {
        #legend("topright",legend=nsims,lty=rep(1,length(nsims)),col=c(1,2*(2:length(nsims)-1)),title=eval(bquote(expression(paste("Sim from ",M[.(paste(strsplit(mod.sim.name,"")[[1]][2:4],sep="",collapse=""))],sep="")))),cex=1.2)
      }
    }
    dev.off()
  }
}

fmri_pl_loglik("M101", 1:3, 500, c("M101","M011","M001"), 496:594, c("phi","sigma2s"), 1:3, sd.fac=c(1,5,10,15))
fmri_pl_loglik("M101", 1:3, 500, c("M101","M011","M001"), 496:594, c("phi","sigma2s"), 1:3)
fmri_pl_loglik("M011", 1:3, 500, c("M101","M011","M001"), 496:594, c("phi","sigma2s"), 1:3, sd.fac=c(1,5,10,15))
fmri_pl_loglik("M011", 1:3, 500, c("M101","M011","M001"), 496:594, c("phi","sigma2s"), 1:3)
fmri_pl_loglik("M001", 1, 500, c("M101","M011","M001"), 1:9, c("sigma2s","sigma2m"), 1:3)
fmri_pl_loglik("M001", 1, 500, c("M101","M011","M001"), 1:9, c("sigma2s","sigma2m"), 1:3, sd.fac=c(1,5,10,15))

fmri_pl_loglik_dens <- function(mod.sim, n.par, n.sim, mods, np, nruns, sd.fac=1, mle=TRUE, smooth=FALSE, mix=FALSE, dmax, blim)
{
  # Load simulated data
  load(paste(dpath,"fmri-ar-sim-",mod.sim,".rdata",sep=""))
  sim = sims[[1]][[n.par]][[n.sim]]
  rm(sims)
  nt = length(sim$y)
  d = dim(sim$true.params$U)[2]
  
  # Switch model labels
  M.ind = c(which(mods == "M101"),which(mods == "M011"))
  mod.names = mods
  mod.names[M.ind[1]] = mods[M.ind[2]]
  mod.names[M.ind[2]] = mods[M.ind[1]]
  if(mod.sim == "M101") mod.sim.name = "M011"
  if(mod.sim == "M011") mod.sim.name = "M101"
  if(mod.sim == "M001") mod.sim.name = "M001"
  
  # model labels
  mlabels = rep(NA,length(mods))
  for(i in 1:length(mods)) mlabels[i] = eval(bquote(expression(M[.(paste(strsplit(mod.names[i],"")[[1]][-1],sep="",collapse=""))])))
  
  # Load log marginal likelihoods
  lmargliks = array(NA, c(length(nruns),length(mods),length(np)))
  for(k in 1:length(np))
  {
    for(j in 1:length(mods))
    {
      for(i in 1:length(nruns))
      {
        if(mods[j] != "M001")
        {
          loaded = try(load(paste(dpath,paste(mod.sim,mods[j],n.par,n.sim,nruns[i],np[k],sd.fac,mle,smooth,mix,sep="-"),".rdata",sep="")),silent=TRUE)
          if(!(class(loaded) == "try-error")) lmargliks[i,j,k] = pf.out$lmarglik
        } else {
          # Calculate marginal likelihood for M001
          y = sim$y
          U = sim$true.params$U
          if(mle)
          {
            fit = lm(y ~ U - 1)
            prior = list()
            prior$b0 = fit$coef
            prior$B0 = matrix(c((sd.fac^2)*1000,0,0,(sd.fac^2)*225),nr=2)
            prior$am0 = ig.mom(summary(fit)$sigma^2, (sd.fac^2)*500)[1]
            prior$bm0 = ig.mom(summary(fit)$sigma^2, (sd.fac^2)*500)[2]
          } else {
            prior = list()
            prior$b0 = rep(0,d)
            prior$B0 = 1e6*diag(d)
            prior$am0 = 1e-6
            prior$bm0 = 1e-6
          }
          lmargliks[i,j,k] = br.lmarglik(y, U, prior$b0, solve(prior$B0), prior$am0, prior$bm0)
        }
      }
    }
  }
  
  # Plot kernel density estimates of of log marginal likelihoods
  ind = which(mods == "M001")
  if(length(ind) > 0)
  {
    dlmargliks = lmargliks[,-ind,,drop=FALSE]
    md = rep(lmargliks[1,ind,1], 2)
  } else {
    dlmargliks = lmargliks
    md = c(Inf,-Inf)
  }
  if(missing(dmax)) dmax = max(apply(dlmargliks, 2:3, function(x) max(density(x)$y,na.rm=TRUE)),na.rm=TRUE)
  if(missing(blim))
  {
    bmax = max(md[2],apply(dlmargliks, 2:3, function(a) max(density(a)$x,na.rm=TRUE)),na.rm=TRUE) 
    bmin = min(md[1],apply(dlmargliks, 2:3, function(a) min(density(a)$x,na.rm=TRUE)),na.rm=TRUE)
  }
  cols = c(1,2,4)
  ncol = ceiling(sqrt(length(np)))
  nrow = ceiling(length(np) / ncol)
  pdf(file=paste(gpath,"fmri_pl_loglik-",paste(mod.sim,n.par,n.sim,sep="-"),".pdf",sep=""),width=5*ncol,height=5*nrow)
  par(mfrow=c(nrow,ncol),mar=c(5,6,5,2)+0.1,mgp=c(4,1,0))
  xlab = c("Log marginal likelihood",rep("",length(np)-1))
  ylab = c("Density",rep("",length(np)-1))
  for(i in 1:length(np))
  {
    if(mods[1] != "M001")
    {
      plot(density(lmargliks[,1,i]),axes=ifelse(i==1,T,F),lwd=2,col=cols[1],main=paste(np[i], " particles",sep=""),xlab=xlab[i],ylab=ylab[i],xlim=c(bmin,bmax),ylim=c(0,dmax),cex.axis=1.5,cex.lab=1.75,cex.main=2.5)
      box()
    } else {
      plot(0,0,lwd=2,axes=ifelse(i==1,T,F),col="white",main=paste(np[i], " particles",sep=""),xlab=xlab[i],ylab=ylab[i],xlim=c(bmin,bmax),ylim=c(0,dmax),cex.axis=1.5,cex.lab=1.75,cex.main=2.5)
      abline(v=lmargliks[1,1,i], col=cols[1])
      box()
    }
    if(length(mods) > 1) for(k in 2:length(mods)) if(mods[k] != "M001") lines(density(lmargliks[,k,i]),lwd=2,col=cols[k]) else abline(v=lmargliks[1,k,i], col=cols[k])
    if(i == 1) legend("topleft",legend=mlabels,lty=rep(1,length(mods)),lwd=rep(2,length(mods)),col=cols,cex=nrow)
    if(i == 1) mtext(eval(bquote(expression(paste(M[.(paste(strsplit(mod.sim.name,"")[[1]][2:4],sep="",collapse=""))],": ",beta," = (",.(paste(sim$true.params$theta[1:d],sep="",collapse=",")),")",", ",phi," = ",.(sim$true.params$theta[d+1]),", ",sigma[s]^2," = ",.(sim$true.params$theta[d+2]),", ",sigma[m]^2," = ",.(sim$true.params$theta[d+3]))))),side=3)
  }
  dev.off()
}

fmri_pl_loglik_dens("M101",498,1,c("M101","M011","M001"),c(500,1000,5000,10000),1:20)
fmri_pl_loglik_dens("M011",500,1,c("M101","M011","M001"),c(500,1000,5000,10000),1:20)
fmri_pl_loglik_dens("M001",6,1,c("M101","M011","M001"),c(500,1000,5000,10000),1:20)

fmri_pl_loglik_comp <- function(mod.sim, n.par, n.sim, mods, np, nruns, sd.fac=1, mle=TRUE, smooth=FALSE, mix=FALSE, dmax, blim)
{
  require(compositions)
  
  # Load simulated data
  load(paste(dpath,"fmri-ar-sim-",mod.sim,".rdata",sep=""))
  sim = sims[[1]][[n.par]][[n.sim]]
  rm(sims)
  nt = length(sim$y)
  d = dim(sim$true.params$U)[2]
  
  # Switch model labels
  M.ind = c(which(mods == "M101"),which(mods == "M011"))
  mod.names = mods
  mod.names[M.ind[1]] = mods[M.ind[2]]
  mod.names[M.ind[2]] = mods[M.ind[1]]
  if(mod.sim == "M101") mod.sim.name = "M011"
  if(mod.sim == "M011") mod.sim.name = "M101"
  if(mod.sim == "M001") mod.sim.name = "M001"
  
  # model labels
  mlabels = rep(NA,length(mods))
  for(i in 1:length(mods)) mlabels[i] = eval(bquote(expression(M[.(paste(strsplit(mod.names[i],"")[[1]][-1],sep="",collapse=""))])))
  
  # Load log marginal likelihoods
  lmargliks = array(NA, c(length(nruns),length(mods),length(np)))
  for(k in 1:length(np))
  {
    for(j in 1:length(mods))
    {
      for(i in 1:length(nruns))
      {
        if(mods[j] != "M001")
        {
          loaded = try(load(paste(dpath,paste(mod.sim,mods[j],n.par,n.sim,nruns[i],np[k],sd.fac,mle,smooth,mix,sep="-"),".rdata",sep="")),silent=TRUE)
          if(!(class(loaded) == "try-error")) lmargliks[i,j,k] = pf.out$lmarglik
        } else {
          # Calculate marginal likelihood for M001
          y = sim$y
          U = sim$true.params$U
          if(mle)
          {
            fit = lm(y ~ U - 1)
            prior = list()
            prior$b0 = fit$coef
            prior$B0 = matrix(c((sd.fac^2)*1000,0,0,(sd.fac^2)*225),nr=2)
            prior$am0 = ig.mom(summary(fit)$sigma^2, (sd.fac^2)*500)[1]
            prior$bm0 = ig.mom(summary(fit)$sigma^2, (sd.fac^2)*500)[2]
          } else {
            prior = list()
            prior$b0 = rep(0,d)
            prior$B0 = 1e6*diag(d)
            prior$am0 = 1e-6
            prior$bm0 = 1e-6
          }
          lmargliks[i,j,k] = br.lmarglik(y, U, prior$b0, solve(prior$B0), prior$am0, prior$bm0)
        }   
      }
    }
  }
  
  # Compute posterior model probabilities
  pmargliks = array(NA, dim(lmargliks))
  for(k in 1:length(np))
  {
    for(j in 1:length(mods))
    {
      for(i in 1:length(nruns))
      {
        pmargliks[i,,k] = postModProbs(lmargliks[i,,k], rep(1/length(mods),length(mods)))
      }
    }
  }
  
  # Compositional plots
  cols = 2:(length(np)+1)
  ncol = ceiling(sqrt(length(np)))
  nrow = ceiling(length(np) / ncol)
  pdf(file=paste(gpath,"fmri_pl_comp-",paste(mod.sim,n.par,n.sim,sep="-"),".pdf",sep=""),width=5*ncol,height=5*nrow)
  par(mfrow=c(nrow,ncol),cex.lab=2)
  for(i in 1:length(np))
  {
    plot(acomp(pmargliks[,,i]),labels=mlabels, lwd=2)
    title(paste(np[i]," particles",sep=""),cex.main=2)
    if(i == 1) mtext(eval(bquote(expression(paste(M[.(paste(strsplit(mod.sim.name,"")[[1]][2:4],sep="",collapse=""))],": ",beta," = (",.(paste(sim$true.params$theta[1:d],sep="",collapse=",")),")",", ",phi," = ",.(sim$true.params$theta[d+1]),", ",sigma[s]^2," = ",.(sim$true.params$theta[d+2]),", ",sigma[m]^2," = ",.(sim$true.params$theta[d+3]))))),side=3,cex=1.2,line=-1)
  }
  dev.off()
}

fmri_pl_loglik_comp("M101",498,1,c("M101","M011","M001"),c(500,1000,5000,10000),1:20)
fmri_pl_loglik_comp("M011",500,1,c("M101","M011","M001"),c(500,1000,5000,10000),1:20)
fmri_pl_loglik_comp("M001",6,1,c("M101","M011","M001"),c(500,1000,5000,10000),1:20)
