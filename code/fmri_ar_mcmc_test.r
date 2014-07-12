source("fmri_ar_mcmc.r")
require(dlm)

# Set data path
dpath = "../data/"

fmri_ar_mcmc_test <- function(mod.sim, mod.est, n.par, n.sim, n.chain, nsims, nburn, nthin, sd.fac=1, diffuse=TRUE, progress=FALSE, print.iter=FALSE)
{
  # Load simulated data
  load(paste(dpath,"fmri-ar-sim-",mod.sim,".rdata",sep=""))
  sim = sims[[1]][[n.par]][[n.sim]]
  nt = length(sim$y)
  d = length(sim$true.params$theta) - 3
  
  # Create FF depending on which model to fit
  if(mod.est == "M101")
  {
    FF = sim$true.params$U[,d]
  } else if(mod.est == "M011") {
    FF = rep(1,nt)
  } else if(mod.est == "M111") {
    FF = 1 + sim$true.params$U[,d]
  } else {
    stop("mod.est must be 'M101','M011', or 'M111'")
  }
  
  # Set prior
  abm = ig.mom(sim$true.params$theta[d+3], sd.fac*500)
  abs = ig.mom(sim$true.params$theta[d+2], sd.fac*500)
  prior = list(b0=sim$true.params$theta[1:d],B0=matrix(c(sd.fac*1000,0,0,sd.fac*225),nr=2), am0=abm[1], bm0=abm[2], phi0=sim$true.params$theta[d+1], Phi0=sd.fac*0.25, as0=abs[1], bs0=abs[2])

  # Use diffuse prior?
  if(diffuse)
  {
    prior$B0 = 1e6*diag(d)
    prior$Phi0 = 1e6
    prior$am0 = prior$bm0 = prior$as0 = prior$bs0 = 1e-6
  }
  
  # Run MCMC
  mcmc.details = list(n.sims = nsims, n.thin = nthin, n.burn = nburn)
  out = fmri.ar.mcmc(sim$y, list(U=sim$true.params$U, FF=FF), prior, mcmc.details=mcmc.details, progress=progress, print.iter=print.iter)
  
  # Calculate MLE
  try(fit <- dlmMLE(sim$y, c(sim$true.params$theta[d+1], log(sim$true.params$theta[d+2]), log(sim$true.params$theta[d+3])), function(par) build.ar1(par, sim$true.params$U, FF)), silent=T)
  if(!(class(fit) == "try-error"))
  {
    s=dlmSmooth(dlmFilter(sim$y, build.ar1(fit$par, sim$true.params$U, FF)))
    theta.mle = c(s$s[nt+1,1:d],fit$par[1],exp(fit$par[2:3]))
    x.mle = s$s[,d+1]
  } else {
    theta.mle = rep(NA,d+3)
    x.mle = rep(NA,nt+1)
  }
  
  # Allocate results in list and save
  out.est = list(out=out,theta.mle=theta.mle,x.mle=x.mle,prior=prior)
  file = paste(dpath,"fmri_ar_mcmc_test-",paste(mod.sim,mod.est,n.par,n.sim,n.chain,nsims,nburn,nthin,sd.fac,diffuse,sep="-"),".rdata",sep="")
  print(file)
  save(out.est, file = file)
}

require(plyr)
data1 = expand.grid(mod.sim="M101",mod.est="M101",n.par=560,n.sim=1:5,n.chain=1:4,nsims=25000,nburn=5000,nthin=5,diffuse=c(TRUE,FALSE))
data2 = expand.grid(mod.sim="M011",mod.est="M011",n.par=560,n.sim=1:5,n.chain=1:4,nsims=25000,nburn=5000,nthin=5,diffuse=c(TRUE,FALSE))
mydata = rbind(data1,data2)
require(doMC)
registerDoMC()
m_ply(mydata, fmri_ar_mcmc_test, .parallel = TRUE)