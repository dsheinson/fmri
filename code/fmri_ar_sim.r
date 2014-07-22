source("fmri_ar_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to simulate data from AR(1) dlm for M101, M011, or M111
# Assume coefficient of last column of U is dynamic
dlm_ar_sim <- function(design.lab, N, beta, phi, sigma2s, sigma2m, mod)
{
  # Load design
  load(paste(dpath,design.lab,sep=""))
  U = fmri.design[[1]]$X
  d = dim(U)[2]
  nt = dim(U)[1]
  stopifnot(length(beta) == d)
  
  # Create FF according to model
  if(mod == "M101")
  {
    FF = U[,d]
    mod.name = "M011"
  } else if(mod == "M011"){
    FF = rep(1,nt)
    if(sigma2m != 0) mod.name = "M101" else mod.name = "M100"
  } else if(mod == "M001"){
    FF = rep(0,nt)
    mod.name = "M001"
  } else {
    stop("mod must be 'M101', 'M011', or 'M001'")
  }
  
  # Simulate fMRI data
  mysims = list(); length(mysims) = N
  for(i in 1:N)
  {
    sim = dlm.sim(nt, beta, phi, sigma2s, sigma2m, U, FF)
    mysims[[i]] = sim
  }
  
  # Plot last simulation
  pdf(file=paste(gpath,"fmri-ar-sim-",paste(paste(beta,sep="",collapse="-"),phi*100,sigma2s,sigma2m,sep="-"),"-",mod,".pdf",sep=""))
  par(mfrow=c(3,1),mar=c(5,6,4,2))
  plot(1:nt,sim$y,type="l",xlab=expression(t),ylab=expression(y[t]),cex.lab=2)
  title(eval(bquote(expression(paste("Simulated from ",M[.(paste(strsplit(mod.name,"")[[1]][2:4],sep="",collapse=""))],sep="")))),cex.main=2)
  plot(1:nt,U[,d],type="l",xlab="",ylab=expression(conv[t]),cex.lab=2)
  title(eval(bquote(expression(paste(beta," = (",.(paste(sim$true.params$theta[1:d],sep="",collapse=",")),")",", ",phi," = ",.(sim$true.params$theta[d+1]),", ",sigma[s]^2," = ",.(sim$true.params$theta[d+2]),", ",sigma[m]^2," = ",.(sim$true.params$theta[d+3]))))),cex.main=2)
  plot(0:nt,sim$x,type="l",xlab="",ylab=expression(x[t]),main="",cex.lab=2)
  dev.off()
  
  # Save data
  cat(N, beta, phi, sigma2s, sigma2m, mod,"\n")
  return(mysims)
}

require(plyr)
M101data = expand.grid(design.lab="fmri-design-3-500-1.rdata", N = 20, phi=c(seq(0.1,0.9,.1),.95,.99), sigma2s = c(1:4,seq(5,25,5)), sigma2m = c(1:4,seq(5,25,5)), mod = "M101", stringsAsFactors=FALSE)
M101sim = mlply(M101data, function(design.lab,N,phi,sigma2s,sigma2m,mod) dlm_ar_sim(design.lab,N,c(750,15),phi,sigma2s,sigma2m,mod))
sims = list(sim = M101sim, data=M101data)
save(sims,file=paste(dpath,"fmri-ar-sim-M101.rdata",sep=""))

M011data = expand.grid(design.lab="fmri-design-3-500-1.rdata", N = 20, phi=c(seq(0.1,0.9,.1),.95,.99), sigma2s = c(1:4,seq(5,25,5)), sigma2m = c(1:4,seq(5,25,5)), mod = "M011", stringsAsFactors=FALSE)
M011sim = mlply(M011data, function(design.lab,N,phi,sigma2s,sigma2m,mod) dlm_ar_sim(design.lab,N,c(750,15),phi,sigma2s,sigma2m,mod))
sims = list(sim = M011sim, data=M011data)
save(sims,file=paste(dpath,"fmri-ar-sim-M011.rdata",sep=""))

M001data = expand.grid(design.lab="fmri-design-3-500-1.rdata", N = 20, phi=0, sigma2s = 0, sigma2m = c(1:4,seq(5,25,5)), mod = "M001", stringsAsFactors=FALSE)
M001sim = mlply(M001data, function(design.lab,N,phi,sigma2s,sigma2m,mod) dlm_ar_sim(design.lab,N,c(750,15),phi,sigma2s,sigma2m,mod))
sims = list(sim = M001sim, data=M001data)
save(sims,file=paste(dpath,"fmri-ar-sim-M001.rdata",sep=""))

M010data = expand.grid(design.lab="fmri-design-3-500-1.rdata", N = 1000, beta = c(0,1,2,3),phi=c(0.25,0.5,0.75,0.95), sigma2s = 15, sigma2m = 0, mod = "M011", stringsAsFactors=FALSE)
M010sim = mlply(M010data, function(design.lab,N,beta,phi,sigma2s,sigma2m,mod) dlm_ar_sim(design.lab,N,c(750,beta),phi,sigma2s,sigma2m,mod))
sims = list(sim = M010sim, data=M010data)
save(sims,file=paste(dpath,"fmri-ar-sim-M010.rdata",sep=""))
