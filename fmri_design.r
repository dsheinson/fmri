source("fmri_sim.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to simulate design matrices and choose most stable
# N is number of designs to simulate before choosing the best
# n is length of experiment in seconds
# E is number of event types
# intercept is a boolean; should an intercept column be included in design matrix?
# multhrf is a boolean; are multiple hrfs included?
# ... other arguments passed to functions in fmri.sim such as TR length in seconds
optim.design <- function(N, n, E, intercept = TRUE, multhrf = FALSE, ...)
{
  hrf = boynton.hrf(n, ...)$hrf
  designs = list(); length(designs) = N
  for(j in 1:N)
  {
    boxcar = rapidEvent.boxcar(n, E, ...)
    basis1 = fmri.convolve(hrf,boxcar$boxcar[,1], ...)
    if(E > 1)
    {
      basis = matrix(NA, nr = length(basis1$t), nc = E)
      basis[,1] = basis1$basis
      for(i in 2:E) basis[,i] = fmri.convolve(hrf,boxcar$boxcar[,i],multhrf=multhrf, ...)$basis
    } else {basis = as.matrix(basis1$basis)}
    if(intercept) X = cbind(1,basis) else X = basis
    eig = eigen(t(X)%*%X)
    designs[[j]] = list(boxcar=boxcar, X=X, t=basis1$t, min.eig = min(eig$values), hrf=hrf)
    cat(j,eig$values,"\n")
  }
  
  # Find optimal design
  ind = which(sapply(designs, function(x) x$min.eig) == max(sapply(designs, function(x) x$min.eig)))
  
  # Plot optimal design
  pdf(file = paste(gpath,"fmri-design-",N,"-",n,"-",E,".pdf",sep=""))
  par(mfrow=c(3,1))
  ymax = max(designs[[ind]]$boxcar$boxcar)
  plot(designs[[ind]]$boxcar$t,designs[[ind]]$boxcar$boxcar[,1],ylim=c(0,ymax),type="l",xlab="",ylab="Boxcar")
  if(E > 1) for(i in 2:E) lines(designs[[ind]]$boxcar$t,designs[[ind]]$boxcar$boxcar[,i],col=i)
  hl = (length(hrf) %/% n)*30
  plot(designs[[ind]]$boxcar$t[1:hl],hrf[1:hl],type="l",xlab="",ylab="HRF")
  plot(basis1$t,designs[[ind]]$X[,intercept+1],type="l",ylim=c(min(designs[[ind]]$X),max(designs[[ind]]$X)),ylab="Convolution",xlab="Time (s)")
  if(E > 1) for(i in 2:E) lines(basis1$t,designs[[ind]]$X[,intercept+i],col=i)
  dev.off()
  
  # return optimal design
  return(designs[[ind]])
}

set.seed(65)
require(plyr)
mydata = expand.grid(N = 3, n = 500, E = 1)
fmri.design = mlply(mydata, optim.design)
save(fmri.design, file = paste(dpath,"fmri-design-1.rdata",sep=""))
