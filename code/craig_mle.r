source("fmri_ar_functions.r")
require(dlm)

# Set data and graphics paths
dpath = "../data/"
gpath = "../graphs/"

# Load basis function
conv = read.csv(paste(dpath,"basis.csv",sep=""),header=FALSE)[,1]
regions = c("region-frontalpole-left","region-IPS-left","region-IPS-right","region-primaryvisual","region-secondaryvisual-left","region-secondaryvisual-right")
mod = c("M011","M101","M001")
nt = length(conv)
Nr = length(regions)
nm = length(mod)

# Load data: 5 by 5 by 5 voxel clusters from 6 different brain regions
y = array(NA, c(125, nt+6, Nr))
for(i in 1:Nr) y[,,i] = as.matrix(read.csv(paste(dpath,regions[i],".csv",sep=""),header=F))

# Fit voxels to M101, M011, M010, M001 models using maximum likelihood estimation
Nv = 125
mle = array(NA, c(Nv, 5, Nr, nm)); dimnames(mle)[[2]] = expression(beta[0],beta[1],phi,sigma[s]^2,sigma[m]^2); dimnames(mle)[[3]] = regions; dimnames(mle)[[4]] = mod
x = array(NA, c(Nv, nt+1, 3, Nr, nm)); dimnames(mle)[[3]] = regions; dimnames(mle)[[4]] = mod
converge = array(0, c(Nv, Nr, nm)); dimnames(mle)[[3]] = regions; dimnames(mle)[[4]] = mod
error = array(0, c(Nv, Nr, nm)); dimnames(mle)[[3]] = regions; dimnames(mle)[[4]] = mod

# Fit voxels to M101, M011, M111, and M001 models
for(k in 1:nm)
{
  # Create U and f
  U = cbind(1,conv)
  if(mod[k] == "M101")
  {
    f = conv
    dyn = 2
  } else if(mod[k] == "M011") {
    f = rep(1,nt)
    dyn = 1
  } else {
    # do nothing
  }
  
  for(i in 1:Nr)
  {
    for(j in 1:Nv)
    {
      fit0 = lm(y[j,-(1:6),i] ~ conv)
      phi.init = cor(fit0$residuals[1:(nt-1)], fit0$residuals[2:nt])
      sigma2.init = summary(fit0)$sigma^2
      if(mod[k] == "M001")
      {
        mle[j,c(1,2,5),i,k] = c(fit0$coef, sigma2.init)
      } else {
        fit <- try(dlmMLE(y[j,-(1:6),i], c(phi.init,log(sigma2.init/2),log(sigma2.init/2)), function(par) build.ar1(par, U, f), method="Nelder-Mead"), silent=TRUE)
        if(class(fit) == "try-error")
        {
          error[j,i,k] = 1
        } else {
          converge[j,i,k] = fit$convergence
          s = dlmSmooth(dlmFilter(y[j,-(1:6),i], build.ar1(fit$par, U, f)))
          mle[j,,i,k] = c(s$s[nt+1,1:2],fit$par[1],exp(fit$par[2:3]))
          x[j,,1,i,k] = s$s[,3] + s$s[,dyn]
          var.x = sapply(dlmSvd2var(s$U.S, s$D.S), function(v) v[3,3])
          cov.xb = sapply(dlmSvd2var(s$U.S, s$D.S), function(v) v[3,dyn])
          var.b = sapply(dlmSvd2var(s$U.S, s$D.S), function(v) v[dyn,dyn])
          x[j,,2,i,k] = x[j,,1,i,k] - 1.96*sqrt(var.x + var.b + 2*cov.xb)
          x[j,,3,i,k] = x[j,,1,i,k] + 1.96*sqrt(var.x + var.b + 2*cov.xb)
        }
      }
      print(c(i,j,k,error[j,i,k],converge[j,i,k]))
    }
  }
}

# Save output
out = list(mle=mle,x=x,converge=converge,error=error)
save(out, file=paste(dpath,"craig_mle.rdata",sep=""))
