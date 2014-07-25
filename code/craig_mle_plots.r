source("fmri_ar_functions.r")
require(dlm)

# Set data and graphics paths
dpath = "../data/"
gpath = "../graphs/"

# Load fitted data
load(paste(dpath,"craig_mle.rdata",sep=""))
attach(out)

# Plot parameters
expr = expression(beta[0],beta[1],phi,sigma[s]^2,sigma[m]^2); Np = length(expr)
par = c("beta","phi","sigma","SNR")
regions = dimnames(mle)[[3]]; Nr = length(regions)
mods = dimnames(mle)[[4]]; Nm = length(mods)
nt = dim(x)[2] - 1
Nv = dim(mle)[1]

# Switch model labels
M.ind = c(which(mods == "M101"),which(mods == "M011"))
mod.names = mods; mod.names[M.ind[1]] = mods[M.ind[2]]; mod.names[M.ind[2]] = mods[M.ind[1]]

# model labels
mlabels = rep(NA,length(mods))
for(i in 1:length(mods)) mlabels[i] = eval(bquote(expression(M[.(paste(strsplit(mod.names[i],"")[[1]][-1],sep="",collapse=""))])))

# Calculate marginal averages of mles
mle.converge = array(NA, dim(out$mle))
for(j in 1:dim(mle.converge)[3])
{
  for(k in 1:dim(mle.converge)[4])
  {
    ind = which(out$converge[,j,k] == 0 & out$error[,j,k] == 0)
    if(length(ind) < Nv)
    {
      ind.err = which(out$error[,j,k] != 0)
      ind.con = which(out$converge[,j,k] != 0)
      print(c(paste("err:",ind.err),paste("con:",ind.con),j,k))
    }
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

# Plot MLEs across models and brain regions
for(k in 1:length(par))
{
  if((par[k] == "phi" | par[k] == "SNR") & "M001" %in% mods) mod = mods[-which(mods == "M001")] else mod = mods
  Nm = length(mod)
  pdf(paste(gpath,"craig_mle-",par[k],".pdf",sep=""),width=5*Nm,height=5*Nr)
  par(mfrow=c(Nr,Nm),mar=c(8,12,9,2)+0.1,mgp=c(7.5,1,0))
  ylab = c("FP","IPS-left","IPS-right","PV","SV-left","SV-right")
  if(length(mod) > 1) for(l in 2:length(mod)) ylab = cbind(ylab,"")
  for(j in 1:Nr)
  {
    for(i in 1:Nm)
    {
      ind.c = which(error[,j,i] == 0)# & converge[,j,i] == 0)
      ind = which(error[,j,i] == 0 & converge[,j,i] == 0)
      error.rate = round(mean(error[,j,i]), 2)
      conv.rate = round(mean(converge[ind.c,j,i] == 0), 2)
      if(par[k] == "beta")
      {
        xlab = expr[1]
        if(length(mod) > 1) for(l in 2:length(mod)) xlab = c(xlab,"")
        require(KernSmooth)
        b = 1.06*c(sd(mle[ind,1,j,i],na.rm=T),sd(mle[ind,2,j,i],na.rm=T))*Nv^(-.2)
        tmp = bkde2D(cbind(mle[ind,1,j,i],mle[ind,2,j,i]), b, range.x=list(x1=c(min(mle[ind,1,j,])-1.5*b[1],max(mle[ind,1,j,]))+1.5*b[1],x2=c(min(mle[ind,2,j,])-1.5*b[2],max(mle[ind,2,j,])+1.5*b[2])))
        if(j == 1)
        {
          image(tmp$x1,tmp$x2,-tmp$fhat, xlab=xlab[i], ylab=ylab[j,i], main=mlabels[i], cex.lab=4,cex.axis=1.9,cex.main=4)
          if(i == 1) mtext(expr[2], side=2, line=2.5, cex=3)
        } else {
          image(tmp$x1,tmp$x2,-tmp$fhat, xlab="", ylab=ylab[j,i], main="", cex.lab=4,cex.axis=1.9)
        }
        #if(mod[i] != "M001") mtext(paste("Error rate: ",error.rate,", Converge rate: ",conv.rate,sep=""),side=3)
        contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
        abline(h=0)
        points(centers[[i]][[j]][,1], centers[[i]][[j]][,2], pch='x', col=4,cex=3.5,lwd=3)
      } else if(par[k] == "phi") {
        xlab = expr[3]
        if(length(mod) > 1) for(l in 2:length(mod)) xlab = c(xlab,"")
        if(j == 1)
        {
          hist(mle[,3,j,i], xlab=xlab[i], ylab=ylab[j,i], main=mlabels[i], cex.lab=4, cex.axis=1.9, cex.main=4)
        } else {
          hist(mle[,3,j,i], xlab="", ylab=ylab[j,i], main="", cex.lab=4, cex.axis=1.9)
        }
        abline(v=centers[[i]][[j]][,3],col=4,lwd=4)
      } else if(par[k] == "sigma") {
        if(mod[i] != "M001")
        {
          xlab = expr[4]
          if(length(mod) > 1) for(l in 2:length(mod)) xlab = c(xlab,"")
          require(KernSmooth)
          b = 1.06*c(sd(mle[,4,j,i],na.rm=T),sd(mle[,5,j,i],na.rm=T))*Nv^(-.2)
          tmp = bkde2D(cbind(mle[ind,4,j,i],mle[ind,5,j,i]), b, range.x=list(x1=c(0,max(mle[ind,4,j,i]) + 1.5*b[1]),x2=c(0,max(mle[ind,5,j,i]) + 1.5*b[2])))
          if(j == 1)
          {
            image(tmp$x1,tmp$x2,-tmp$fhat, xlab=xlab[i], ylab=ylab[j,i], main=mlabels[i], cex.lab=4,cex.axis=1.9,cex.main=4)
            if(i == 1) mtext(expr[5], side=2, line=2.5, cex=3)
          } else {
            image(tmp$x1,tmp$x2,-tmp$fhat, xlab="", ylab=ylab[j,i], main="", cex.lab=4,cex.axis=1.9)
          }
          contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
          points(centers[[i]][[j]][,4], centers[[i]][[j]][,5], pch='x', col=4,cex=3.5,lwd=3)
        } else {
          if(length(mod) > 1) for(l in 2:length(mod)) xlab = c(xlab,"")
          if(j == 1)
          {
            hist(mle[,5,j,i], xlab=expr[5], ylab="", main=mlabels[i], cex.lab=4, cex.axis=1.9, cex.main=4)
          } else {
            hist(mle[,5,j,i], xlab="", ylab="", main="", cex.lab=4, cex.axis=1.9)
          }
          if(j < 5) abline(v=centers[[i]][[j]][,5],col=4,lwd=4) else abline(v=centers[[i]][[j]][,3],col=4,lwd=4)
        }
      } else if(par[k] == "SNR"){
        if(mod[i] != "M001")
        {
          if(mod[i] == "M011") xlab = expression(log(sigma[s]^2/sigma[m]^2)) else xlab = expression(sigma[s]^2/sigma[m]^2)
          if(mod[i] == "M011") snrs = log(mle[,4,j,i]/mle[,5,j,i]) else snrs = mle[,4,j,i]/mle[,5,j,i]
          if(j == 1)
          {
            hist(snrs, xlab=xlab, ylab=ylab[j,i], main=mlabels[i], cex.lab=4, cex.axis=1.9, cex.main=4)
          } else {
            hist(snrs, xlab="", ylab=ylab[j,i], main="", cex.lab=4, cex.axis=1.9)
          }
          abline(v=ifelse(mod[i] == "M011",0,1),lwd=4)
          if(mod[i] == "M011")
          {
            abline(v=log(centers[[i]][[j]][,4]/centers[[i]][[j]][,5]),col=4,lwd=4)
          } else {
            abline(v=centers[[i]][[j]][,4]/centers[[i]][[j]][,5],col=4,lwd=4)            
          }
        }
      } else {
        stop("par must contain only 'beta', 'phi', 'sigma', and/or 'SNR'")
      }
    }
  }
  dev.off()
}

# Plot smoothed dynamic slopes
mod = mods
for(i in 1:Nr)
{
  for(k in which(mod != "M001"))
  {
    if(mod[k] == "M011") dyn = 1 else dyn = 2
    for(j in 1:ceiling(Nv / 25))
    {
      js = (j-1)*25
      pdf(paste(gpath,"craig_state-",regions[i],"-",mod[k],"-",j,".pdf",sep=""),width=20,height=20)
      par(mfcol=c(5,5),mfg=c(5,1),mar=c(7,6,4,1)+0.1)
      ymin = min(x[(js+1):(js+25),,,i,k],na.rm=T)
      ymax = max(x[(js+1):(js+25),,,i,k],na.rm=T)
      xlab = c(expression(t),rep("",24))
      ylab = c(eval(bquote(expression(paste(beta[.(dyn-1)]," + ",x[t])))),rep("",24))
      axes = c(T,rep(F,24))
      for(l in (js+1):(js+25))
      {
        ind = ifelse(l == (js+25),25,l %% 25)
        if(error[l,i,k] == 0 & converge[l,i,k] == 0)
        {
          col = ifelse(cluster[l,i,k] == 1, 2, 4)
          plot(0:nt, x[l,,1,i,k], type="l", col=col, ylim=c(ymin,ymax), axes=axes[ind], xlab=xlab[ind], ylab=ylab[ind], main=substitute(paste(hat(phi)," = ",phihat,sep=""),list(phihat=round(mle[l,3,i,k],3))), cex.lab=2.5, cex.main=2.5, cex.axis=1.4)
          lines(0:nt, x[l,,2,i,k], col=col, lty=2)
          lines(0:nt, x[l,,3,i,k], col=col, lty=2)
          abline(h=0, col='gray', lwd=2)
          box()
          if(l == (js+1)) legend("topright",c(expression(x[t]),"95% CI"), lty=c(1,2), cex=1.75)
        } else {
          plot(0:5,0:5,main="",xlab="",ylab="",axes=FALSE,col='white')
          legend(c(1,4), paste(c("Error: ","Converge: "),c(error[l,i,k],converge[l,i,k]),sep=""),cex=1.75)
          box()
        }
      }
      dev.off()
    }
  }
}
