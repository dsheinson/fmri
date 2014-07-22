source("fmri_ar_functions.r")
source("pf_functions.r")

# Set data and graphics paths
dpath = "../data/"
gpath = "../graphs/"

# Load fitted data
load(paste(dpath,"craig_mle.rdata",sep=""))
attach(out)

# Region and model labels
regions = c("FP","IPS-left","IPS-right","PV","SV-left","SV-right")
reg.names = dimnames(mle)[[3]]
Nr = length(regions)
mods = dimnames(mle)[[4]]
Nm = length(mods)
Nv = dim(mle)[1]

# Switch model labels
M.ind = c(which(mods == "M101"),which(mods == "M011"))
mod.names = mods
mod.names[M.ind[1]] = mods[M.ind[2]]
mod.names[M.ind[2]] = mods[M.ind[1]]

# model labels
mlabels = rep(NA,length(mods))
for(i in 1:length(mods)) mlabels[i] = eval(bquote(expression(M[.(paste(strsplit(mod.names[i],"")[[1]][-1],sep="",collapse=""))])))

# Load data: 5 by 5 by 5 voxel clusters from 6 different brain regions
conv = read.csv(paste(dpath,"basis.csv",sep=""),header=FALSE)[,1]
nt = length(conv)
voxels = array(NA, c(125, nt+6, Nr))
for(i in 1:Nr) voxels[,,i] = as.matrix(read.csv(paste(dpath,reg.names[i],".csv",sep=""),header=F))

# Calculate marginal averages of mles
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
      to.clust = apply(mle.converge[,ind,j,k], 1, function(x) all(!
                                                                    is.na(x)))
      km = kmeans(mle.converge[to.clust,ind,j,k], centers=2)
      cluster[to.clust,j,k] = km$cluster
      centers[[k]][[j]] = km$centers
    }
  }
}

# pl parameters
n.run = 1
np = 500
sd.fac = 1
mle = FALSE
smooth = FALSE
mix = FALSE

# Load log marginal likelihoods
lmargliks = array(NA, c(Nv,Nr,length(mods)))
for(k in 1:length(mods))
{
  for(j in 1:Nr)
  {
    for(i in 1:Nv)
    {
      if(mods[k] != "M001")
      {
        loaded = try(load(paste(dpath,paste(i,j,mods[k],n.run,np,sd.fac,mle,smooth,mix,sep="-"),".rdata",sep="")),silent=TRUE)
        if(!(class(loaded) == "try-error")) lmargliks[i,j,k] = pf.out$lmarglik
      } else {
        # Calculate marginal likelihood for M001
        y = voxels[i,-(1:6),j]
        U = cbind(1,conv)
        if(mle)
        {
          fit = lm(y ~ U - 1)
          prior = list()
          prior$b0 = fit$coef
          prior$B0 = matrix(c((sd.fac^2)*1000,0,0,(sd.fac^2)*225),nr=2)
          prior$am0 = ig.mom(summary(fit)$sigma^2, (sd.fac^2)*500)[1]
          prior$bm0 = ig.mom(summary(fit)$sigma^2, (sd.fac^2)*500)[2]
        } else {
          m = ifelse(j < 5, 5, 3)
          prior = list()
          prior$b0 = centers[[k]][[j]][cluster[i,j,k],1:2]
          prior$B0 = matrix(c((sd.fac^2)*1000,0,0,(sd.fac^2)*225),nr=2)
          prior$am0 = ig.mom(centers[[k]][[j]][cluster[i,j,k],m], (sd.fac^2)*500)[1]
          prior$bm0 = ig.mom(centers[[k]][[j]][cluster[i,j,k],m], (sd.fac^2)*500)[2]
        }
        lmargliks[i,j,k] = br.lmarglik(y, U, prior$b0, solve(prior$B0), prior$am0, prior$bm0)    
      }
      print(c(i,j,k))
    }
  }
}

# Plot distribution of log likelihoods for each model across regions
pdf(file=paste(gpath,"craig_pl-loglik-",paste(n.run,np,sd.fac,mle,smooth,mix,sep="-"),".pdf",sep=""),width=15,height=10)
par(mfrow=c(2,3),mar=c(7,6,4,2)+0.1,mgp=c(4,1,0))
for(i in 1:Nr)
{
  ind = which(mods == "M001")
  if(length(ind) > 0)
  {
    dlmargliks = lmargliks[,i,-ind,drop=FALSE]
    md = rep(lmargliks[1,i,ind], 2)
  } else {
    dlmargliks = lmargliks
    md = c(Inf,-Inf)
  }
  dmax = max(apply(dlmargliks, 2:3, function(x) max(density(x,na.rm=T)$y,na.rm=TRUE)),na.rm=TRUE)
  bmax = max(md[2],apply(dlmargliks, 2:3, function(a) max(density(a,na.rm=T)$x,na.rm=TRUE)),na.rm=TRUE) 
  bmin = min(md[1],apply(dlmargliks, 2:3, function(a) min(density(a,na.rm=T)$x,na.rm=TRUE)),na.rm=TRUE)
  cols = 2:(length(mods)+1)
  if(i == 1) xlab = "Log marginal likelihood" else xlab=""
  if(i == 1) ylab = "Density" else ylab=""
  plot(density(lmargliks[,i,1],na.rm=T),lwd=2,col=cols[1],main=regions[i],xlab=xlab,ylab=ylab,xlim=c(bmin,bmax),ylim=c(0,dmax),cex.axis=1.5,cex.lab=1.75,cex.main=2.5)
  if(length(mods) > 1) for(k in 2:length(mods)) lines(density(lmargliks[,i,k],na.rm=T),lwd=2,col=cols[k])
  if(i == 1) legend("topright",legend=mlabels,lty=rep(1,length(mods)),lwd=rep(2,length(mods)),col=cols,cex=1.3)
}
dev.off()

# Compositional plots across regions
require(compositions)
pdf(file=paste(gpath,"craig_pl-comp-",paste(n.run,np,sd.fac,mle,smooth,mix,sep="-"),".pdf",sep=""),width=15,height=10)
par(mfrow=c(2,3),mar=c(7,6,4,2)+0.1,mgp=c(4,1,0))
pmargliks = array(NA, dim(lmargliks))
for(i in 1:Nr)
{
  # Compute posterior model probabilities
  for(j in 1:Nv)
  {
    pmargliks[j,i,] = postModProbs(lmargliks[j,i,], rep(1/length(mods),length(mods)))
  }
  
  # Compositional plot
  ac = acomp(pmargliks[,i,])
  for(j in 1:Nv) for(k in 1:length(mods)) if(is.BDL(ac[j,k])) ac[j,k] = 1e-6
  plot(ac, plotMissings=FALSE, labels=mlabels, lwd=2, cex.lab=2)
  title(regions[i],cex.main=2)
}
dev.off()

# Plot smoothed dynamic slopes with cluster identification and postmodprobs
mod = mods
for(i in 1:Nr)
{
  for(k in which(mod != "M001"))
  {
    if(mod[k] == "M011") dyn = 1 else dyn = 2
    for(j in 1:ceiling(Nv / 25))
    {
      js = (j-1)*25
      pdf(paste(gpath,"craig_state-",regions[i],"-",mod[k],"-",j,"-",mle,".pdf",sep=""),width=20,height=20)
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
          col = ifelse(cluster[l,i,k] == 1, 1, 6)
          plot(0:nt, x[l,,1,i,k], type="l", col=col, ylim=c(ymin,ymax), axes=axes[ind], xlab=xlab[ind], ylab=ylab[ind], main=substitute(paste(hat(phi)," = ",phihat,sep=""),list(phihat=round(mle.converge[l,3,i,k],3))), cex.lab=2.5, cex.main=2.5, cex.axis=1.4)
          lines(0:nt, x[l,,2,i,k], col=col, lty=2)
          lines(0:nt, x[l,,3,i,k], col=col, lty=2)
          abline(h=0, col='gray', lwd=2)
          box()
          pd = numeric(length(mods))
          cols = 1:length(mods) + 1
          for(m in 1:length(mods))
          {
            if(m == 1)
            {
              pd[m] = nt*pmargliks[l,i,m] 
              segments(pd[m],ymax+0.03*(ymax-ymin),pd[m],ymax-0.05*(ymax-ymin),lwd=2)
              polygon(c(pd[m],pd[m],0,0),c(ymax-0.05*(ymax-ymin),ymax+0.03*(ymax-ymin),ymax+0.03*(ymax-ymin),ymax-0.05*(ymax-ymin)),col=cols[m])
            } else {
              pd[m] = pd[m-1] + nt*pmargliks[l,i,m]
              segments(pd[m],ymax+0.03*(ymax-ymin),pd[m],ymax-0.05*(ymax-ymin),lwd=2)
              polygon(c(pd[m],pd[m],pd[m-1],pd[m-1]),c(ymax-0.05*(ymax-ymin),ymax+0.03*(ymax-ymin),ymax+0.03*(ymax-ymin),ymax-0.05*(ymax-ymin)),col=cols[m])
            }
          }
#          for(m in 1:length(mods)) if(m == 1) mtext(mlabels[m], side=3, at = pd[m]/2, line = -1.5) else mtext(mlabels[m], side=3, at = pd[m-1] + (pd[m]-pd[m-1])/2, line = -1.5)
          abline(h=ymax-0.05*(ymax-ymin),lwd=2)
          if(l == (js+1)) legend("bottomright",c("Cluster 1","Cluster 2","95% CI"), lty=c(1,1,2), cex=1.5, col=c(1,6,col))
          if(l == (js+6)) legend("bottomright",mlabels, fill=cols, cex=1.5)          
        } else {
          plot(0:5,0:5,main="",xlab="",ylab="",axes=FALSE,col='white')
          legend(c(1,4), paste(c("Error: ","Converge: "),c(error[l,i,k],converge[l,i,k]),sep=""),cex=1.5)
          box()
        }
      }
      dev.off()
    }
  }
}