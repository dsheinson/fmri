# Function to approximate the log marginal likelihood using particle filtering given a list returned by rm_pf()
pf.lmarglik <- function(out, unnorm=F)
{
  nt = dim(out$weight)[2]-1
  if(unnorm) log.marglik = sum(log(apply(exp(out$n.weights),2,mean)*apply(exp(out$p.weights)*out$weight[,1:nt],2,sum))) else log.marglik = sum(log(apply(exp(out$increment)*out$weight[,1:nt],2,sum)))
  return(log.marglik)
}

# Function calculate the true log marginal likelihood of a bayesian regression model
br.lmarglik <- function(y, U, b0, B0.prec, am0, bm0)
{
  nt = length(y)
  Bn.prec = t(U)%*%U + B0.prec
  bn = solve(Bn.prec)%*%(t(U)%*%y + B0.prec%*%b0)
  amn = am0 + nt/2
  bmn = bm0 + .5*(sum(y^2) + t(b0)%*%B0.prec%*%b0 - t(bn)%*%Bn.prec%*%bn)
  loglik = -(nt/2)*log(2*pi) + 0.5*determinant(B0.prec)$mod - 0.5*determinant(Bn.prec)$mod + am0*log(bm0) - amn*log(bmn) + lgamma(amn) - lgamma(am0)
  return(loglik)
}

# Function to calculate posterior model probabilities given a set of log marginal likelihoods and prior model probabilities
postModProbs <- function(lmarglik, priorModProbs)
{
  stopifnot(length(lmarglik) == length(priorModProbs))
  stable = FALSE
  while(!stable)
  {
    denom = sum(exp(lmarglik)*priorModProbs,na.rm=T)
    if(denom == 0)
    {
      lmarglik = lmarglik + 100
    } else {
      stable = TRUE
    }
  }
  return(exp(lmarglik)*priorModProbs / sum(exp(lmarglik)*priorModProbs,na.rm=T))
}

# pf.quantile - function that calculates quantiles of filtered distributions; returns a 3-D array with dimensions number of observations, number of parameters and/or states, and number of quantiles desired
# Arguments:
# out - an 3-D dimensional array of particles with dimensions number of parameters and/or states, number of particles, and number of observations; this can be the state or theta component of the output list from functions bf, apf, or kd_pf
# wts - a matrix of wts with rows corresponding to particles and columns corresponding to observations; can be the weight component of the output list from functions bf, apf, or kd_pf; number of rows must match the 2nd dimension of out and number of columns must match the third dimension of out
# ftheta - a function that maps a vector from one scale to another, e.g. exp(theta) to map theta from log(scale) to original scale; can vary with parameter
# probs - a vector of quantiles to be calculated
# normwt - TRUE or FALSE corresponding to the normwt argument in wtd.quantile
pf.quantile = function(out, wts, ftheta, probs=.5, normwt=TRUE)
{
  require(Hmisc)
  nq = length(probs)
  np = dim(out)[1]
  tt = dim(out)[3]
  if(dim(wts)[1] != dim(out)[2] | dim(wts)[2] != tt) stop("dimensions of out and wts don't match")
  quantiles = array(NA,dim=c(tt,np,nq))
  for(i in 1:tt)
  {
    for(j in 1:np)
    {
      quantiles[i,j,] = wtd.quantile(ftheta(out[j,,i], j), wts[,i], normwt=normwt, probs=probs)
    }
  }
  return(quantiles)
}

# pf.mix.quantile - function that calculates quantiles of filtered distributions for states/parameters using mixture distributions, where sufficient statistics define the components of the mixture at each particle; returns a 3-D array with dimensions number of observations, number of params/states, and number of quantiles desired
# Arguments:
# out - a (#params)-length list of 3-D dimensional arrays of sufficient statistics for particles with dimensions number of sufficient statistics for parameter, number of particles, and number of observations; this can be the suff.state or suff.theta component of the output list from function pl
# wts - a matrix of wts with rows corresponding to particles and columns corresponding to observations; can be the weight component of the output list from function pl; number of rows must match the 2nd dimension of out and number of columns must match the third dimension of out
# probs - a vector of quantiles to be calculated
# F - a (#params)-length list of functions that return the cdf of the mixture distribution for a parameter given a quantile, an (dim(suff stat)) by (#particles) matrix of sufficient statistics and a (#particles)-length vector mixture weights
# int - 2-element vector giving interval range over which 'uniroot' searches for quantiles
pf.mix.quantile = function(out, wts, F, probs=c(.025,.975), int, print.iter=FALSE)
{
  npar = length(out)  
  nq = length(probs)
  tt = dim(wts)[2]
  quantiles = array(NA, c(tt,npar,nq))
  for(k in 1:npar)
  {
    stopifnot(dim(wts)[1] == dim(out[[k]])[2])
    for(i in 1:tt)
    {
      for(j in 1:nq)
      {
        G = function(x) F[[k]](x, out[[k]][,,i], wts[,i]) - probs[j]
        quantiles[i,k,j] = uniroot(G,int[[k]])$root
        if(print.iter) print(c(i,j,k))
      }
    }
  }
  return(quantiles)
}