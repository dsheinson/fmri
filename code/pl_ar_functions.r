# Functions to run particle learning for AR(1) DLM model

# theta = (beta, phi, sigma2s, sigma2m) (d + 3)-length vector
# y is a scalar, univariate observation
# x = (xt, t) (last element to to track time, since U and FF change over time)
# U is an nt by d matrix of known covariates
# FF is an nt-length vector of known state coefficients

dlpred.pl <- function(y, x, theta, U, FF)
{
  d = dim(U)[2]
  U = t(U[x[3]+1,])
  FF = FF[x[3]+1]
  mu = (U%*%theta[1:d]) + FF*theta[d+1]*x[1]
  tau = sqrt((FF^2)*theta[d+2] + theta[d+3])
  log.pred = dnorm(y, mu, tau, log=T)
  return(log.pred)
}

revo.pl <- function(y, x, theta, U, FF)
{
  d = dim(U)[2]
  U = t(U[x[3]+1,])
  FF = FF[x[3]+1]
  pred.var = theta[d+2]*FF^2 + theta[d+3]
  tau = sqrt((theta[d+3]*theta[d+2])/pred.var)
  mu = (theta[d+2] / pred.var)*(y - (U%*%theta[1:d]))*FF + (theta[d+3] / pred.var)*theta[d+1]*x[1]
  return(c(rnorm(1,mu,tau),x[1],x[3]+1))
}

rprior.pl <- function(prior)
{
  if(is.null(prior$am0)) sigma2m = 0 else sigma2m = 1 / rgamma(1, prior$am0, prior$bm0)  
  d = length(prior$b0)
  cl = chol(prior$B0)
  beta = prior$b0 + t(cl)%*%rnorm(d, 0, sqrt(sigma2m))
  sigma2s = 1 / rgamma(1, prior$as0, prior$bs0)
  cl = chol(prior$Phi0)
  phi = prior$phi0 + t(cl)%*%rnorm(1, 0, sqrt(sigma2s))
  suff.theta = c(prior$b0, as.numeric(prior$B0), ifelse(is.null(prior$am0), NA, prior$am0), 0, prior$phi0, prior$Phi0, prior$as0, 0, prior$bm0 + .5*t(prior$b0)%*%solve(prior$B0)%*%prior$b0, prior$bs0 + .5*(prior$phi0^2)/prior$Phi0)
  return(list(x=c(0,0,0), theta = c(beta,phi,sigma2s,sigma2m),suff.x=NA,suff.theta=suff.theta))
}

rmove.pl <- function(suff.theta, prior)
{
  d = length(prior$b0)
  bm = suff.theta[d+d^2+7]
  if(is.null(prior$am0)) sigma2m = 0 else sigma2m = 1 / rgamma(1, suff.theta[d + d^2 + 1], bm)
  beta = suff.theta[1:d] + t(chol(matrix(suff.theta[(d+1):(d + d^2)],nr=d)))%*%rnorm(d, 0, sqrt(sigma2m))
  bs = suff.theta[d+d^2+8]
  sigma2s = 1 / rgamma(1, suff.theta[d + d^2 + 5], bs)
  phi = suff.theta[d + d^2 + 3] + t(chol(matrix(suff.theta[d + d^2 + 4])))%*%rnorm(1, 0, sqrt(sigma2s))
  return(c(beta,phi,sigma2s,sigma2m))
}

smap.theta.pl <- function(suff.theta, y, x, U, FF, prior)
{
  d = dim(U)[2]
  U = t(U[x[3],])
  FF = FF[x[3]]
  tUU = t(U)%*%U
  tUy = t(U)%*%(y - FF*x[1])
  Binv.curr = solve(matrix(suff.theta[(d+1):(d+d^2)],nr=d))
  Binv = Binv.curr + tUU
  Binvb = Binv.curr%*%suff.theta[1:d] + tUy
  B = solve(Binv)
  b = B%*%Binvb
  suff.theta[1:d] = b
  suff.theta[(d+1):(d+d^2)] = as.numeric(B)
  suff.theta[d+d^2+1] = suff.theta[d+d^2+1] + .5
  suff.theta[d+d^2+2] = suff.theta[d+d^2+2] + (y - FF*x[1])^2  
  suff.theta[d+d^2+7] = .5*(suff.theta[d+d^2+2] + t(prior$b0)%*%solve(prior$B0)%*%prior$b0 - t(suff.theta[1:d])%*%solve(matrix(suff.theta[(d+1):(d+d^2)],nr=d))%*%suff.theta[1:d]) + prior$bm0

  Phiinv.curr = 1 / suff.theta[d+d^2+4]
  Phiinv = Phiinv.curr + x[2]^2
  Phiinvphi = Phiinv.curr*suff.theta[d+d^2+3] + x[1]*x[2]
  Phi = 1 / Phiinv
  phi = Phi*Phiinvphi
  suff.theta[d+d^2+3] = phi
  suff.theta[d+d^2+4] = Phi
  suff.theta[d+d^2+5] = suff.theta[d+d^2+5] + .5
  suff.theta[d+d^2+6] = suff.theta[d+d^2+6] + x[1]^2
  suff.theta[d+d^2+8] = .5*(suff.theta[d+d^2+6] + (prior$phi0^2)/prior$Phi0 - (suff.theta[d+d^2+3]^2)/suff.theta[d+d^2+4]) + prior$bs0

  return(suff.theta)
}

smap.state.pl <- function(suff.x, y, theta)
{
  return(NA)
}

dlevo.pl <- function(x.new, x.curr, theta)
{
  d = length(theta) - 3
  return(dnorm(x.new[1], theta[d+1]*x.curr[1], sqrt(theta[d+2]), log=T))
}