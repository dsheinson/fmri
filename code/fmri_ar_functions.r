# Function to simulate nt + 1 univariate states (x_t) and nt univariate observations (y_t) from an AR(1) dlm
# y_t = U_t*beta + F_t*x_t + v_t, v_t iid N(0,sigma2m)
# x_t = phi*x_t-1 + w_t, w_t iid N(0,sigma2s)
# U is an nt by d matrix of regression covariates
# FF is an nt length vector of coefficients for the state
dlm.sim <- function(nt, beta, phi, sigma2s, sigma2m, U, FF, x0=0)
{
  y = numeric(nt)
  x = numeric(nt+1); x[1] = x0
  v = rnorm(nt, 0, sqrt(sigma2m))
  w = rnorm(nt, 0, sqrt(sigma2s))
  for(i in 1:nt)
  {
    x[i+1] = phi*x[i] + w[i]
    y[i] = t(U[i,])%*%beta + FF[i]*x[i+1] + v[i]
  }
  return(list(y=y,x=x,true.params=list(theta=c(beta,phi,sigma2s,sigma2m),U=U,FF=FF)))
}

# Function to build a 'dlm' object for a regression model with AR(1) error structure
# No constraint on autocorrelation parameter phi
# Assume univariate data
# par = c(phi, log(sigma2s), log(sigma2m))
# U is an nt by d matrix of regression covariates
# f is an nt length vector of coefficients for the state

build.ar1 <- function(par, U, f)
{
  require(dlm)
  d = dim(U)[2]
  nt = dim(U)[1]
  FF = matrix(rep(0,d+1),nr=1)
  JFF = matrix(1:(d+1),nr=1)
  X = cbind(U,f)
  V = exp(par[3])
  GG = diag(d+1); GG[d+1,d+1] = par[1]
  W = 0*diag(d+1); W[d+1,d+1] = exp(par[2])
  m0 = rep(0,d+1)
  C0 = 1e6*diag(d+1); C0[d+1,d+1] = 0
  mod = dlm(FF=FF,V=matrix(exp(par[3])),GG=GG,W=W,m0=m0,C0=C0,JFF=JFF,X=X)
  return(mod)
}

###################
# Utility functions
###################

ig.mom <- function(mu, sigma2)
{
  b = mu*((mu^2)/sigma2 + 1)
  a = b / mu + 1
  return(c(a, b))
}
