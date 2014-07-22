source("fmri_ar_functions.r")

########## Gibbs sampling functions ####################
# fmri.ar.mcmc runs an MCMC algorithm for a DLM where the state follows an AR(1) process
# Inputs:
# y is an nt length vector
# psi is a list with components U an nt by d matrix of beta covariates and FF an nt-length vector of state covariates
# prior is a list with components b0, B0 (mean and covariance on normal prior for beta), phi0, Phi0 (mean and covariance on normal prior for phi), am0, bm0 (shape and rate for inverse-gamma prior on sigma2m), and as0, bs0 (shape and rate for inverse-gamma prior on sigma2m)
# initial (optional) is a list with components x an nt + 1 length vector of states and theta a d+3 length vector of fixed parameters
# mcmc.details (optional) a list with scalar components n.sims (total number of MCMC iterations), n.burn (burn in iterations), and n.thin (save every n.thin-th iteration after burn-in)
# steps (optional) a character vector list which parameters to sample in the MCMC
# progress a boolean indicating if a progress bar should be displayed
# print.iter a boolean indicating if the iteration number should be printed (only if progress = FALSE)

fmri.ar.mcmc <- function(y, psi, prior, initial, mcmc.details, steps, progress=TRUE, print.iter=FALSE) {
  # Get dimensions of data and parameters
  nt = length(y)
  d = length(prior$b0)
  
  # Set up initial values
  if(missing(initial)) {
    initial = list()
    initial$x = rnorm(nt+1)
    initial$theta = c(rep(0,d+1),1,1)
  }
  x = initial$x
  theta = initial$theta
  
  # MCMC details
  if(missing(mcmc.details)) {
    n.thin <- 1   # save every n.thin(th) iteration
    n.sims <- 10
    n.burn <- 0
  } else {
    n.thin = mcmc.details$n.thin
    n.sims = mcmc.details$n.sims
    n.burn = mcmc.details$n.burn
  }
  if (missing(steps)) {
    steps=c('betasigma2m','phisigma2s','x')
  } 
  
  # save structures
  n.iter <- (n.sims - n.burn)%/%n.thin
  keep.theta <- matrix(NA, n.iter, d+3)
  keep.x  <- matrix(NA, nr=n.iter, nc=nt + 1)
  
  # Run mcmc
  if(progress) pb = txtProgressBar(0,n.sims,style=3)
  for (i in 1:n.sims) {
    if(print.iter & !progress & (i %% n.thin == 0))
    {
      print(i)
    } else if(progress) {
      setTxtProgressBar(pb,i)
    }
    
    if('betasigma2m' %in% steps) theta[c(1:d,d+3)] = sample.betasigma2m(y, x, theta, psi, prior)
    if('phisigma2s' %in% steps) theta[(d+1):(d+2)] = sample.phisigma2s(y, x, theta, psi, prior)
    if('x' %in% steps) x = sample.states(y, x, theta, psi, prior)
    
    # Only save every n.thin iteration
    if (ii <- save.iteration(i,n.burn,n.thin)) {
      keep.theta[ii,] = theta
      keep.x[ii,] = x
    }
  }
  
  return(list(theta = keep.theta, x=keep.x, mcmc.details=list(n.sims=n.sims,n.thin=n.thin,n.burn=n.burn), initial=initial))
}

# Functions to sample from full conditional distributions
# Inputs:
# y is a nt length vector of data
# x is a nt+1 length vector of states
# theta a d+1 length vector of fixed parameters
# psi is a list with components U a nt by d matrix of beta covariates and FF a nt length vector of state covariates
# prior is a list with components b0, B0 (mean and covariance on normal prior for beta), phi0, Phi0 (mean and covariance on normal prior for phi), am0, bm0 (shape and rate for inverse-gamma prior on sigma2m), and as0, bs0 (shape and rate for inverse-gamma prior on sigma2m)

sample.betasigma2m <- function(y, x, theta, psi, prior)
{
  nt = length(y)
  d = length(prior$b0)
  
  Uy = UU = yy = 0
  for(k in 1:nt)
  {
    Uy = Uy + matrix(psi$U[k,],nc=1)%*%(y[k] - psi$FF[k]*x[k+1])
    UU = UU + psi$U[k,]%*%t(psi$U[k,])
    yy = yy + (y[k] - psi$FF[k]*x[k+1])^2
  }
  
  B0.prec = solve(prior$B0)
  B0invb0 = B0.prec%*%prior$b0
  Bn.prec = UU + B0.prec
  Bn = solve(Bn.prec)
  bn = Bn%*%(Uy + B0invb0)

  amn = nt/2 + prior$am0
  bmn = .5*(yy + t(prior$b0)%*%B0invb0 - t(bn)%*%Bn.prec%*%bn) + prior$bm0
  
  sigma2m = 1 / rgamma(1, amn, bmn)
  beta = bn + t(chol(Bn))%*%rnorm(d,0,sqrt(sigma2m))
  return(c(beta,sigma2m))
}

sample.phisigma2s <- function(y, x, theta, psi, prior)
{
  tt = length(x)
  nt = tt - 1
  x0 = x[1:nt]
  x1 = x[2:tt]
  
  Phi0.prec = solve(prior$Phi0)
  Phi0invphi0 = Phi0.prec%*%prior$phi0
  Phin.prec = sum(x0^2) + Phi0.prec
  Phin = solve(Phin.prec)
  phin = Phin%*%(sum(x0*x1) +Phi0invphi0)
  
  asn = nt/2 + prior$as0
  bsn = .5*(sum(x1^2) + t(prior$phi0)%*%Phi0invphi0 - t(phin)%*%Phin.prec%*%phin) + prior$bs0
  
  sigma2s = 1 / rgamma(1, asn, bsn)
  phi = phin + t(chol(Phin))%*%rnorm(1,0,sqrt(sigma2s))
  return(c(phi,sigma2s))
}

sample.states <- function(y, x, theta, psi, prior)
{
  # Forward filtering backward sampling
  d = length(prior$b0)
  return(ffbs(y, psi$U, theta[1:d], psi$FF, theta[d+1], theta[d+3], theta[d+2], 0, 0))
}

###################
# Utility Functions
###################

save.iteration <- function(current.iter, n.burn, n.thin) {
  if (current.iter<n.burn | ((current.iter-n.burn)%%n.thin)!=0) {
    return(0)
  } else {
    return((current.iter-n.burn)/n.thin)
  }
}


ffbs <- function(y, U, beta, FF, G, V, W, m0, C0)
{
  nt = length(y)
  d = length(beta)

  # Initialize output arrays
  m = C = numeric(nt+1)
  f = Q = A = R = numeric(nt)
  m[1] = m0; C[1] = C0
  
  # Forward-filtering
  for(i in 1:nt)
  {    
    A[i] = G*m[i]; R[i] = (G^2)*C[i] + W
    f[i] = t(U[i,])%*%beta + FF[i]*A[i]; Q[i] = (FF[i]^2)*R[i] + V
    
    e = y[i] - f[i]; Qinv = 1 / Q[i]
    RFQinv = R[i]*FF[i]*Qinv
    m[i+1] = A[i] + RFQinv*e
    C[i+1] = R[i] - RFQinv*FF[i]*R[i]
  }
  
  # Backward-sampling
  x = numeric(nt+1)
  x[nt+1] = m[nt+1] + rnorm(1,0,sqrt(C[nt+1]))
  for(i in nt:1)
  {
    Rinv = 1 / R[i]
    CGRinv = C[i]*G*Rinv
    h = m[i] + CGRinv*(x[i+1] - A[i])
    H = C[i] - CGRinv*G*C[i]
    x[i] = h + rnorm(1,0,sqrt(H))
  }
  return(x)
}