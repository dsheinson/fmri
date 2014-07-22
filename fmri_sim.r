# boynton.hrf - returns the Boynton (1996) proposed hrf with parameters T0, the lag (seconds) between stimulus presentation and the initial rise in BOLD response, and shape parameters n and lambda
# d is an integer that specifies the number of hrf values to give per second (temporal resolution) 
# N is a scalar giving the length of the returned hrf (in seconds)
# t.zero - true or false, should hrf value of 0 at 0 seconds be included in returned list?
# returns a list with components t and hrf, vectors of the same length, with t[i] giving the time in seconds of hrf[i]
boynton.hrf <- function(N, d = 5, T0 = 0, n = 4, lambda = 2, t.zero = TRUE)
{
  Nt = ceiling(N*d)
  hrf = rep(NA, Nt)
  for(t in 1:Nt)
  {
    if(t/d > T0) hrf[t] = (((t/d - T0)^(n-1))*exp(-(t/d)/lambda)) / ((lambda^n)*factorial(n-1)) else hrf[t] = 0
  }
  t = (1:Nt)/d
  if(t.zero)
  {
    hrf = c(0, hrf)
    t = c(0, t)
  }
  return(list(t=t, hrf=hrf))
}

# rapidEvent.boxcar - generates neural activation boxcar functions with event times jittered according to a geometric distribution
# Inputs:
# n - numeric, length of experiment in seconds
# E - number of different events
# d - integer, number of boxcar elements given per second (temporal resolution)
# TR - numeric, number of seconds between TRs
# nmax - positive integer, maximum number of TRs to wait before presenting next stimulus
# duration - E-length vector with numeric elements number of seconds neural activiation lasts for each boxcar after presenting stimulus
# p - scalar between 0 and 1, parameter of geometric distribution from which to sample number of TRs to wait before presenting next stimulus
# boxcar.height - E-length vector of heights for each boxcar function
# Returns: n*d by E matrix of boxcar vectors of '0's and 'boxcar.height's
rapidEvent.boxcar <- function(n, E = 2, d=5, TR=2, nmax=10, duration=2, p=.5, boxcar.height=c(1,1), t.zero=TRUE)
{
  Nt = ceiling(n*d)
  boxcar = matrix(0, nr = Nt, nc = E)
  TR.curr = 0
  index = 0
  while(TR*(TR.curr+1) + duration < n)
  {
    g = rgeom(1,p)
    while(g > nmax) g = rgeom(1,p)
    if(TR*(TR.curr + g + 1) + duration < n)
    {
      e.ind = sample(1:E, 1)
      boxcar[(index + floor(TR*g*d) + 1):(index + floor(TR*g*d) + max(1,floor(duration*d))),e.ind] = boxcar.height[e.ind]
    }
    TR.curr = TR.curr + g + 1 + ceiling(duration/TR)
    index = floor(TR.curr*TR*d)
  }
  t = (1:Nt)/d
  if(t.zero)
  {
    boxcar = rbind(0, boxcar)
    t = c(0, t)
  }
  return(list(t=t, boxcar=boxcar))
}

# fmri.convolve - convolves a boxcar function of neural activiation with a haemodynamic response function (hrf) or combination of hrfs
# Inputs:
# hrf - either a numeric vector giving a single hrf or a matrix of multiple hrfs along the columns
# boxcar - a numeric vector of length equal to that of hrf (if vector) or number of rows of hrf (if matrix) giving neural activations
# TR - scalar, number of seconds per TR
# d - integer, number of elements in hrf and boxcar vectors per second
# multhrf - a boolean; if FALSE, a single hrf is assumed over the length of boxcar; if TRUE, different hrfs are assumed at intervals specified by 'cut'
# cut - a numeric vector of length equal to the number of columns of hrf - 1, or exactly the number of columns of hrf if time 0 is included in cut vector; ignored if multhrf = FALSE
# t.zero - true or false, is t = 0 seconds included at beginning of hrf and boxcar?
# Returns: vector of numerical convolution of boxcar with hrf (or multiple hrfs) of length = floor(length(boxcar) / tpr)
fmri.convolve <- function(hrf,boxcar,TR=2,d=5,multhrf=FALSE,cut,t.zero=TRUE)
{
  n = length(boxcar)
  hrf = as.matrix(hrf)
  if(dim(hrf)[1] != n) stop("hrf and boxcar must be same length")
  if(multhrf)
  {
    if(missing(cut)) cut = seq(0,n,len=dim(hrf)[2])
    if(!(is.numeric(cut) & length(cut) > 0)) stop("cut must be a numeric vector with at least one element")
    cut = sort(cut)
    if(!(0 %in% cut)) cut = c(0,cut)
    if(length(cut) != dim(hrf)[2]) stop("number of hrfs must equal number of cut points + 1")
  }
  basis = rep(0,n)
  pb = txtProgressBar(0,n,style=3)
  for(i in 1:n)
  {
    setTxtProgressBar(pb,i)
    for(j in 0:(i-1))
    {
      index = 1
      if(multhrf)
      {
        for(k in 2:length(cut)) 
        {
          if(j < cut[k] & j >= cut[k - 1])
          {
            index = k - 1
            break
          } else if(j >= cut[k] & k == length(cut)) {
            index = k
            break
          } else {
            # do nothing
          }
        }
      }
      basis[i] = basis[i] + boxcar[j+1]*hrf[i-j,index]
    }
  }
  tpr = TR*d
  ind = seq(t.zero + tpr,n,tpr)
  if(tpr %% 1 != 0)
  {
    basis.new = rep(NA, length(ind))
    for(i in 1:length(ind)) basis.new[i] = (basis[floor(ind[i])] + basis[ceiling(ind[i])])/2
  } else { basis.new = basis[ind] }
  return(list(t = (ind-t.zero)/d, basis = basis.new))
}