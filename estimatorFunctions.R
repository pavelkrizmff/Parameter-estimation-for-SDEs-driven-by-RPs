##########################################################################
### Estimating H and sigma: 
### returns estimates of Hurst parameter (H) and noise intensity (sigma) 
### based on the discretely sampled trajectory represented by the vector "x" 
##########################################################################


estimateHSigma <- function(x)
  {N <- length(x)-1
  quadrVar <- sum((x[2:(N+1)]-x[1:N])^2)
  step2 <- 2*(1:(N/2))
  xStep2 <- x[step2]
  quadrVarStep2 <- sum((xStep2[2:(N/2)]-xStep2[1:(N/2-1)])^2)
  estimH <- (log(quadrVarStep2/quadrVar)/(2*log(2)))+ 1/2
  estimSigma <- exp(-1/(2*log(2))*( log(N/2)*log(quadrVar) - log(N)*log(quadrVarStep2)) )
  c(estimH,estimSigma)
}



##########################################################################
### Estimating drift:
### function "estimateLambda" returns estimate of the drift intensity lambda based on the following inputs:
### "x" = numeric vector representing discretely sampled trajectory of the solution to SDE
### sigma = noise intensity (use true value or plug-in the estimate)
### H = Hurst parameter (use true value or plug-in the estimate)
### delta = deceleration parameter (use delta = 1 if sigma and H are true values and delta < 1 if sigma and H are estimated)
###  driftFction = drift function "f" for the underlying SDE
##########################################################################
myD <- function(H){(1/(H+1))*sqrt(2*(2*H-1)/H)}

myW <- function(xSub,h,H,sigma)  #xSub = the subsequence of x, h=time step, H=Hurst param., sigma = noise intensity
{
  n <- length(xSub)-1
  quadrVar <- sum((xSub[2:(n+1)]-xSub[1:n])^2)
  quadrVarNormalized <- sum(((xSub[2:(n+1)]-xSub[1:n])^2)/((sigma^2)*(h^(2*H)))-1) 
  ((1/h)^(1-H))/(4*myD(H))*h*quadrVarNormalized
}


estimateLambda <- function(x,sigma,H,delta,driftFction) #
{
  N <- length(x)-1   
  k <- floor(N/(N^delta))  #jump in nr. of observations
  h <- k/N  #jump in time units
  
  xSub1 <- x[seq(1,N/2 + 1,k)]
  xSub2 <- x[seq(N/2 + 1,N + 1,k)]
  U1 <- x[N/2 + 1]- x[1] - sigma*myW(xSub1,h,H,sigma)
  U2 <- x[N + 1] - x[N/2 + 1] - sigma*myW(xSub2,h,H,sigma)
  Z1 <- h*sum(sapply(xSub1,driftFction))
  Z2 <- h*sum(sapply(xSub2,driftFction))
  (Z1*U1+Z2*U2)/(Z1^2 + Z2^2)
  }
