# codes to simulate functional time series with trend

TTrend <- function(s,t){
  #it defines a functional trend
  z <- 25*t*sin(2*pi*s) 
  return(z)
}


brow<-function(TT,h,N, sigma=1){ 
  # it simulates a Brownian motion
  #TT= [0,TT]
  m<-length(seq(0,TT,h))
  Aux= matrix(rnorm(n=N*(m-1) , sd=sqrt(h*sigma) ), m-1,N )
  Bt=rbind(rep(0,N),apply(Aux, 2, cumsum)  )
  return(Bt) }


data.ind <- function(N=100,m=50,TT=1){
  # it simulates a sequence of functional white noise
  # args:
  #   N: number of functions 
  #   m: num of points in a grid
  #   T: max of time to be simulated
  # values: a matrix

  data <- NULL
  h <- TT/m
  tt=seq(0,1,by=h)
  data <- brow(TT,h,N) 
  return(data) }


Tkernel <-  function(t,u,norm.op=0.5 ){
  # it defines a kernel for FAR(1) model
  # args: 
  #   t,u: arg values for the function
  #   norm.op: norm of the corresponding operator (should be <1)
  # values:
  #   z: value of the function at (t,u)
  C <- norm.op/0.746824
  z <- C * exp( -(t^2 + u^2) / 2 )
  return(z)}


sim.FTS <- function(N=100, m=50, norm.op=0.5, rangeval=NULL, burn=50){
  # it simulates a functional time series with trend: Y_n (s) = T_n(s) + X_n(s)
  # args:
  #   N: number of functions to be simulated (sample size)
  #   m: num of points for one function
  #   norm.op: norm of the corresponding operator (should be <1)
  #   rangeval: a numeric vector of length 2 defining the interval over which 
  #             the functional data object can be evaluated; default value 
  #             is if(is.null(breaks)) 0:1 else range(breaks).
  # values: a list 
  #   s: a vector of number where functions are evaluated
  #   DataM: a matrix containing the values of the functions 
  #   TrendM: a matrix woth the vlaeus of the trend used 
  #   DataX: a Matrix of the stationary time series data
  #   KerenlM: a matrix corresponding to the kernel of the coeff operator
  
  h=1/m
  
  if( is.null(rangeval) ){
    rangeval=c(0,1)
    s=seq(rangeval[1], rangeval[2] ,by=h)  
  }else{
    s=seq(rangeval[1], rangeval[2] ,by=h) 
  }
  
  KK <- outer(s,s, function(t,s) Tkernel(t,s, norm.op = norm.op))
  fw <- data.ind(N=N+burn+1,m=m,TT=rangeval[2]-rangeval[1] )
  
  X0=fw[,1]
  fw=fw[,-1]
  DatosFAR <- matrix(0,length(s),N+burn)
  DatosFAR[,1] <- X0 
  for(i in 2:(N+burn)){
    DatosFAR[,i] <- (t(KK) %*% DatosFAR[,(i-1)])*h + fw[,i]
  }
  DatosFAR <- DatosFAR[,(burn+1):(N+burn)] # functional TS
  # now we simulate the trend
  t <- seq(0,1, length.out = N)
  trendM <- outer(s,t,function(s,t) TTrend(s,t))
  DatosM <- DatosFAR + trendM
  return( list(s=s,DataM=DatosM, TrendM=trendM, DataX=DatosFAR,KernelM=KK)) 
}
