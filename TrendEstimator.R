
require(plot3D)
require(fda)
require(mgcv)

Trend.hat <- function(Data, nbases=10, nbaset=15,rangeval= c(0,1),lambdas=NULL,lambdat=NULL,
                      method=c('1D','2D'),plot=TRUE,...){
  # method= it determines how to compute lambda;
  # 1D= if lambda==NULL it uses a univ time series and fda package
  # 2D= if lambda==NULL it uses 'te' function in gam function
  # if lambda!=NULL and method 1D it uses fda package other case it uses gamm function
  method=match.arg(method)
  m <- dim(Data)[1]
  N <- dim(Data)[2]
  s <- seq(rangeval[1],rangeval[2],length.out = m)
  t <- seq(0,1,length.out = N)
  Mm <- mesh(s,t)
  mx <- Mm$x; my <- Mm$y
  #
  #MeanC <- rowMeans(Data)
  #Data <- Data - MeanC
  #
  bases <- create.bspline.basis(rangeval = rangeval, nbasis = nbases, norder = 4)
  baset <- create.bspline.basis(rangeval = c(0,1), nbasis = nbaset, norder = 4)
  
  if(method=='1D'){
    if( is.null(lambdas)){
      auxs <- rowMeans(Data)
      s0 <- s
      lambdas <- gam(auxs~te(s0, bs='cr'),method='REML')$sp
    }
    if( is.null(lambdat)){
      auxt <- colMeans(Data)
      t0 <- t
      lambdat <- gam(auxt~te(t0, bs='cr'),method='REML')$sp
    }
    sp.f=c(lambdas, lambdat)
    fdPs <- fdPar(bases, 2, lambda = lambdas)
    fdPt <- fdPar(baset, 2, lambda = lambdat)
    trend.hat <- smooth.bibasis2(s,t,Data, fdPs, fdPt)
    Trend.hat.Mat <- eval.bifd(s,t,trend.hat$bifdobj) 
  }else{
    sp <- c(lambdas,lambdat)
    trend.hat <- gam(c(Data)~te(c(mx),c(my), k=c(nbases,nbaset),bs="cr"),sp=sp,method='REML')
    Trend.hat.Mat <- matrix(trend.hat$fitted.values,m,N)
    sp.f <- trend.hat$sp
  }
  
    if(plot){
    Mm <- mesh(s,t)
    mx <- Mm$x; my <- Mm$y
    surf3D( mx,my, Trend.hat.Mat, colkey = FALSE ,cex.main=2, xlab="s", ylab="t",
           theta = 45,phi=30,facets = FALSE,bty="g", ticktype = "detailed"
           ,scale=TRUE,expand=.7, col=rev(heat.colors(dim(Data)[2])),...)
    }
  return(list(result=trend.hat, TMat=Trend.hat.Mat, sp=sp.f )  )
}

