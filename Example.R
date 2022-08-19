
# load simFTS.R. It contains code to simulate functional time series with trend
source('simFTS.R')

# simulated toy data ---- ------------------------------------------------------
# set parameters
N= 200 # sample size
m= 50  # num of values for each function
norm.op=0.5 # norm of the coeff operator for the stationary component

set.seed(1985)
Data <- sim.FTS(N=N,m=m,norm.op=norm.op, rangeval = c(0,1))

# plot data simulated
library(plot3D)

t <- seq(0,1,length.out = N) # time is normalized to (0,1)
Mm <- mesh(Data$s,t)
mx <- Mm$x; my <- Mm$y
surf3D( mx,my, Data$DataM, colkey = FALSE ,cex.main=2, xlab="s", ylab="t",
        theta = 45,phi=30,facets = FALSE,bty="g", ticktype = "detailed",
        scale=TRUE,expand=.7, col=rev(heat.colors(N)),main='Data (Yn)')

matplot(Data$s, Data$DataX, type = 'l', ylab = '', xlab = 's',
        main='Stationary Functional Time Series')

surf3D(mx,my, Data$TrendM, colkey = FALSE ,cex.main=2, xlab="s", ylab="t", 
        theta = 45,phi=30,facets = FALSE,bty="g", ticktype = "detailed", 
        main= 'Truth T(s,t)', scale=TRUE,expand=.7, col=rev(heat.colors(N)))


# estimate T(s,t) ---- ------------------------------------------------------
source('TrendEstimator.R')
source('smooth.bifda2.R')

hatT <-  Trend.hat(Data = Data$DataM, main='Trend Estimated', plot = FALSE)

par(mfrow=c(1,2))
surf3D(mx,my, Data$TrendM, colkey = FALSE ,cex.main=2, xlab="s", ylab="t", 
       theta = 45,phi=30,facets = FALSE,bty="g", ticktype = "detailed", 
       main= 'Truth T(s,t)', scale=TRUE,expand=.7, col=rev(heat.colors(N)))

surf3D(mx,my, hatT$TMat, colkey = FALSE ,cex.main=2, xlab="s", ylab="t", 
       theta = 45,phi=30,facets = FALSE,bty="g", ticktype = "detailed", 
       main= 'Estimated T(s,t)', scale=TRUE,expand=.7, col=rev(heat.colors(N)))
par(mfrow=c(1,1))

# trend forecast -----------------------------------------------------------
source('forecast.trend.R')
dim(Data$TrendM)

ff <-  Forecast_trend(data=Data$TrendM[,-c(199,200)], h=2)
# plot the true values and the forecast 
plot(Data$TrendM[,199], type = 'l', col=1)
lines(Data$TrendM[,200], col=2)
lines(ff[,1], col=1, lty = 2)
lines(ff[,2], col=3, lty = 2)


