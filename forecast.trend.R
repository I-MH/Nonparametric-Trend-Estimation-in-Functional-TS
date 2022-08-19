
Forecast_trend <- function(data, h=1){
  # this function is to forecast trend 
  # it uses Taylor expansion with linear term
  # args: 
  #   data = a matrix representing the evaluations of the trend, columns represents time
  #   h = an integer representing steps ahead to be forecast
  # values: a matrix representing the foecast values.   
  m <- dim(data)[1]
  N <- dim(data)[2]
  
  forecast_t <- matrix(0, m, h+1)
  forecast_t[,1] <- data[,N]    # this is an auxiliar column when h>1
  Deriv.data <- data[,N]-data[,N-1]
  forecast_t[,2] <- data[,N] + Deriv.data
  
  if(h>1){
    for (j in 2:h) {
      Deriv.data <- forecast_t[,j]-forecast_t[,j-1]
      forecast_t[,j+1] <- forecast_t[,j] + Deriv.data
    }
  }
  forecast_t <- forecast_t[,-1]
  return(forecast_t)
}


# see Example.R for an example