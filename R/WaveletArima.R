#' WaveletArima 
#'
#' @param y 
#' @param wf 
#' @param J 
#' @param boundary 
#' @param NForecast 
#' @param max.p 
#' @param map.q 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' 
WaveletArima<- function(y,wf="haar",J,boundary,NForecast,...)
{
  n1 <- length(y)-NForecast
  train_o <- y[1:n1]
  test_o <- y[-c(1:n1)]
  
  modwt <- lapply(waveslim::modwt(train_o, wf=wf, n.levels=J, boundary=boundary), function(ls) ls[1:n1])
  WS <- do.call(cbind,modwt)
  #sum(abs(train_o-rowSums(WS))); plot.ts(cbind(train_o,rowSums(WS)))
  
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of ARIMA model to the Wavelet Coef                #
  #-----------------------------------------------------------#
  for(WVLevel in 1:ncol(WS))
  {
    ts <- NULL
    ts <- WS[,WVLevel]
    WaveletARMAFit <- forecast::auto.arima(x=as.ts(ts),...)
    
    WaveletARIMAPredict <- WaveletARMAFit$fitted
    WaveletARIMAForecast <- forecast::forecast(WaveletARMAFit,h=NForecast)
    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletARIMAPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletARIMAForecast$mean))
  }
  
  #-----------------------------------------------------------#
  # Forecast combination                                      #
  #-----------------------------------------------------------#
  FinalPrediction <- rowSums(AllWaveletPrediction,na.rm = T)
  FinalForecast <- rowSums(AllWaveletForecast,na.rm = T)
  
  Accuracy_Train <- accuracy(FinalPrediction, train_o)
  Accuracy_Test <- accuracy(FinalForecast, test_o)
  
  return(list(Accuracy.Train = Accuracy_Train,
              Accuracy.Test = Accuracy_Test,
    
              FinalForecast=FinalForecast,
              FinalPrediction=FinalPrediction))
}
