#' WaveletAnn
#'
#' @param J 
#' @param wf 
#' @param boundary 
#' @param NForecast 
#' @param y 
#' @param ... 
#'
#' @return
#' @export
#' @importFrom forecast accuracy
#'
#' @examples
#' 
WaveletAnn<- function(y,wf="haar",J,boundary,NForecast,...)
{
  n1 <- length(y)-NForecast
  train_o <- y[1:n1]
  test_o <- y[-c(1:n1)]
  
  modwt <- waveslim::modwt(train_o, wf=wf, n.levels=J, boundary=boundary)
  WS <- do.call(cbind,modwt)
  
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of ANN model to the Wavelet Coef                  #
  #-----------------------------------------------------------#
  for(WVLevel in 1:ncol(WS))
  {
    ts <- NULL
    ts <- WS[,WVLevel]
    WaveletANNFit <- forecast::nnetar(y=as.ts(ts),...)

    WaveletANNPredict <- WaveletANNFit$fitted
    WaveletANNForecast <- forecast::forecast(WaveletANNFit,h=NForecast)
    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletANNPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletANNForecast$mean))
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

#' WaveletAnn with external regressors
#'
#' @param y 
#' @param xreg 
#' @param J 
#' @param wf 
#' @param boundary 
#' @param NForecast 
#' @param flag.comb 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' 
WaveletAnn.xreg<- function(y,xreg,wf="haar",J,boundary,NForecast,flag.comb="BG",...)
{
  n1 <- length(y)-NForecast
  train_o <- y[1:n1]
  test_o <- y[-c(1:n1)]
  
  xreg <- as.matrix(xreg)
  ndim<-ncol(xreg)
  
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of ANN model to the Wavelet Coef                  #
  #-----------------------------------------------------------#
  for(i in 1:ndim){

    WaveletANNFit <- forecast::nnetar(y=as.ts(train_o), xreg = xreg[1:n1,i],...)
    fit.xreg <- auto.arima(xreg[1:n1,i])
    WaveletANNPredict <- as.matrix(WaveletANNFit$fitted)
    
    #WaveletANNForecast <- forecast::forecast(WaveletANNFit,xreg = xreg[-c(1:n1),i], h=NForecast)
    WaveletANNForecast <- forecast::forecast(WaveletANNFit,xreg=forecast(fit.xreg,h=NForecast)$mean, h=NForecast)

    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletANNPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletANNForecast$mean))
  }
  
  # ###additional method
  # WaveletANNFit <- forecast::nnetar(y=as.ts(train_o), xreg = xreg[1:n1,],...)
  # 
  # WaveletANNPredict <- as.matrix(WaveletANNFit$fitted)
  # WaveletANNForecast <- forecast::forecast(WaveletANNFit,xreg = xreg[-c(1:n1),], h=NForecast)$mean

  #-----------------------------------------------------------#
  # Forecast combination                                      #
  #-----------------------------------------------------------#
  data.comb <- foreccomb(train_o, AllWaveletPrediction, newobs=test_o, newpreds=AllWaveletForecast)
  data.comb$Forecasts_Train <- ts(data.comb$Forecasts_Train)
  comb.fit <- do.call(paste0("comb_",flag.comb),list(data.comb))
  
  # wt <- matrix(comb.fit$Weights)
  # FinalPrediction <- comb.fit$Fitted
  # FinalForecast <- as.numeric(AllWaveletForecast%*%wt)
  # sum(abs(comb.fit$Forecasts_Test-FinalForecast)) 

  return(list(Accuracy.Train = comb.fit$Accuracy_Train,
              Accuracy.Test = comb.fit$Accuracy_Test,
              
              FinalForecast=comb.fit$Forecasts_Test,
              FinalPrediction=comb.fit$Fitted,
              
              Forecast = AllWaveletForecast,
              Comb = comb.fit))
}

