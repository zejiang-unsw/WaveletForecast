#' WaveletArima 
#'
#' @param y 
#' @param wf 
#' @param J 
#' @param boundary 
#' @param NForecast 
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


#' WaveletArima with external regressors
#'
#' @param y 
#' @param xreg 
#' @param NForecast 
#' @param method 
#' @param flag.comb 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' 
WaveletArima.reg<- function(y,xreg,NForecast,method=T,flag.comb="BG",...)
{
  n1 <- length(y)-NForecast
  train_o <- y[1:n1]
  test_o <- y[-c(1:n1)]
  
  ndim <- ncol(xreg)
  
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  if(method){
    #-----------------------------------------------------------#
    # Fitting of Arima model to predictor individually            #
    #-----------------------------------------------------------#
    for(i in 1:ndim){
      
      WaveletArimaFit <- forecast::auto.arima(y=as.ts(train_o), xreg = xreg[1:n1,i],...)
      fit.xreg <- auto.arima(xreg[1:n1,i]) # work with both seasonal and non-seasonal data
      WaveletArimaPredict <- as.matrix(WaveletArimaFit$fitted)
      
      #use the previous observation as instances or by forecast?
      #WaveletArimaForecast <- forecast::forecast(WaveletArimaFit,xreg=xreg[(n1-NForecast+1):n1,i], h=NForecast)
      WaveletArimaForecast <- forecast::forecast(WaveletArimaFit,xreg=forecast(fit.xreg,h=NForecast)$mean, h=NForecast)
      
      AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletArimaPredict)
      AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletArimaForecast$mean))
    }
    
    #pred_na <- c(1:max(apply(AllWaveletPrediction,2,function(x) max(which(is.na(x)))))) #NA from nnetar
    #data.comb <- foreccomb(train_o[-pred_na], AllWaveletPrediction[-pred_na,], newobs=test_o, newpreds=AllWaveletForecast)
    
    data.comb <- foreccomb(train_o, AllWaveletPrediction, newobs=test_o, newpreds=AllWaveletForecast)
    #data.comb$Forecasts_Train <- ts(data.comb$Forecasts_Train) #needed when there is NA
    comb.fit <- do.call(paste0("comb_",flag.comb),list(data.comb))
    
    
    Accuracy.Train = comb.fit$Accuracy_Train
    Accuracy.Test = comb.fit$Accuracy_Test
    FinalForecast = comb.fit$Forecasts_Test
    FinalPrediction=comb.fit$Fitted
    # wt <- matrix(comb.fit$Weights)
    # FinalForecast <- as.numeric(AllWaveletForecast%*%wt)
    # sum(abs(comb.fit$Forecasts_Test-FinalForecast)) 
    
  } else {
    #-----------------------------------------------------------#
    # Fitting of Arima model to predictor together                #
    #-----------------------------------------------------------#
    
    WaveletArimaFit <- forecast::auto.arima(y=as.ts(train_o), xreg = xreg[1:n1,],...)
    WaveletArimaPredict <- WaveletArimaFit$fitted
    
    if(T) { #use forecasts as instance
      fit.xreg <- lapply(1:ndim, function(i) auto.arima(xreg[1:n1,i])) # work with both seasonal and non-seasonal data
      xreg.n <- sapply(fit.xreg, function(ls) forecast(ls,h=NForecast)$mean)
      
      if(NForecast==1) WaveletArimaForecast <- forecast::forecast(WaveletArimaFit,xreg=t(xreg.n),h=NForecast)$mean 
      else WaveletArimaForecast <- forecast::forecast(WaveletArimaFit,xreg=xreg.n,h=NForecast)$mean
      
    } else { #use past observation as instance
      
      if(NForecast==1) WaveletArimaForecast <- forecast::forecast(WaveletArimaFit,xreg=t(xreg[(n1-NForecast+1):n1,]),h=NForecast)$mean
      else WaveletArimaForecast <- forecast::forecast(WaveletArimaFit,xreg=xreg[(n1-NForecast+1):n1,],h=NForecast)$mean
    }
    
    
    Accuracy.Train = accuracy(as.ts(WaveletArimaPredict), train_o)
    Accuracy.Test = accuracy(as.ts(WaveletArimaForecast), test_o)
    FinalForecast = WaveletArimaForecast
    FinalPrediction=WaveletArimaPredict
  }
  
  
  return(list(Accuracy.Train = Accuracy.Train,
              Accuracy.Test = Accuracy.Test,
              
              FinalForecast = FinalForecast,
              FinalPrediction=FinalPrediction,
              pic=pic))
}


