#' WaveletSS
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
WaveletSS<- function(y,wf="haar",J,boundary,NForecast,...)
{
  n1 <- length(y)-NForecast
  train_o <- y[1:n1]
  test_o <- y[-c(1:n1)]
  
  modwt <- lapply(waveslim::modwt(train_o, wf=wf, n.levels=J, boundary=boundary), function(ls) ls[1:n1])
  WS <- do.call(cbind,modwt)
  #sum(abs(train_o-rowSums(WS))); plot.ts(cbind(train_o,rowSums(WS)))
  
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of State space model to the Wavelet Coef          #
  #-----------------------------------------------------------#
  for(WVLevel in 1:ncol(WS))
  {
    #ts <- NULL
    ts <- ts(WS[,WVLevel], start=start(y), frequency=frequency(y))
    
    # state space model
    WaveletFit <- StructTS(ts, type = "trend")
    
    WaveletPredict <- fitted(WaveletFit)
    WaveletForecast <- forecast::forecast(WaveletFit,h=NForecast)
    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletForecast$mean))
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

#' WaveletSS with external regressors
#'
#' @param y 
#' @param xreg 
#' @param NForecast 
#' @param order 
#' @param freq 
#' @param method 
#' @param flag.comb 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
#' 
WaveletSS.xreg<- function(y,xreg,NForecast, order=1, freq=4, par, wt, method=F,flag.comb="BG",...)
{
  n1 <- length(y)-NForecast
  train_o <- y[1:n1]
  test_o <- y[-c(1:n1)]
  #stopifnot(frequency(y)==freq) 
  stopifnot(length(y)==nrow(as.matrix(xreg)))
  ndim <- ncol(xreg)
  
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  if(method){
    #-----------------------------------------------------------#
    # Fitting of SS model to predictor individually             #
    #-----------------------------------------------------------#
    ssm <- function(parm, X) {
      #dlmModPoly(order, dV=parm[1]) + dlmModSeas(freq) + dlmModReg(X=X, dV=parm[2], dW=c(parm[3],0))
      dlmModPoly(order, dV=parm[1]) + dlmModReg(X=X, dV=parm[2], dW=c(parm[3],0))
    }
    #par=rep(0,3)
    
    for(i in 1:ndim){
      
      WaveletFit <- dlmMLE(y=as.ts(train_o), parm=par, X=xreg[1:n1,i], build=ssm) 
      if(WaveletFit$convergence==0) print("converge") else print("Did not converge")
      
      #fit.xreg <- auto.arima(xreg[1:n1,i]) # work with both seasonal and non-seasonal data
      #modFit <- ssm(parm=WaveletFit$par, X=forecast(fit.xreg,h=NForecast)$mean)
      #modFit <- ssm(parm=WaveletFit$par, X=xreg[1:n1,i])
      modFit <- ssm(parm=WaveletFit$par, X=xreg[,i])
      
      WaveletFilter <- dlmFilter(c(train_o, rep(NA, NForecast)), modFit)
      WaveletPredict <- WaveletFilter$f[1:(n1)]
      WaveletForecast <- WaveletFilter$f[-c(1:(n1))]    
      
      AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletPredict)
      AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletForecast))
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
    # Fitting of SS model to predictor together                 #
    #-----------------------------------------------------------#
    ssm <- function(parm, X) {
      #dlmModPoly(order, dV=parm[1]) + dlmModSeas(freq) + dlmModReg(X=X, dV=parm[2], dW=c(parm[3], rep(0,ncol(X))))
      #dlmModPoly(order, dV=parm[1]) + dlmModReg(X=X, dV=parm[2], dW=c(parm[3], rep(1,ncol(X))))
      dlmModPoly(order, dV=parm[1]) + dlmModSeas(freq) + dlmModReg(X=X, dV=parm[2], dW=c(parm[3], wt))
      #dlmModPoly(order, dV=parm[1]) + dlmModReg(X=X, dV=parm[2], dW=parm[3:(3+ncol(X))])
    }
    #par=rep(0,3)
    
    WaveletFit <- dlmMLE(y=train_o, parm=par, X=xreg[1:n1,], build=ssm) 
    #if(WaveletFit$convergence==0) print("converge") else print("Did not converge")
    if(WaveletFit$convergence!=0) print("Did not converge")
    
    modFit <- ssm(parm=WaveletFit$par, X=xreg)
    #modFit <- ssm(parm=WaveletFit$par, X=xreg[1:n1,])
    
    WaveletFilter <- dlmFilter(c(train_o, rep(NA, NForecast)), modFit)
    WaveletPredict <- WaveletFilter$f[1:(n1)]
    WaveletForecast <- WaveletFilter$f[-c(1:(n1))]    
    
    Accuracy.Train = accuracy(WaveletPredict, train_o)
    Accuracy.Test =  accuracy(WaveletForecast, test_o)
    FinalForecast = WaveletForecast
    FinalPrediction=WaveletPredict
  }
  
  
  return(list(Accuracy.Train = Accuracy.Train,
              Accuracy.Test = Accuracy.Test,
              
              FinalForecast = FinalForecast,
              FinalPrediction=FinalPrediction))
}


