#' WaveletBMA
#'
#' @param y 
#' @param xreg 
#' @param NForecast 
#' @param transfo 
#' @param lambda 
#'
#' @return
#' @export
#' @importFrom forecast accuracy
#' 
#' @examples
#' 
WaveletBMA<- function(y, xreg, NForecast, transfo=c("log","sqrt"), lambda=0.5)
{ 
  n1 <- length(y)-NForecast
  
  ##!!!do we need transform for rainfall
  if(transfo=="log"){
    y <- log(y+lambda)
  } else if(transfo=="sqrt"){
    y <- sqrt(y+lambda)
  }
  if(sum(is.na(y))) warning(paste0("There are ",sum(is.na(y)), " missing values after transform!"))
  train_o <- y[1:n1] 
  test_o <- y[-c(1:n1)]
  
  xreg <- as.matrix(xreg) 
  
  FinalForecast <- rep(NA,NForecast); AllPrediction <- NULL
  #-----------------------------------------------------------#
  # Fitting of ANN model to the Wavelet Coef                  #
  #-----------------------------------------------------------#
  for(i in 1:NForecast){
    tmp <- rep(NA, n1)
    
    WaveletBMAFit <- bicreg(x=xreg[1:(n1-i),], y=train_o[-c(1:i)])
    #fit.xreg <- lapply(1:ncol(xreg), function(i) auto.arima(xreg[1:n1,i]))
    colnames(xreg) = WaveletBMAFit$input.names
    WaveletBMAPredict <- predict(WaveletBMAFit, newdata=xreg[1:(n1-i),])
    
    # xreg.n <- sapply(fit.xreg, function(x) forecast(x,h=NForecast)$mean)
    # colnames(xreg.n) = colnames(xreg)
    WaveletBMAForecast <- predict(WaveletBMAFit,newdata=data.frame(t(xreg[n1,])))
    
    tmp[-c(1:i)] <- WaveletBMAPredict$mean
    AllPrediction <- cbind(AllPrediction, tmp)

    FinalForecast[i] <- WaveletBMAForecast$mean 
  }
  

  if(transfo=="log"){
    FinalForecast <- exp(FinalForecast)-lambda
    FinalPrediction <- rowMeans(apply(AllPrediction, 2, function(x) exp(x)-lambda))
  } else if(transfo=="sqrt") {
    FinalForecast <- FinalForecast^2-lambda
    FinalPrediction <- rowMeans(apply(AllPrediction, 2, function(x) x^2-lambda))
  } else {
    
    FinalPrediction <- rowMeans(AllPrediction, na.rm = T) 
  }
  

  
  Accuracy_Train <- accuracy(FinalPrediction, train_o)
  Accuracy_Test <- accuracy(FinalForecast, test_o)

  return(list(Accuracy.Train = Accuracy_Train,
              Accuracy.Test = Accuracy_Test,
              
              FinalForecast=FinalForecast,
              FinalPrediction=FinalPrediction
              ))
}

