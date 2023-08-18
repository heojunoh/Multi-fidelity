IMSPEselect3 <- function(newx, fit, level){
  ### Updating Ki using true output (not MC sample) ###
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  
  fit1 <- fit$fit.RNAmf1$fit1
  fit2 <- fit$fit.RNAmf1$fit2
  fit3 <- fit$fit3
  constant <- fit$constant
  kernel <- fit$kernel
  
  newx <- matrix(newx, nrow=1)
  
  X1 <- t(t(fit1$X)*attr(fit1$X,"scaled:scale")+attr(fit1$X,"scaled:center"))
  X2 <- matrix(t(t(fit2$X)*attr(fit2$X,"scaled:scale")+attr(fit2$X,"scaled:center"))[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)
  X3 <- matrix(t(t(fit3$X)*attr(fit3$X,"scaled:scale")+attr(fit3$X,"scaled:center"))[,-ncol(fit3$X)], ncol=ncol(fit3$X)-1)
  
  if(constant){
    y1 <- fit1$y
    y2 <- fit2$y
    y3 <- fit3$y
  }else{
    y1 <- fit1$y+attr(fit1$y,"scaled:center")
    y2 <- fit2$y+attr(fit2$y,"scaled:center")
    y3 <- fit3$y+attr(fit3$y,"scaled:center")
  }
  
  
  ### Choose level 1 ###
  if(level == 1){
    # y1.select <- f1(newx) ### need to change for different function
    y1.select <- fl(newx, l=1)
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
  }
  
  ### Choose level 2 ###
  if(level == 2){
    # y1.select <- f1(newx) ### need to change for different function
    # y2.select <- f2(newx)
    y1.select <- fl(newx, l=1)
    y2.select <- fl(newx, l=3)
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
    X2 <- rbind(X2, newx)
    y2 <- c(y2, y2.select)
  }
  
  ### Choose level 3 ###
  if(level == 3){
    # y1.select <- f1(newx) ### need to change for different function
    # y2.select <- f2(newx)
    # y3.select <- f3(newx)
    y1.select <- fl(newx, l=1)
    y2.select <- fl(newx, l=3)
    y3.select <- fl(newx, l=5)
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
    X2 <- rbind(X2, newx)
    y2 <- c(y2, y2.select)
    X3 <- rbind(X3, newx)
    y3 <- c(y3, y3.select)
  }
  
  fit <- RNAmf2(X1, y1, X2, y2, X3, y3, kernel=kernel, constant=constant)
  
  return(list(fit=fit))
}

