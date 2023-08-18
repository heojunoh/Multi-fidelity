IMSPEselect2 <- function(newx, fit, level){
  ### Updating Ki using true output (not MC sample) ###
  # newx; new point which will be added, should be n=1
  # fit1; first layer's emulator
  # fit2; second layer's emulator f_M(X2, f1(x2))
  
  fit1 <- fit$fit1
  fit2 <- fit$fit2
  constant <- fit$constant
  kernel <- fit$kernel
  
  newx <- matrix(newx, nrow=1)
  
  X1 <- t(t(fit1$X)*attr(fit1$X,"scaled:scale")+attr(fit1$X,"scaled:center"))
  X2 <- matrix(t(t(fit2$X)*attr(fit2$X,"scaled:scale")+attr(fit2$X,"scaled:center"))[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)
  
  if(constant){
    y1 <- fit1$y
    y2 <- fit2$y
  }else{
    y1 <- fit1$y+attr(fit1$y,"scaled:center")
    y2 <- fit2$y+attr(fit2$y,"scaled:center")
  }
  
  
  ### Choose level 1 ###
  if(level == 1){
    d1 <- data.frame(newx*0.5+0.25, rep(0.05, 1)) # scale X to [-1,1]
    write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                      splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                      intern = TRUE)
    d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    y1.select <- d2$V4
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
  }
  
  ### Choose level 2 ###
  if(level == 2){
    d1 <- data.frame(newx*0.5+0.25, rep(0.05, 1)) # scale X to [-1,1]
    write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                      splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                      intern = TRUE)
    d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    y1.select <- d2$V4
    
    d1 <- data.frame(newx*0.5+0.25, rep(0.025, 1)) # scale X to [-1,1]
    write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                      splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                      intern = TRUE)
    d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    y2.select <- d2$V4
    
    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
    X2 <- rbind(X2, newx)
    y2 <- c(y2, y2.select)
  }
  
  fit <- RNAmf(X1, y1, X2, y2, kernel=kernel, constant=constant)
  
  return(list(fit=fit))
}

