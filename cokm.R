##########################
### Le Gartiet package ###
##########################

install.packages("devtools")
library(devtools)
install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(lhs)
library(plgp)


IMSPEcokm <- function(x, newx, fit){
  sumcokm <- summary(fit)
  
  k2 <- k3 <- rep(0, length(x))
  
  # k_n_1
  r1 <- covar.sep(x, rbind(fit$cok[[1]]@X, newx), d=(sumcokm$Cov.Val[[1]])^2, g=eps)
  R1 <- covar.sep(rbind(fit$cok[[1]]@X, newx), d=(sumcokm$Cov.Val[[1]])^2, g=eps)
  k1 <- diag(sumcokm$Var.Val[[1]]*( 1 - r1 %*% solve(R1) %*% t(r1) ))
  
  for(i in 1:length(x)){
    # k_n_2
    # F and f are 1; constant trend
    # H, R, F
    f2 <- 1
    F2 <- matrix(f2, ncol=1, nrow=nrow(fit$cok[[2]]@X))
    H2 <- rbind(fit$cok[[2]]@AR.F %*% f2, predict(fit, newx, "SK")$mux[[1]]) # y1(D2)
    h2 <- predict(fit, newx, "SK")$mean * f2   # predictive mean*f2
    r2 <- covar.sep(x[i], rbind(fit$cok[[2]]@X, newx), d=(sumcokm$Cov.Val[[2]])^2, g=eps)
    R2 <- covar.sep(rbind(fit$cok[[2]]@X, newx), d=(sumcokm$Cov.Val[[2]])^2, g=eps)
    
    sig2rho1 <- sumcokm$Rho.Val[[1]]^2 + solve(t(H2) %*% solve(R2) %*% H2)
    sig2delta2 <- sumcokm$Var.Val[[2]]*( 1 - cbind(h2, r2) %*% solve( rbind(cbind(0, t(H2)), cbind(H2, R2)) ) %*% t(cbind(h2, r2)) )
    
    k2[i] <- k1[i] * sig2rho1 + sig2delta2
    
    # k_n_3
    # F and f are 1; constant trend
    # H, R, F
    f3 <- 1
    F3 <- matrix(f3, ncol=1, nrow=nrow(fit$cok[[3]]@X))
    H3 <- rbind(fit$cok[[3]]@AR.F %*% f3, predict(fit, newx, "SK")$mux[[2]]) # y2(D3)
    h3 <- predict(fit, newx, "SK")$mean * f3   # predictive mean*f3
    r3 <- covar.sep(x[i], rbind(fit$cok[[3]]@X, newx), d=(sumcokm$Cov.Val[[3]])^2, g=eps)
    R3 <- covar.sep(rbind(fit$cok[[3]]@X, newx), d=(sumcokm$Cov.Val[[3]])^2, g=eps)
    
    sig2rho2 <- sumcokm$Rho.Val[[2]]^2 + solve(t(H3) %*% solve(R3) %*% H3)
    sig2delta3 <- sumcokm$Var.Val[[3]]*( 1 - cbind(h3, r3) %*% solve( rbind(cbind(0, t(H3)), cbind(H3, R3)) ) %*% t(cbind(h3, r3)) )
    
    k3[i] <- k2[i] * sig2rho2 + sig2delta3
  }
  
  return(list(k1=mean(k1), k2=mean(k2), k3=mean(k3)))
}


# IMSPEcokm1select <- function(x, newx, fit){
#   
#   y1.select <- fl(newx, l=1)
#   fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
#   
#   fit.new <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
#                            # coef.trend = list(0,c(0,0),c(0,0)),
#                            lower=0.001, upper=0.1,
#                            response = list(rbind(fit$cok[[1]]@y, y1.select),fit$cok[[2]]@y,fit$cok[[3]]@y), nlevel = 3)
#   
#   return(list(fit=fit.new))
# }


### This is updating model everytime.
IMSPEcokmselect <- function(x, newx, fit, level){
  
  ### choose level 1 ###
  if(level == 1){
    y1.select <- fl(newx, l=1)
    fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
    
    fit.new <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
                        # coef.trend = list(0,c(0,0),c(0,0)),
                        lower=0.001, upper=0.1,
                        response = list(rbind(fit$cok[[1]]@y, y1.select),fit$cok[[2]]@y,fit$cok[[3]]@y), nlevel = 3)
  }
  
  ### choose level 2 ###
  if(level == 2){
    y1.select <- fl(newx, l=1)
    y2.select <- fl(newx, l=3)
    fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
    fit$Dnest$ind[[1]] <- c(fit$Dnest$ind[[1]], nrow(fit$Dnest$PX))
    fit$Dnest$n[1] <- fit$Dnest$n[1]+1
    
    fit.new <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
                        # coef.trend = list(0,c(0,0),c(0,0)),
                        lower=0.001, upper=0.1,
                        response = list(rbind(fit$cok[[1]]@y, y1.select),rbind(fit$cok[[2]]@y, y2.select),fit$cok[[3]]@y), nlevel = 3)
  }
  
  ### choose level 3 ###
  if(level == 3){
    y1.select <- fl(newx, l=1)
    y2.select <- fl(newx, l=3)
    y3.select <- fl(newx, l=5)
    fit$Dnest$PX <- rbind(fit$Dnest$PX, newx)
    fit$Dnest$ind[[1]] <- c(fit$Dnest$ind[[1]], nrow(fit$Dnest$PX))
    fit$Dnest$ind[[2]] <- c(fit$Dnest$ind[[2]], nrow(fit$Dnest$PX))
    fit$Dnest$n[1] <- fit$Dnest$n[1]+1
    fit$Dnest$n[2] <- fit$Dnest$n[2]+1
    
    fit.new <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit$Dnest, covtype="gauss",
                        # coef.trend = list(0,c(0,0),c(0,0)),
                        lower=0.001, upper=0.1,
                        response = list(rbind(fit$cok[[1]]@y, y1.select),rbind(fit$cok[[2]]@y, y2.select),rbind(fit$cok[[3]]@y, y3.select)), nlevel = 3)
  }
  return(list(fit=fit.new))
}





######
n1 <- 9; n2 <- 7; n3 <- 5
set.seed(3)
X1 <- maximinLHS(n1, 1)
X2 <- maximinLHS(n2, 1)
X3 <- maximinLHS(n3, 1)

NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))

X1 <- NestDesign$PX
X2 <- ExtractNestDesign(NestDesign,2)
X3 <- ExtractNestDesign(NestDesign,3)

y1 <- fl(X1, l=1)
y2 <- fl(X2, l=3)
y3 <- fl(X3, l=5)

fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                         # coef.trend = list(0,c(0,0),c(0,0)),
                         lower=0.001, upper=0.1,
                         response = list(y1,y2,y3), nlevel = 3)



y1.select <- fl(g[which.min(Icand1cokm)], l=1)
y2.select <- fl(g[which.min(Icand1cokm)], l=3)
y3.select <- fl(g[which.min(Icand1cokm)], l=5)
fit.muficokm$Dnest
fit.muficokm$Dnest$PX <- rbind(fit.muficokm$Dnest$PX, g[which.min(Icand1cokm)])
fit.muficokm$Dnest$ind[[1]] <- c(fit.muficokm$Dnest$ind[[1]], nrow(fit.muficokm$Dnest$PX))
fit.muficokm$Dnest$n[1] <- fit.muficokm$Dnest$n[1]+1
fit.muficokm$Dnest$ind[[2]] <- c(fit.muficokm$Dnest$ind[[2]], nrow(fit.muficokm$Dnest$PX))
fit.muficokm$Dnest$n[2] <- fit.muficokm$Dnest$n[2]+1

fit.muficokm$Dnest






MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit.muficokm$Dnest, covtype="gauss",lower=0.001, upper=0.1,
         response = list(rbind(fit.muficokm$cok[[1]]@y, y1.select),fit.muficokm$cok[[2]]@y,fit.muficokm$cok[[3]]@y), nlevel = 3)
MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit.muficokm$Dnest, covtype="gauss",lower=0.001, upper=0.1,
         response = list(rbind(fit.muficokm$cok[[1]]@y, y1.select),rbind(fit.muficokm$cok[[2]]@y, y2.select),fit.muficokm$cok[[3]]@y), nlevel = 3)
MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit.muficokm$Dnest, covtype="gauss",lower=0.001, upper=0.1,
         response = list(rbind(fit.muficokm$cok[[1]]@y, y1.select),rbind(fit.muficokm$cok[[2]]@y, y2.select),rbind(fit.muficokm$cok[[3]]@y, y3.select)), nlevel = 3)





fit.muficokm$Dnest
NestDesign



MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit.muficokm$Dnest, covtype="gauss",lower=0.001, upper=0.1,
         response = list(rbind(fit.muficokm$cok[[1]]@y, y1.select),rbind(fit.muficokm$cok[[2]]@y, y2.select),rbind(fit.muficokm$cok[[3]]@y, y3.select)), nlevel = 3)
MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",lower=0.001, upper=0.1,
         response = list(y1,y2,y3), nlevel = 3)





km(~1, 
   design =fit.muficokm$Dnest$PX, 
   response = rbind(fit.muficokm$cok[[1]]@y, y1.select),
   covtype="gauss",lower=0.001, upper=0.1)
















