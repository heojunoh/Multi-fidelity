##########################
### Le Gratiet package ###
##########################

install.packages("devtools")
library(devtools)
install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(lhs)
library(plgp)

eps <- sqrt(.Machine$double.eps)

### synthetic function ###
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*0.8^l*5 + (term1+term2)^3 + exp(-2*term1*term2) #* (term1+term2)^3 # term1 + error term + interaction term
}

### training data ###
n1 <- 11; n2 <- 9; n3 <- 7
set.seed(1)
X1 <- maximinLHS(n1, 1)
X2 <- maximinLHS(n2, 1)
X3 <- maximinLHS(n3, 1)

NestDesign <- NestedDesignBuild(design = list(X1,X2,X3)) # Just make them nested for my setting.

X1 <- NestDesign$PX
X2 <- ExtractNestDesign(NestDesign,2)
X3 <- ExtractNestDesign(NestDesign,3)

y1 <- fl(X1, l=1)
y2 <- fl(X2, l=3)
y3 <- fl(X3, l=5)


### test data ###
x <- seq(0,1,0.01)

fit.muficokm <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = NestDesign, covtype="gauss",
                         # coef.trend = list(0,c(0,0),c(0,0)), # This is for mean function.
                         lower=0.001, upper=0.1, # You should adjust this lower and upper bound depending on the dataset or function. Sometimes cov is singular.
                         response = list(y1,y2,y3), nlevel = 3)
pred.muficokm <- predict(fit.muficokm, x, "SK", cov.compute=TRUE) # I use SK but there is UK too. See dicekriging package for it.

sqrt(mean((pred.muficokm$mean-fl(x, l=Inf))^2)) # RMSE




### cok, ZD, response, nlevel, Dnest, nuggets ###
fit.muficokm$cok
fit.muficokm$ZD
### c; cov between newdata and initial data
### Variance of each GP delta_l, rho ###
### cov; theta, var; sigma^2_l, trend; beta; constant trend
sum.muficokm <- summary(fit.muficokm)

sum.muficokm$Cov.Val
sum.muficokm$Var.Val
sum.muficokm$Rho.Val
sum.muficokm$Trend.Val


pred.muficokm$mean
pred.muficokm$mux[[1]]




y1.select <- fl(36, l=1)
fit.muficokm$Dnest$PX <- rbind(fit.muficokm$Dnest$PX, 36)

fit.muficokmnew <- MuFicokm(formula = list(~1,~1,~1), MuFidesign = fit.muficokm$Dnest, covtype="gauss",
                    # coef.trend = list(0,c(0,0),c(0,0)),
                    lower=0.001, upper=0.1,
                    response = list(rbind(fit.muficokm$cok[[1]]@y, y1.select),fit.muficokm$cok[[2]]@y,fit.muficokm$cok[[3]]@y), nlevel = 3)

mean(predict(fit.muficokm, x, "SK", cov.compute=TRUE)$sig2)
IMSPEcokm(g, 36, fit.muficokm)
mean(predict(fit.muficokmnew, x, "SK")$varx[[1]])
mean(predict(fit.muficokmnew, x, "SK")$varx[[2]])
mean(predict(fit.muficokmnew, x, "SK")$varx[[3]])

