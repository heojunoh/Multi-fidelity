### nonlinear Example ###
library(lhs)
library(laGP)

install.packages("devtools")
library(devtools)
install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(lhs)
library(plgp)

library(ggplot2)
library(gridExtra)
library(dplyr)
library(patchwork) # To display 2 charts together
library(hrbrthemes)
extrafont::font_import() 
extrafont::loadfonts()

### synthetic function ###
f1 <- function(x)
{
  sin(8*pi*x)
}

f2 <- function(x)
{ 
  (x-sqrt(2))*(sin(8*pi*x))^2
}

n1 <- 13; n2 <- 8
eps <- sqrt(.Machine$double.eps)

set.seed(2)
X1 <- maximinLHS(n1, 1)
X2 <- maximinLHS(n2, 1)

NestDesign <- NestedDesignBuild(design = list(X1,X2))

X1 <- NestDesign$PX
X2 <- ExtractNestDesign(NestDesign,2)

y1 <- f1(X1)
y2 <- f2(X2)

### test data ###
x <- seq(0,1,length.out=100)

### RNAM ###
fit.RNAM <- RNAM(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
predy <- predRNAM(fit.RNAM, x)$mu
predsig2 <- predRNAM(fit.RNAM, x)$sig2

### cokm ###
fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                         lower=eps, upper=0.1,
                         coef.trend = list(0,c(0,0)), response = list(y1,y2), nlevel = 2)
pred.muficokm <- predict(fit.muficokm, x, "SK")






p13 <- ggplot(data.frame(x=rep(x,3), y=c(predy, f1(x), f2(x)), 
                         Linetype=factor(c(rep("RNA", 100), rep("f1", 100), rep("f2", 100)), 
                                         levels=c("RNA","f1","f2")), 
                         xx=c(rep(NA,100), X1, rep(NA, 100-nrow(X1)), X2, rep(NA, 100-nrow(X2))),
                         yy=c(rep(NA,100), y1, rep(NA, 100-length(y1)), y2, rep(NA, 100-nrow(y2)))), aes(x=x), color=group) +
  theme_ipsum() + 
  theme(axis.title.x = element_text(size=20,margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=20,margin = margin(t = 10), hjust=0.5),
        text=element_text(size=16,  family="serif")
  )+
  labs(x="Proposed", y = "y")+ 
  geom_point(aes(x=xx, y=yy, group=Linetype, shape=Linetype, color=Linetype), size=3)+
  geom_line(aes(y=y, group=Linetype, color=Linetype))+ 
  scale_color_manual(name = "Line type", values = c(f1 = "red", f2 = "green", RNA="blue")) +
  scale_shape_manual(name = "Line type", values = c(f1 = 17, f2 = 16, RNA=NA))+
  geom_ribbon(data=data.frame(x=x, predy=predy), aes(ymin=predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), ymax=predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2))), alpha=0.1, 
              color = "grey", linetype = "blank")+
  scale_y_continuous(
    limits = c(-2,1.3)
  )

p13


p14 <- ggplot(data.frame(x, pred.muficokm$mean), aes(x=x), color=group) +
  theme_ipsum() + 
  theme(axis.title.x = element_text(size=20,margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(margin = margin(t = 10), hjust=0.5),
        text=element_text(size=16,  family="serif")
  )+
  labs(x="CoKriging", y = "")+ 
  geom_line(aes(y=f2(x)), color="red")+ 
  geom_line(aes(y=pred.muficokm$mean), color="blue") + 
  geom_point(data=data.frame(X1,y1),aes(x=X1, y=y1),col="red", shape=2, size=1.5) +
  geom_point(data=data.frame(X2,y2),aes(x=X2, y=y2),col="green", shape=19, size=1.5) +
  geom_ribbon(aes(ymin=pred.muficokm$mean-1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2)), ymax=pred.muficokm$mean+1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2))), alpha=0.1, 
              color = "grey", linetype = "blank")+
  scale_y_continuous(
    limits = c(-2,1.3),
    # # Features of the first axis
    # name = "predictive posterior mean",
    # 
    # # Add a second axis and specify its features
    # sec.axis = sec_axis(~./coeff)#, name="Uncertainty reduction")
  )

p14




p18 <- grid.arrange(p13, p14, nrow=1, widths=c(1,1))






