### nonlinear Example ###
library(lhs)
library(laGP)
source("GP.R")
source("KOH.R")
source("closed.R")
source("score.R")
source("IMSPE1.R")

install.packages("devtools")
library(devtools)
install_github("cran/MuFiCokriging")
library(MuFiCokriging)
library(lhs)
library(plgp)

library(ggplot2)
library(dplyr)
install.packages("patchwork")
library(patchwork) # To display 2 charts together
install.packages("hrbrthemes")
library(hrbrthemes)
library(gridExtra)
library(grid)

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

### closed ###
fit.RNAmf <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
predy <- predRNAmf(fit.RNAmf, x)$mu
predsig2 <- predRNAmf(fit.RNAmf, x)$sig2

### cokm ###
fit.muficokm <- MuFicokm(formula = list(~1,~1), MuFidesign = NestDesign, covtype="gauss",
                         lower=eps, upper=0.1,
                         coef.trend = list(0,c(0,0)), response = list(y1,y2), nlevel = 2)
pred.muficokm <- predict(fit.muficokm, x, "SK")


NARmean <- c(-0.64842572,-0.64401556,-0.63456949,-0.63451209,-0.59112004,
             -0.61897197,-0.60394271,-1.23960041,-0.58514688,-0.58001279,
             -0.58652095,-0.57287946,-0.56885379,-0.57381196,-0.55603445,
             -0.5555246 ,-0.54523126,-0.52288103,-0.53652772,-0.51696515,
             -0.5013992 ,-0.50263074,-0.50920289,-0.50144504,-0.48966226,
             -0.47428948,-0.47099004,-0.44352044,-0.79849489,-0.44764871,
             -0.4407007 ,-0.42793549,-0.42667296,-0.44340277,-0.40074396,
             -0.41572591,-0.400692  ,-0.28004983,-0.38577961,-0.41225845,
             -0.3884197 ,-0.37814481,-0.35728287,-0.37283161,-0.35919018,
             -0.33645263,-0.34439913,-0.33378732,-0.3283503 ,-0.3160327 ,
             -0.30883575,-0.30158666,-0.29479207,-0.29638845,-0.28041887,
             -0.28888706,-0.55385623,-0.25586533,-0.2621796 ,-0.27564454,
             -0.24045041,-0.22661068,-0.23557057,-0.21720306,-0.2147197 ,
             -0.23602073,-0.19965745,-0.2041364 ,-0.19005733,-0.18146192,
             -0.16178178,-0.16284502,-0.14741647,-0.0273453 ,-0.16159134,
             -0.12272423,-0.13568855,-0.12988432,-0.39566698,-0.11501883,
             -0.10761918,-0.10641518,-0.11008389,-0.10335874,-0.06756379,
             -0.06229448,-0.06814745,-0.00774291,-0.06113242,-0.06255629,
             -0.02063464,-0.04125828,-0.03314737,-0.03288397,-0.0302045 ,
             -0.00871872,-0.21223052,-0.00659858, 0.02236903, 0.00681209)

NARsig2 <- c(0.18493265,0.19002683,0.19631839,0.1874243 ,0.18459319,
             0.19034009,0.18231581,0.02686901,0.1789189 ,0.1785596 ,
             0.18516496,0.18567901,0.17486536,0.18496372,0.17614526,
             0.17547028,0.16895601,0.16826385,0.17556962,0.1745893 ,
             0.16702749,0.17504175,0.16629112,0.17274002,0.17103381,
             0.16074665,0.16601098,0.15653752,0.04205821,0.15973028,
             0.16478193,0.14847821,0.15888696,0.16506179,0.15407477,
             0.15591323,0.1603599 ,0.10447765,0.15924793,0.15494279,
             0.14958567,0.15330886,0.15364302,0.15793245,0.15757494,
             0.15543123,0.15505941,0.14954022,0.14649163,0.14862542,
             0.15069285,0.14043152,0.1433886 ,0.14542771,0.138974  ,
             0.14226575,0.11011231,0.1410255 ,0.14275524,0.14683637,
             0.1395714 ,0.13072613,0.13601204,0.14086881,0.13766053,
             0.15750309,0.13341537,0.13342081,0.13829342,0.12701006,
             0.12892809,0.13094207,0.1287667 ,0.08239719,0.13573912,
             0.12999785,0.12984666,0.12250858,0.01422206,0.12863803,
             0.12865252,0.13099914,0.13003553,0.12894386,0.12231389,
             0.12223731,0.12629793,0.10373919,0.11711038,0.11658316,
             0.12025417,0.129037  ,0.12541308,0.12842163,0.12549758,
             0.11931937,0.01631042,0.12460786,0.11780417,0.12039394)


Icurrent <- mean(predsig2) # current IMSPE
Icurrent

### ALC ###
integvarlist <- integvar(x, fit.RNAmf, mc.sample=5)
which.min(integvarlist$intvar1)
which.min(integvarlist$intvar2)
c(Icurrent - integvarlist$intvar1[which.min(integvarlist$intvar1)], Icurrent - integvarlist$intvar2[which.min(integvarlist$intvar2)])
which.max(c(Icurrent - integvarlist$intvar1[which.min(integvarlist$intvar1)], Icurrent - integvarlist$intvar2[which.min(integvarlist$intvar2)])/c(1,(1+6)))

Iselect1 <- IMSPEselect1(x[which.min(integvarlist$intvar1)], fit.RNAmf, level=1)
predy1 <- predRNAmf(Iselect1$fit, x)$mu
predsig21 <- predRNAmf(Iselect1$fit, x)$sig2


### ALMC ###
integvarlist2 <- integvar(x[which.max(predsig2)], fit.RNAmf, mc.sample=5)
which.max(predsig2)
x[which.max(predsig2)]
c(Icurrent - integvarlist2$intvar1[which.min(integvarlist2$intvar1)], Icurrent - integvarlist2$intvar2[which.min(integvarlist2$intvar2)])
which.max(c(Icurrent - integvarlist2$intvar1[which.min(integvarlist2$intvar1)], Icurrent - integvarlist2$intvar2[which.min(integvarlist2$intvar2)])/c(1,(1+6)))

Iselect2 <- IMSPEselect1(x[which.max(predsig2)], fit.RNAmf, level=1)
predy2 <- predRNAmf(Iselect2$fit, x)$mu
predsig22 <- predRNAmf(Iselect2$fit, x)$sig2


### ALM ###
which.max(pred.GP(fit.RNAmf$fit1, x)$sig2)
which.max(predsig2)
x[which.max(pred.GP(fit.RNAmf$fit1, x)$sig2)]
x[which.max(predsig2)]
which.max(c(max(pred.GP(fit.RNAmf$fit1, x)$sig2), max(predsig2))/c(1,(1+6)))

Iselect3 <- IMSPEselect1(x[which.max(pred.GP(fit.RNAmf$fit1, x)$sig2)], fit.RNAmf, level=1)
predy3 <- predRNAmf(Iselect3$fit, x)$mu
predsig23 <- predRNAmf(Iselect3$fit, x)$sig2




plot(x, Icurrent-integvarlist$intvar1, type="l", lwd=2, col=3, ylim=range(Icurrent-integvarlist$intvar1))
plot(x, Icurrent-integvarlist$intvar2, type="l", lwd=2, col=3, ylim=range(Icurrent-integvarlist$intvar2))

which.min(integvarlist$intvar1)
which.min(integvarlist$intvar2)

### Fast update; Equation 6.6. in Surrogates ###
### ALC; How much can be improved. Equation 6.6. in Surrogates ###
alcfast <- c(Icurrent - integvarlist$intvar1[which.min(integvarlist$intvar1)], Icurrent - integvarlist$intvar2[which.min(integvarlist$intvar2)])
alcfast

### cost; 1, 2, 3 ###
which.max(alcfast/c(2,(2+8)))
alcfast/c(2,(2+8))


### Plotting the chosen point ###
# plot(x, predy, type="l", lwd=2, col=3, ylim=c(-2,1)) 
# lines(x, predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
# lines(x, predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), col=3, lty=2)
# 
# curve(f2(x),add=TRUE, col=1,lwd=2,lty=2) # high fidelity(TRUE); Black
# 
# points(X1, y1, pch="1", col="red")
# points(X2, y2, pch="2", col="red")
# 
# text(g[which.min(Icand1fast)], predy[which.min(Icand1fast)], expression("1*"), col="red")
# text(g[which.min(Icand2fast)], predy[which.min(Icand2fast)], expression("2*"), col="red")
# 
# 
# lines(g, Icand1fast, lwd=2, col=3, ylim=range(Icand1fast))
# lines(g, Icand2fast, lwd=2, col=3, ylim=range(Icand2fast))


### Plotting the model, design point, and IMSPE ###


# p13 <- ggplot(data.frame(x=rep(x,3), y=c(predy, f1(x), f2(x)), 
#                          Linetype=factor(c(rep("RNA", 100), rep("f1", 100), rep("f2", 100)), 
#                                          levels=c("RNA","f1","f2")), 
#                          xx=c(rep(NA,100), X1, rep(NA, 100-nrow(X1)), X2, rep(NA, 100-nrow(X2))),
#                          yy=c(rep(NA,100), y1, rep(NA, 100-length(y1)), y2, rep(NA, 100-nrow(y2)))), aes(x=x), color=group) +
#   theme_ipsum() + 
#   theme(axis.title.x = element_text(size=20,margin = margin(t = 10), hjust=0.5),
#         axis.title.y = element_text(size=20,margin = margin(t = 10), hjust=0.5),
#         text=element_text(size=16,  family="serif")
#   )+
#   labs(x="Proposed", y = "y")+ 
#   geom_point(aes(x=xx, y=yy, group=Linetype, shape=Linetype, color=Linetype), show.legend=FALSE, size=2)+
#   geom_line(aes(y=y, group=Linetype, color=Linetype), show.legend=FALSE)+ 
#   scale_color_manual(name = "Line type", values = c(f1 = "red", f2 = "green", RNA="blue")) +
#   scale_shape_manual(name = "Line type", values = c(f1 = 2, f2 = 19, RNA=NA))+
#   geom_ribbon(data=data.frame(x=x, predy=predy), aes(ymin=predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), ymax=predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2))), alpha=0.1, 
#               color = "grey", linetype = "blank")+
#   scale_y_continuous(
#     limits = c(-2,1.3)
#   )
# 
# p13


# p14 <- ggplot(data.frame(x=rep(x,3), y=c(pred.muficokm$mean, f1(x), f2(x)), 
#                          Linetype=factor(c(rep("Cokriging", 100), rep("f1", 100), rep("f2", 100)), 
#                                          levels=c("Cokriging","f1","f2")), 
#                          xx=c(rep(NA,100), X1, rep(NA, 100-nrow(X1)), X2, rep(NA, 100-nrow(X2))),
#                          yy=c(rep(NA,100), y1, rep(NA, 100-length(y1)), y2, rep(NA, 100-nrow(y2)))), aes(x=x), color=group) +
#   theme_ipsum() + 
#   theme(axis.title.x = element_text(size=20,margin = margin(t = 10), hjust=0.5),
#         axis.title.y = element_text(size=20,margin = margin(t = 10), hjust=0.5),
#         text=element_text(size=16,  family="serif")
#   )+
#   labs(x="Cokriging", y = "y")+ 
#   geom_point(aes(x=xx, y=yy, group=Linetype, shape=Linetype, color=Linetype), show.legend=FALSE, size=2)+
#   geom_line(aes(y=y, group=Linetype, color=Linetype), show.legend=FALSE)+ 
#   scale_color_manual(name = "Line type", values = c(f1 = "red", f2 = "green", RNA="blue")) +
#   scale_shape_manual(name = "Line type", values = c(f1 = 2, f2 = 19, RNA=NA))+
#   geom_ribbon(data=data.frame(x=x, predy=pred.muficokm$mean), aes(ymin=predy-1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2)), ymax=predy+1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2))), alpha=0.1, 
#               color = "grey", linetype = "blank")+
#   scale_y_continuous(
#     limits = c(-2,1.3)
#   )
# 
# p14


# p15 <- ggplot(data.frame(x=rep(x,3), y=c(NARmean, f1(x), f2(x)), 
#                          Linetype=factor(c(rep("NARGP", 100), rep("f1", 100), rep("f2", 100)), 
#                                          levels=c("NARGP","f1","f2")), 
#                          xx=c(rep(NA,100), X1, rep(NA, 100-nrow(X1)), X2, rep(NA, 100-nrow(X2))),
#                          yy=c(rep(NA,100), y1, rep(NA, 100-length(y1)), y2, rep(NA, 100-nrow(y2)))), aes(x=x), color=group) +
#   theme_ipsum() + 
#   theme(axis.title.x = element_text(size=20,margin = margin(t = 10), hjust=0.5),
#         axis.title.y = element_text(size=20,margin = margin(t = 10), hjust=0.5),
#         text=element_text(size=16,  family="serif")
#   )+
#   labs(x="NARGP", y = "y")+ 
#   geom_point(aes(x=xx, y=yy, group=Linetype, shape=Linetype, color=Linetype), show.legend=FALSE, size=2)+
#   geom_line(aes(y=y, group=Linetype, color=Linetype), show.legend=FALSE)+ 
#   scale_color_manual(name = "Line type", values = c(f1 = "red", f2 = "green", RNA="blue")) +
#   scale_shape_manual(name = "Line type", values = c(f1 = 2, f2 = 19, RNA=NA))+
#   geom_ribbon(data=data.frame(x=x, predy=NARmean), aes(ymin=predy-1.96*sqrt(NARsig2*length(y2)/(length(y2)-2)), ymax=predy+1.96*sqrt(NARsig2*length(y2)/(length(y2)-2))), alpha=0.1, 
#               color = "grey", linetype = "blank")+
#   scale_y_continuous(
#     limits = c(-2,1.3)
#   )


p13 <- ggplot(data.frame(x, predy), aes(x=x), color=group) +
  theme_bw() + 
  theme(axis.title.x = element_text(size=15,margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=15,margin = margin(t = 10), hjust=0.5),
        panel.border = element_blank(),
        legend.position="none"
  )+
  labs(x="RNAmf", y = "y")+ 
  geom_line(aes(y=f1(x), color="low-fidelity"), linetype="dashed")+ 
  geom_line(aes(y=f2(x), color="high-fidelity"))+ 
  geom_line(aes(y=predy, color="prediction"), size=1.5) +
  scale_color_manual(name = "", values = c("low-fidelity" = scales::hue_pal()(3)[1], "high-fidelity" = "black", "prediction"=scales::hue_pal()(3)[3])) +
  scale_shape_manual(name = "", values = c("low-fidelity" = 1, "high-fidelity" = 2, "prediction"=NA)) +
  geom_point(data=data.frame(X1,y1),aes(x=X1, y=y1, color="low-fidelity", shape="low-fidelity"), size=3) +
  geom_point(data=data.frame(X2,y2),aes(x=X2, y=y2, color="high-fidelity", shape="high-fidelity"), size=3) +
  geom_point(data=data.frame(xx = 0, yy= 5),aes(x=xx, y=yy, color="prediction", shape="prediction"), size=3) +
  geom_ribbon(aes(ymin=predy-1.96*sqrt(predsig2*length(y2)/(length(y2)-2)), ymax=predy+1.96*sqrt(predsig2*length(y2)/(length(y2)-2))), alpha=0.1, 
              color = "grey", linetype = "blank")+
  scale_y_continuous(
    limits = c(-1.8,1.3)) 


p13


p14 <- ggplot(data.frame(x, pred.muficokm$mean), aes(x=x), color=group) +
  theme_bw() + 
  theme(axis.title.x = element_text(size=13,margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=13,margin = margin(t = 10), hjust=0.5),
        panel.border = element_blank(),
        legend.position="none"
  )+
  labs(x="CoKriging", y = "y")+ 
  geom_line(aes(y=f1(x), color="low-fidelity"), linetype="dashed")+ 
  geom_line(aes(y=f2(x), color="high-fidelity"))+ 
  geom_line(aes(y=pred.muficokm$mean, color="prediction"), size=1.5) +
  scale_color_manual(name = "", values = c("low-fidelity" = scales::hue_pal()(3)[1], "high-fidelity" = "black", "prediction"=scales::hue_pal()(3)[3])) +   
  scale_shape_manual(name = "", values = c("low-fidelity" = 1, "high-fidelity" = 2, "prediction"=NA)) +  
  geom_point(data=data.frame(X1,y1),aes(x=X1, y=y1, color="low-fidelity", shape="low-fidelity"), size=3) +
  geom_point(data=data.frame(X2,y2),aes(x=X2, y=y2,color="high-fidelity", shape="high-fidelity"), size=3) +
  geom_point(data=data.frame(xx = 0, yy= 5),aes(x=xx, y=yy, color="prediction", shape="prediction"), size=3) +
  geom_ribbon(aes(ymin=pred.muficokm$mean-1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2)), ymax=pred.muficokm$mean+1.96*sqrt(pred.muficokm$sig2*length(y2)/(length(y2)-2))), alpha=0.1, 
              color = "grey", linetype = "blank")+
  scale_y_continuous(
    limits = c(-1.8,1.3))  

p14


p15 <- ggplot(data.frame(x, NARmean), aes(x=x), color=group) +
  theme_bw() + 
  theme(axis.title.x = element_text(size=13,margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=13,margin = margin(t = 10), hjust=0.5),
        panel.border = element_blank(),
        legend.position="none"
  )+
  labs(x="NARGP", y = "y")+ 
  geom_line(aes(y=f1(x), color="low-fidelity"), linetype="dashed")+ 
  geom_line(aes(y=f2(x), color="high-fidelity"))+ 
  geom_line(aes(y=NARmean, color="prediction"), size=1.5) +
  scale_color_manual(name = "", values = c("low-fidelity" = scales::hue_pal()(3)[1], "high-fidelity" = "black", "prediction"=scales::hue_pal()(3)[3])) +   
  scale_shape_manual(name = "", values = c("low-fidelity" = 1, "high-fidelity" = 2, "prediction"=NA)) +
  geom_point(data=data.frame(X1,y1),aes(x=X1, y=y1, color="low-fidelity", shape="low-fidelity"), size=3) +
  geom_point(data=data.frame(X2,y2),aes(x=X2, y=y2,color="high-fidelity", shape="high-fidelity"), size=3) +
  geom_point(data=data.frame(xx = 0, yy= 5),aes(x=xx, y=yy, color="prediction", shape="prediction"), size=3) +
  geom_ribbon(aes(ymin=NARmean-1.96*sqrt(NARsig2*length(y2)/(length(y2)-2)), ymax=NARmean+1.96*sqrt(NARsig2*length(y2)/(length(y2)-2))), alpha=0.1, 
              color = "grey", linetype = "blank")+
  scale_y_continuous(
    limits = c(-1.8,1.3))

p15


# p16 <- ggplot(data.frame(x, UR = c(Icurrent-integvarlist$intvar1, Icurrent-integvarlist$intvar2), level = c(rep("Low",100), rep("High",100)))
#               , aes(x=x), color=group) +
#   theme_ipsum() + 
#   theme(axis.title.x = element_text(size=20, margin = margin(t = 10), hjust=0.5),
#         axis.title.y = element_text(size=20, margin = margin(t = 10, r=10), hjust=0.5),
#         text=element_text(size=16,  family="serif"))+
#   labs(x="x", y = "Uncertainty Reduction")+ 
#   geom_line(aes(y=UR, group=level, color=level)) + 
#   scale_color_manual(values = c("green", "red"))+
#   geom_point(data=data.frame(x, UR=Icurrent-integvarlist$intvar1),aes(x=x[which.max(UR)], y=max(UR)), col="red", shape=2, size=2) +
#   geom_point(data=data.frame(x, UR=Icurrent-integvarlist$intvar2),aes(x=x[which.max(UR)], y=max(UR)), col="green", shape=19, size=2) +
#   scale_y_continuous(
#     limits = c(0,0.016),
#   )
# 
# p16


pp17 <- ggplot(data.frame(x, predy1), aes(x=x), color=group) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        text=element_text(size=16,  family="serif")
  )+
  labs(x="ALM", y = "y")+ 
  geom_line(aes(y=f1(x), color="low-fidelity"), linetype="dashed")+ 
  geom_line(aes(y=f2(x), color="high-fidelity"))+ 
  geom_line(aes(y=10, color="new point"))+ 
  geom_line(aes(y=predy1, color="prediction"), size=1.5) +
  scale_color_manual(name = "", values = c("low-fidelity" = scales::hue_pal()(3)[1], "high-fidelity" = "black", "prediction"=scales::hue_pal()(3)[3], "new point"=NA), breaks =c("low-fidelity", "high-fidelity", "prediction", "new point")) +
  scale_shape_manual(name = "", values = c("low-fidelity" = 1, "high-fidelity" = 2, "prediction"=32, "new point"=15), breaks =c("low-fidelity", "high-fidelity", "prediction", "new point")) +
  geom_point(data=data.frame(X1,y1),aes(x=X1, y=y1, color="low-fidelity", shape="low-fidelity"), size=3) +
  geom_point(data=data.frame(X2,y2),aes(x=X2, y=y2, color="high-fidelity", shape="high-fidelity"), size=3) +
  geom_point(data=data.frame(xx = 0, yy= 5),aes(x=xx, y=yy, color="prediction", shape="prediction"), size=3) +
  geom_point(data=data.frame(xx=x[which.max(pred.GP(fit.RNAmf$fit1, x)$sig2)],yy=f1(x[which.max(pred.GP(fit.RNAmf$fit1, x)$sig2)])),aes(x=xx, y=yy, color="low-fidelity", shape="new point"), size=3) +
  # geom_point(data=data.frame(xx=x[which.max(predsig2)],yy=f2(x[which.max(predsig2)])),aes(x=xx, y=yy, color="high-fidelity", shape="new point"), size=3) +
  geom_ribbon(aes(ymin=predy1-1.96*sqrt(predsig21*length(y2)/(length(y2)-2)), ymax=predy1+1.96*sqrt(predsig21*length(y2)/(length(y2)-2))), alpha=0.1, 
              color = "grey", linetype = "blank")+
  scale_y_continuous(
    limits = c(-1.8,1.3)) 

pp17


pp18 <- ggplot(data.frame(x, predy2), aes(x=x), color=group) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        text=element_text(size=16,  family="serif")
  )+
  labs(x="ALC", y = "y")+ 
  geom_line(aes(y=f1(x), color="low-fidelity"), linetype="dashed")+ 
  geom_line(aes(y=f2(x), color="high-fidelity"))+ 
  geom_line(aes(y=10, color="new point"))+ 
  geom_line(aes(y=predy2, color="prediction"), size=1.5) +
  scale_color_manual(name = "", values = c("low-fidelity" = scales::hue_pal()(3)[1], "high-fidelity" = "black", "prediction"=scales::hue_pal()(3)[3], "new point"=NA), breaks =c("low-fidelity", "high-fidelity", "prediction", "new point")) +
  scale_shape_manual(name = "", values = c("low-fidelity" = 1, "high-fidelity" = 2, "prediction"=32, "new point"=15), breaks =c("low-fidelity", "high-fidelity", "prediction", "new point")) +
  geom_point(data=data.frame(X1,y1),aes(x=X1, y=y1, color="low-fidelity", shape="low-fidelity"), size=3) +
  geom_point(data=data.frame(X2,y2),aes(x=X2, y=y2, color="high-fidelity", shape="high-fidelity"), size=3) +
  geom_point(data=data.frame(xx = 0, yy= 5),aes(x=xx, y=yy, color="prediction", shape="prediction"), size=3) +
  geom_point(data=data.frame(xx=x[which.min(integvarlist$intvar1)],yy=f1(x[which.min(integvarlist$intvar1)])),aes(x=xx, y=yy, color="low-fidelity", shape="new point"), size=3) +
  # geom_point(data=data.frame(xx=x[which.min(integvarlist$intvar2)],yy=f2(x[which.min(integvarlist$intvar2)])),aes(x=xx, y=yy, color="high-fidelity", shape="new point"), size=3) +
  geom_ribbon(aes(ymin=predy2-1.96*sqrt(predsig22*length(y2)/(length(y2)-2)), ymax=predy2+1.96*sqrt(predsig22*length(y2)/(length(y2)-2))), alpha=0.1, 
              color = "grey", linetype = "blank")+
  scale_y_continuous(
    limits = c(-1.8,1.3)) 

pp18


pp19 <- ggplot(data.frame(x, predy3), aes(x=x), color=group) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        text=element_text(size=16,  family="serif")
  )+
  labs(x="ALMC", y = "y")+ 
  geom_line(aes(y=f1(x), color="low-fidelity"), linetype="dashed")+ 
  geom_line(aes(y=f2(x), color="high-fidelity"))+ 
  geom_line(aes(y=10, color="new point"))+ 
  geom_line(aes(y=predy3, color="prediction"), size=1.5) +
  scale_color_manual(name = "", values = c("low-fidelity" = scales::hue_pal()(3)[1], "high-fidelity" = "black", "prediction"=scales::hue_pal()(3)[3], "new point"=NA), breaks =c("low-fidelity", "high-fidelity", "prediction", "new point")) +
  scale_shape_manual(name = "", values = c("low-fidelity" = 1, "high-fidelity" = 2, "prediction"=32, "new point"=15), breaks =c("low-fidelity", "high-fidelity", "prediction", "new point")) +
  geom_point(data=data.frame(X1,y1),aes(x=X1, y=y1, color="low-fidelity", shape="low-fidelity"), size=3) +
  geom_point(data=data.frame(X2,y2),aes(x=X2, y=y2, color="high-fidelity", shape="high-fidelity"), size=3) +
  geom_point(data=data.frame(xx = 0, yy= 5),aes(x=xx, y=yy, color="prediction", shape="prediction"), size=3) +
  geom_point(data=data.frame(xx=x[which.max(predsig2)],yy=f1(x[which.max(predsig2)])),aes(x=xx, y=yy, color="low-fidelity", shape="new point"), size=3) +
  # geom_point(data=data.frame(xx=x[which.max(predsig2)],yy=f2(x[which.max(predsig2)])),aes(x=xx, y=yy, color="high-fidelity", shape="new point"), size=3) +
  geom_ribbon(aes(ymin=predy3-1.96*sqrt(predsig23*length(y2)/(length(y2)-2)), ymax=predy3+1.96*sqrt(predsig23*length(y2)/(length(y2)-2))), alpha=0.1, 
              color = "grey", linetype = "blank")+
  scale_y_continuous(
    limits = c(-1.8,1.3)) 

pp19




# p.comptime <- (p13 | p14 | p15 | p17 | p18 | p19 + plot_layout(guides = "collect", ncol = 3) & theme(legend.position = "right")) / p15
# grid.arrange(p13, p14, p15, p16, nrow=2, widths=c(1,1,1), layout_matrix=rbind(c(1,2,3),c(4,4,4)))

(p13 + p14 + p15 ) / (pp17 + pp18 + pp19 + plot_layout(guides = "collect") & theme(legend.position = "bottom"))





### ALC ###
integvarlist <- integvar(x, fit.RNAmf, mc.sample=30)
which.min(integvarlist$intvar1)
which.min(integvarlist$intvar2)
c(Icurrent - integvarlist$intvar1[which.min(integvarlist$intvar1)], Icurrent - integvarlist$intvar2[which.min(integvarlist$intvar2)])
which.max(c(Icurrent - integvarlist$intvar1[which.min(integvarlist$intvar1)], Icurrent - integvarlist$intvar2[which.min(integvarlist$intvar2)])/c(1,(1+4)))

Iselect1 <- IMSPEselect1(x[which.min(integvarlist$intvar1)], fit.RNAmf, level=1)
predy1 <- predRNAmf(Iselect1$fit, x)$mu
predsig21 <- predRNAmf(Iselect1$fit, x)$sig2


### ALMC ###
integvarlist2 <- integvar(x[which.max(predsig2)], fit.RNAmf, mc.sample=5)
which.max(predsig2)
x[which.max(predsig2)]
c(Icurrent - integvarlist2$intvar1[which.min(integvarlist2$intvar1)], Icurrent - integvarlist2$intvar2[which.min(integvarlist2$intvar2)])
which.max(c(Icurrent - integvarlist2$intvar1[which.min(integvarlist2$intvar1)], Icurrent - integvarlist2$intvar2[which.min(integvarlist2$intvar2)])/c(1,(1+4)))

Iselect2 <- IMSPEselect1(x[which.max(predsig2)], fit.RNAmf, level=1)
predy2 <- predRNAmf(Iselect2$fit, x)$mu
predsig22 <- predRNAmf(Iselect2$fit, x)$sig2


### ALM ###
which.max(pred.GP(fit.RNAmf$fit1, x)$sig2)
which.max(predsig2)
x[which.max(pred.GP(fit.RNAmf$fit1, x)$sig2)]
x[which.max(predsig2)]
which.max(c(max(pred.GP(fit.RNAmf$fit1, x)$sig2), max(predsig2))/c(1,(1+4)))

Iselect3 <- IMSPEselect1(x[which.max(pred.GP(fit.RNAmf$fit1, x)$sig2)], fit.RNAmf, level=1)
predy3 <- predRNAmf(Iselect3$fit, x)$mu
predsig23 <- predRNAmf(Iselect3$fit, x)$sig2



p16 <- ggplot(data.frame(x, postvar = c(pred.GP(fit.RNAmf$fit1, x)$sig2, predsig2/(1+4), rep(NA, 100)), level = c(rep("low-fidelity",100), rep("high-fidelity",100), rep("predictive variance",100)))
              , aes(x=x), color=group) +
  theme_bw() +
  theme(axis.title.x = element_text(size=20, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        legend.position="none"
  )+
  labs(x="ALM", y = "ALM")+
  geom_line(aes(y=postvar, group=level, color=level)) +
  scale_color_manual(values = c("low-fidelity"= scales::hue_pal()(3)[1], "high-fidelity"= "black", "predictive variance"="grey"))+
  geom_point(data=data.frame(x, postvar=pred.GP(fit.RNAmf$fit1, x)$sig2),aes(x=x[which.max(postvar)], y=max(postvar)), col=scales::hue_pal()(3)[1], shape=2, size=2) +
  geom_point(data=data.frame(x, postvar=predsig2/(1+4)),aes(x=x[which.max(postvar)], y=max(postvar)), col="black", shape=19, size=2) +
  scale_y_continuous(
    limits = c(0,0.3),
  )

p16


p17 <- ggplot(data.frame(x, UR = c(Icurrent-integvarlist$intvar1, (Icurrent-integvarlist$intvar2)/(1+4), rep(NA, 100)), level = c(rep("low-fidelity",100), rep("high-fidelity",100), rep("predictive variance",100)))
              , aes(x=x), color=group) +
  theme_bw() +
  theme(axis.title.x = element_text(size=20, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        legend.position="none"
        )+
  labs(x="ALC", y = "ALC")+
  geom_line(aes(y=UR, group=level, color=level)) +
  scale_color_manual(values = c("low-fidelity"= scales::hue_pal()(3)[1], "high-fidelity"= "black", "predictive variance"="grey"))+
  geom_point(data=data.frame(x, UR=Icurrent-integvarlist$intvar1),aes(x=x[which.max(UR)], y=max(UR)), col=scales::hue_pal()(3)[1], shape=2, size=2) +
  geom_point(data=data.frame(x, UR=(Icurrent-integvarlist$intvar2)/(1+4)),aes(x=x[which.max(UR)], y=max(UR)), col="black", shape=19, size=2) +
  geom_line(data=data.frame(xx=seq(0,1,0.01),yy=rep(10,101)),aes(x=xx, y=yy), col="predictive variance") +
  scale_y_continuous(
    limits = c(0,0.016),
  )

p17


p18 <- ggplot(data.frame(x, UR = c(Icurrent-integvarlist$intvar1, (Icurrent-integvarlist$intvar2)/(1+4), rep(NA, 100)), level = c(rep("low-fidelity",100), rep("high-fidelity",100), rep("predictive variance",100)))
              , aes(x=x), color=group) +
  theme_bw() +
  theme(axis.title.x = element_text(size=20, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        legend.position="none",
        ) +
  labs(x="ALMC", y = "ALMC")+
  geom_line(aes(y=UR, group=level, color=level)) +
  geom_line(data=data.frame(x, postvar=predsig2), mapping= aes(y=postvar/20, color="predictive variance")) +
  geom_vline(xintercept = x[which.max(predsig2)], color="grey", linetype="dashed") +
  scale_color_manual(values = c("low-fidelity"= scales::hue_pal()(3)[1], "high-fidelity"= "black", "predictive variance"="grey"))+
  geom_point(data=data.frame(x, postvar=predsig2),aes(x=x[which.max(postvar)], y=(Icurrent-integvarlist$intvar1)[which.max(postvar)]), col=scales::hue_pal()(3)[1], shape=2, size=2) +
  geom_point(data=data.frame(x, postvar=predsig2),aes(x=x[which.max(postvar)], y=(Icurrent-integvarlist$intvar2)[which.max(postvar)]/(1+4)), col="black", shape=19, size=2) +
  scale_y_continuous(
    name = expression("ALMC"),
    limits = c(0,0.016)
    # ylim.prim <- c(0, 0.016),
    # ylim.sec <- c(0, 0.3),
    # sec.axis = sec_axis(~ .  * 20, name = "predictive variance")
    )

p18



# p.comptime <- (p13 | p14 | p15 | p17 | p18 | p19 + plot_layout(guides = "collect", ncol = 3) & theme(legend.position = "right")) / p15
# grid.arrange(p13, p14, p15, p16, nrow=2, widths=c(1,1,1), layout_matrix=rbind(c(1,2,3),c(4,4,4)))

(pp17 + pp18 + pp19 + plot_layout(guides = "collect") & theme(legend.position = "bottom")) / (p16 + p17 + (p18+guides(colour = "none")) + plot_layout(guides = "collect") & theme(legend.position = "bottom"))





