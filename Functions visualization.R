#######################################
### Changed functions visualization ###
#######################################
par(mfrow=c(2,2))

# 3-layers
fl <- function(x, l){
  term1 <- sin(2*pi*x)
  term2 <- 0.2 * sin(8*pi*x)
  
  term1 + term2*5*0.8^l + (term1+term2)^3 + exp(-2*term1*term2) 
}
x <- seq(0, 1, length = 100)
y <- fl(x, 1)
plot(x, y, main="3-layers", type="l", col="red")
curve(fl(x,l=3),add=TRUE, col="orange",lty=1) # high fidelity(TRUE); Black
curve(fl(x,l=5),add=TRUE, col="green",lty=1) # high fidelity(TRUE); Black
curve(fl(x,l=Inf),add=TRUE, col="black",lty=1) # high fidelity(TRUE); Black

# Currin
Currin2 <- function(x,y){
  term1 <- (1-exp(-1/(2*y)))*(2300*x^3+1900*x^2+2092*x+60)/(100*x^3+500*x^2+4*x+20)
  term2 <- exp(-1.4*x)*cos(3.5*pi*y)
  
  # term1 + term2 * 16*2^(-l) 
  (term1 * (1+0.8^1))
}
x <- y <- seq(0, 1, length = 30)
z <- outer(x, y, Currin2)
persp(x, y, z, shade=.4, theta = 50, phi = 15, ticktype='detailed', main="Currin2")

# Gramacy1
Gramacy1 <- function(x,y){
# 10*x1*exp(-x1^2-x2^2)
  10*x*exp(-x^2-y^2)* (1+0.8^1) 
}
x <- y <- seq(-2, 4, length = 30)
z <- outer(x, y, Gramacy1)
persp(x, y, z, shade=.4, theta = 30, phi = 15, ticktype='detailed', main="Gramacy1")

# Gramacy2
Gramacy2 <- function(x,y,l=Inf){
  term1 <- exp(-(x-1)^2) + exp(-0.8*(x+1)^2) - 0.05*sin(8*(x+0.1))
  term2 <- exp(-(y-1)^2) + exp(-0.8*(y+1)^2) - 0.05*sin(8*(y+0.1))
  
  # - term1 * term2
  - term1 * term2 * (1+0.8^1)
}
x <- y <- seq(-2, 2, length = 30)
z <- outer(x, y, Gramacy2)
persp(x, y, z, shade=.4, theta = 30, phi = 15, ticktype='detailed', main="Gramacy2")


########################################
### Original functions visualization ###
########################################
par(mfrow=c(1,3))

# Currin
Currin2 <- function(x,y){
  term1 <- (1-exp(-1/(2*y)))*(2300*x^3+1900*x^2+2092*x+60)/(100*x^3+500*x^2+4*x+20)
  term2 <- exp(-1.4*x)*cos(3.5*pi*y)
  
  # term1 + term2 * 16*2^(-l) 
  term1 + 16*term2 * 0.5
}
x <- y <- seq(0, 1, length = 30)
z <- outer(x, y, Currin2)
persp(x, y, z, shade=.4, theta = 50, phi = 15, ticktype='detailed', main="Currin2")

# Gramacy1
Gramacy1 <- function(x,y){
  # 10*x1*exp(-x1^2-x2^2)
  10*x*exp(-x^2-y^2) 
}
x <- y <- seq(-2, 4, length = 30)
z <- outer(x, y, Gramacy1)
persp(x, y, z, shade=.4, theta = 30, phi = 15, ticktype='detailed', main="Gramacy1")

# Gramacy2
Gramacy2 <- function(x,y,l=3){
  term1 <- exp(-(x-1)^2) + exp(-0.8*(x+1)^2) - 0.05*sin(8*(x+0.1))
  term2 <- exp(-(y-1)^2) + exp(-0.8*(y+1)^2) - 0.05*sin(8*(y+0.1))
  
  # - term1 * term2
  - term1 * term2
}
x <- y <- seq(-2, 2, length = 30)
z <- outer(x, y, Gramacy2)
persp(x, y, z, shade=.4, theta = 30, phi = 15, ticktype='detailed', main="Gramacy2")




# Branin
branin <- function(x, y){
  (-1.275*x^2/pi^2+5*x/pi+y-6)^2 + (10-5/(4*pi))*cos(x)+ 10
  }

braninm <- function(x, y)
{ 
  10*sqrt((-1.275*(x+2)^2/pi^2+5*(x+2)/pi+(y+2)-6)^2 + (10-5/(4*pi))*cos((x+2))+ 10) + 2*(x-0.5) - 3*(3*y-1) - 1
  }

braninl <- function(x, y)
{ 10*sqrt((-1.275*(1.2*x+0.4)^2/pi^2+5*(1.2*x+0.4)/pi+(1.2*y+0.4)-6)^2 + (10-5/(4*pi))*cos((1.2*x+0.4))+ 10) + 2*(1.2*x+1.9) - 3*(3*(1.2*y+2.4)-1) - 1 - 3*y + 1
}

x <- seq(-5, 10, length = 30)
y <- seq(0, 15, length = 30)
z <- outer(x, y, branin)
zm <- outer(x, y, braninm)
zl <- outer(x, y, braninl)

persp(x, y, zl, shade=.01, theta = 330, phi = 25, ticktype='detailed', col="red")
par(new=TRUE)
persp(x, y, zm, shade=.01, theta = 330, phi = 25, ticktype='detailed', col="orange")
par(new=TRUE)
persp(x, y, z, shade=.01, theta = 330, phi = 25, ticktype='detailed', col="green", main="Branin")



library(lattice)
library(akima)


SurfaceData <- data.frame(
  x=rep(x,each=30,times=3),
  y=rep(y,30*3),
  z=c(t(z), t(zm), t(zl)),
  type=factor(rep(c("High","Medium","Low"),each=900))
)

SurfaceData$type <- factor(SurfaceData$type, levels = c("High", "Medium", "Low"))

wireframe(z~x*y,data=SurfaceData,group=type,
          col.groups=rgb(c(255,0,0), 
                         c(0,255,0), 
                         c(0,0,255),alpha = 120,maxColorValue = 255) ,
          scales = list(arrows=FALSE, col="black",font=10),
          xlab = list("x1",rot=30),
          ylab = list("x2",rot=-30),
          zlab = list("f",rot=90),
          zlim = c(-195,310),
          key=list(space="bottom",
            text=list(c("High","Medium","Low"),col=c(1,1,1)),
                   lines=list(lty=c(1,1,1),lwd=c(6,6,6),col=rgb(c(255,0,0), 
                                                   c(0,255,0), 
                                                   c(0,0,255),maxColorValue = 255))),
          par.settings = list(axis.line = list(col = "transparent")),
          par.box = list(col=c(1,0,0)),col=1
          )


### computation time of fitting and prediction ### 
### proposed; (0.132, 0.325, 0.595) ### 
### Cokm; (0.033, 0.057, 0.059) ### 
### NARGP; (4.357, 38.488(reconstraining ), 21.158) ### 



