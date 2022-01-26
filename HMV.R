# Hull-White/Vasicek, Ornstein-Uhlenbeck process

## HMV : dXt = µ(θ − Xt)dt + σdWt
## OU : dXt = −µXtdt + σdWt with µ; σ > 0.

## Hull-White/Vasicek Models
## dX(t) = 4 * (2.5 - X(t)) * dt + 1 *dW(t), X0=10
set.seed(1234)
X <- HWV(N=1000,M=10,mu = 4, theta = 2.5,sigma = 1,x0=10)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")


## Ornstein-Uhlenbeck Process
## dX(t) = -4 * X(t) * dt + 1 *dW(t) , X0=2
set.seed(1234)
X <- OU(N=1000,M=10,mu = 4,sigma = 1,x0=10)
plot(X,plot.type="single")
lines(as.vector(time(X)),rowMeans(X),col="red")



# Monthly Interest Rates

data(Irates)
rates <- Irates[,"r1"]
rates <- window(rates, start=1964.471, end=1989.333)
## CKLS modele vs CIR modele
## CKLS : dX(t) = (theta1+theta2* X(t))* dt + theta3 * X(t)^theta4 * dW(t)
fx <- expression(theta[1]+theta[2]*x)
gx <- expression(theta[3]*x^theta[4])
fitmod <- fitsde(rates,drift=fx,diffusion=gx,pmle="euler",
                 start = list(theta1=1,theta2=1,
                              theta3=1,theta4=1),optim.method = "L-BFGS-B")
theta <- coef(fitmod)
N <- length(rates)
res <- snssde1d(drift=fx,diffusion=gx,M=1000,t0=time(rates)[1],
                T=time(rates)[N],
                Dt=deltat(rates),x0=rates[1],N=N)

plot(res,type="n",ylim=c(0,35))
lines(rates,col=2,lwd=2)
lines(time(res),apply(res$X,1,mean),col=3,lwd=2)
lines(time(res),apply(res$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(res),apply(res$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topleft",c("real data","mean path",
                   paste("bound of", 95," confidence")),inset = .01,
       col=2:4,lwd=2,cex=0.8)
