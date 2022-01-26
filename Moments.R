# MCM.sde --- >Parallel Monte-Carlo Methods for SDE’s

## Example 1 : (1 dim)
## dX(t) = 3*(1-X(t)) dt + 0.5 * dW(t), X(0)=5, t in [0,10]
## set the model 1d
f <- expression(3*(1-x));g <- expression(0.5)
mod1d <- snssde1d(drift=f,diffusion=g,x0=5,T=10,M=50)
## function of the statistic(s) of interest.
sde.fun1d <- function(data, i){
  d <- data[i, ]
  return(c(mean(d),Mode(d),var(d)))
}
mc.sde1d = MCM.sde(model=mod1d,statistic=sde.fun1d,R=100,exact=list(Me=1,Mo=1,Va=0.5^2/6),
                   names=c("Me(10)","Mo(10)","Va(10)"))
mc.sde1d
plot(mc.sde1d,index=1)
plot(mc.sde1d,index=2)
plot(mc.sde1d,index=3)


## Example 2 : with Parallel computing

mod1d <- snssde1d(drift=f,diffusion=g,x0=5,T=10,M=1000)
## On Windows or Unix
mc.sde1d = MCM.sde(model=mod1d,statistic=sde.fun1d,R=1000,exact=list(Me=1,Mo=1,Va=0.5^2/6),
                   names=c("Me(10)","Mo(10)","Va(10)"),parallel="snow",ncpus=parallel::detectCores())
mc.sde1d
## On Unix only
mc.sde1d = MCM.sde(model=mod1d,statistic=sde.fun1d,R=1000,exact=list(Me=1,Mo=1,Va=0.5^2/6),
                   names=c("Me(10)","Mo(10)","Va(10)"),parallel="multicore",ncpus=parallel::detectCores())
mc.sde1d


## Example 3: (2 dim)
## dX(t) = 1/mu*(theta-X(t)) dt + sqrt(sigma) * dW1(t),
## dY(t) = X(t) dt + 0 * dW2(t)
## Not run:
## Set the model 2d
mu=0.75;sigma=0.1;theta=2
x0=0;y0=0;init=c(x=0,y=0)
f <- expression(1/mu*(theta-x), x)
g <- expression(sqrt(sigma),0)
OUI <- snssde2d(drift=f,diffusion=g,M=1000,Dt=0.01,x0=init)
## function of the statistic(s) of interest.
sde.fun2d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),var(d$x),var(d$y),cov(d$x,d$y)))
}
## Monte-Carlo at time = 5
mc.sde2d_a = MCM.sde(model=OUI,statistic=sde.fun2d,R=100,time=5,
                     parallel="snow",ncpus=parallel::detectCores())
mc.sde2d_a
## Monte-Carlo at time = 10
mc.sde2d_b = MCM.sde(model=OUI,statistic=sde.fun2d,R=100,time=10,
                     parallel="snow",ncpus=parallel::detectCores())
mc.sde2d_b
## Compared with exact values at time 5 and 10
E_x <- function(t) theta+(x0-theta)*exp(-t/mu)
V_x <- function(t) 0.5*sigma*mu *(1-exp(-2*(t/mu)))
E_y <- function(t) y0+theta*t+(x0-theta)*mu*(1-exp(-t/mu))
V_y <- function(t) sigma*mu^3*((t/mu)-2*(1-exp(-t/mu))+0.5*(1-exp(-2*(t/mu))))
cov_xy <- function(t) 0.5*sigma*mu^2 *(1-2*exp(-t/mu)+exp(-2*(t/mu)))
## at time=5
mc.sde2d_a = MCM.sde(model=OUI,statistic=sde.fun2d,R=100,time=5,
                     exact=list(m1=E_x(5),m2=E_y(5),S1=V_x(5),S2=V_y(5),C12=cov_xy(5)),
                     parallel="snow",ncpus=parallel::detectCores())
mc.sde2d_a
plot(mc.sde2d_a,index=1)
plot(mc.sde2d_a,index=2)
## at time=10

mc.sde2d_b = MCM.sde(model=OUI,statistic=sde.fun2d,R=100,time=10,
                     exact=list(m1=E_x(10),m2=E_y(10),S1=V_x(10),S2=V_y(10),C12=cov_xy(10)),
                     parallel="snow",ncpus=parallel::detectCores())
mc.sde2d_b
plot(mc.sde2d_b,index=1)
plot(mc.sde2d_b,index=2)
## End(Not run)
## Example 4: (3 dim)
## dX(t) = sigma*(Y(t)-X(t)) dt + 0.1 * dW1(t)
## dY(t) = (rho*X(t)-Y(t)-X(t)*Z(t)) dt + 0.1 * dW2(t)
## dZ(t) = (X(t)*Y(t)-bet*Z(t)) dt + 0.1 * dW3(t)
## W1(t), W2(t) and W3(t) are three correlated Brownian motions with Sigma
## Not run:
## Set the model 3d
sigma=10;rho=28; bet=8/3
f <- expression(sigma*(y-x),rho*x-y-x*z,x*y-bet*z)
g <- expression(0.1,0.1,0.1)
# correlation matrix
Sigma <-matrix(c(1,0.3,0.5,0.3,1,0.2,0.5,0.2,1),nrow=3,ncol=3)
mod3d <- snssde3d(x0=rep(0,3),drift=f,diffusion=g,M=1000,Dt=0.01,corr=Sigma)
## function of the statistic(s) of interest.
sde.fun3d <- function(data, i){
  d <- data[i,]
  return(c(mean(d$x),mean(d$y),mean(d$z)))
}
## Monte-Carlo at time = 10
mc.sde3d = MCM.sde(mod3d,statistic=sde.fun3d,R=100,parallel="snow",ncpus=parallel::detectCores())
mc.sde3d
## End(Not run)



# MEM.sde Moment Equations Methods for SDE

library(deSolve)
## Example 1: 1-dim
## dX(t) = mu * X(t) * dt + sigma * X(t) * dW(t)
## Symbolic ODE's of mean and variance
f <- expression(mu*x)
g <- expression(sigma*x)
res1 <- MEM.sde(drift=f,diffusion=g,type="ito")
res2 <- MEM.sde(drift=f,diffusion=g,type="str")
res1
res2
## numerical approximation of mean and variance
para <- c(mu=2,sigma=0.5)
t <- seq(0,1,by=0.001)
init <- c(m=1,S=0)
res1 <- MEM.sde(drift=f,diffusion=g,solve=TRUE,init=init,parms=para,time=t)
res1
matplot.0D(res1$sol.ode,main="Mean and Variance of X(t), type Ito")
plot(res1$sol.ode,select=c("m","S"))
## approximation at time = 0.75
summary(res1,at=0.75)
##
res2 <- MEM.sde(drift=f,diffusion=g,solve=TRUE,init=init,parms=para,time=t,type="str")
res2
matplot.0D(res2$sol.ode,main="Mean and Variance of X(t), type Stratonovich")
plot(res2$sol.ode,select=c("m","S"))
## approximation at time = 0.75
summary(res2,at=0.75)
## Comparison:
plot(res1$sol.ode, res2$sol.ode,ylab = c("m(t)"),select="m",xlab = "Time",
     col = c("red", "blue"))
plot(res1$sol.ode, res2$sol.ode,ylab = c("S(t)"),select="S",xlab = "Time",
     col = c("red", "blue"))


## Example2: 2-dim
## dX(t) = 1/mu*(theta-X(t)) dt + sqrt(sigma) * dW1(t),
## dY(t) = X(t) dt + 0 * dW2(t)
## Not run:
para=c(mu=0.75,sigma=0.1,theta=2)
init=c(m1=0,m2=0,S1=0,S2=0,C12=0)
t <- seq(0,10,by=0.001)
f <- expression(1/mu*(theta-x), x)
g <- expression(sqrt(sigma),0)
res2d <- MEM.sde(drift=f,diffusion=g,solve=TRUE,init=init,parms=para,time=t)
res2d
## Exact moment
mu=0.75;sigma=0.1;theta=2;x0=0;y0=0
E_x <- function(t) theta+(x0-theta)*exp(-t/mu)
V_x <- function(t) 0.5*sigma*mu *(1-exp(-2*(t/mu)))
E_y <- function(t) y0+theta*t+(x0-theta)*mu*(1-exp(-t/mu))
V_y <- function(t) sigma*mu^3*((t/mu)-2*(1-exp(-t/mu))+0.5*(1-exp(-2*(t/mu))))
cov_xy <- function(t) 0.5*sigma*mu^2 *(1-2*exp(-t/mu)+exp(-2*(t/mu)))
##
summary(res2d,at=5)
E_x(5);E_y(5);V_x(5);V_y(5);cov_xy(5)
matplot.0D(res2d$sol.ode,select=c("m1"))
curve(E_x,add=TRUE,col="red")
## plot
plot(res2d$sol.ode)
matplot.0D(res2d$sol.ode,select=c("S1","S2","C12"))
plot(res2d$sol.ode[,"m1"], res2d$sol.ode[,"m2"], xlab = "m1(t)",
     ylab = "m2(t)", type = "l",lwd = 2)
hist(res2d$sol.ode,select=c("m1","m2"), col = c("darkblue", "red", "orange", "black"))


## Example3: 2-dim with correlation
## Heston model
## dX(t) = mu*X(t) dt + sqrt(Y(t))*X(t) * dW1(t),
## dY(t) = lambda*(theta-Y(t)) dt + sigma*sqrt(Y(t)) * dW2(t)
## with E(dw1dw2)=rho
f <- expression( mu*x, lambda*(theta-y) )
g <- expression( sqrt(y)*x, sigma*sqrt(y) )
RHO <- expression(rho)
res2d <- MEM.sde(drift=f,diffusion=g,corr=RHO)
res2d
## Numerical approximation
RHO <- expression(0.5)
para=c(mu=1,lambda=3,theta=0.5,sigma=0.1)
ini=c(m1=10,m2=2,S1=0,S2=0,C12=0)
res2d = MEM.sde(drift=f,diffusion=g,solve=TRUE,parms=para,init=ini,time=seq(0,1,by=0.01))
res2d
matplot.0D(res2d$sol.ode,select=c("m1","m2"))
matplot.0D(res2d$sol.ode,select=c("S1","S2","C12"))



## Example4: 3-dim
## dX(t) = sigma*(Y(t)-X(t)) dt + 0.1 * dW1(t)
## dY(t) = (rho*X(t)-Y(t)-X(t)*Z(t)) dt + 0.1 * dW2(t)
## dZ(t) = (X(t)*Y(t)-bet*Z(t)) dt + 0.1 * dW3(t)
## with E(dw1dw2)=rho1, E(dw1dw3)=rho2 and E(dw2dw3)=rho3
f <- expression(sigma*(y-x),rho*x-y-x*z,x*y-bet*z)
g <- expression(0.1,0.1,0.1)
RHO <- expression(rho1,rho2,rho3)
## Symbolic moments equations
res3d = MEM.sde(drift=f,diffusion=g,corr=RHO)
res3d
## Numerical approximation
RHO <- expression(0.5,0.2,-0.7)
para=c(sigma=10,rho=28,bet=8/3)
ini=c(m1=1,m2=1,m3=1,S1=0,S2=0,S3=0,C12=0,C13=0,C23=0)
res3d = MEM.sde(drift=f,diffusion=g,solve=T,parms=para,init=ini,time=seq(0,1,by=0.01))
res3d
summary(res3d,at=0.25)
summary(res3d,at=0.50)
summary(res3d,at=0.75)
plot(res3d$sol.ode)
matplot.0D(res3d$sol.ode,select=c("m1","m2","m3"))
matplot.0D(res3d$sol.ode,select=c("S1","S2","S3"))
matplot.0D(res3d$sol.ode,select=c("C12","C13","C23"))
plot3D(res3d$sol.ode[,"m1"], res3d$sol.ode[,"m2"],res3d$sol.ode[,"m3"], xlab = "m1(t)",
       ylab = "m2(t)",zlab="m3(t)", type = "l",lwd = 2,box=F)


# moment ----> Monte-Carlo statistics of SDE’s

## Example 1:
## dX(t) = 2*(3-X(t)) *dt + dW(t)
set.seed(1234)
f <- expression( 2*(3-x) )
g <- expression( 1 )
mod <- snssde1d(drift=f,diffusion=g,M=10000,T=10)
## Monte-Carlo statistics of 5000 trajectory of X(t) at final time T of 'mod'
summary(mod)
kurtosis(mod)
skewness(mod)
mean(mod)
Median(mod)
Mode(mod)
moment(mod,order=4)
cv(mod)
bconfint(mod,level = 0.95) ## of mean

