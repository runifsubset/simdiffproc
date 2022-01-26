## Maximum Pseudo-Likelihood Estimation of 1-D SDE

## The (S3) generic function "fitsde" of estimate drift 
## and diffusion parameters by the method of
## maximum pseudo-likelihood 
## of the 1-dim stochastic differential equation.

##### Example 1:
## Modele GBM (BS)
## dX(t) = theta1 * X(t) * dt + theta2 * x * dW(t)
## Simulation of data

set.seed(1234)
X <- GBM(N =1000,theta=4,sigma=1)
## Estimation: true theta=c(4,1)
fx <- expression(theta[1]*x)
gx <- expression(theta[2]*x)
fres <- fitsde(data=X,drift=fx,diffusion=gx,start = list(theta1=1,theta2=1),lower=c(0,0))
fres
summary(fres)
coef(fres)
logLik(fres)
AIC(fres)
BIC(fres)
vcov(fres)
confint(fres,level=0.95)


## Approximate densities and random generation 
## for first passage time in 1-D SDE

## Example 1: Ito SDE
## dX(t) = -4*X(t) *dt + 0.5*dW(t)
## S(t) = 0 (constant boundary)
set.seed(1234)
# SDE 1d
f <- expression( -4*x )
g <- expression( 0.5 )
mod <- snssde1d(drift=f,diffusion=g,x0=2,M=1000)
# boundary
St <- expression(0)
# random
out <- fptsde1d(mod, boundary=St)
out
summary(out)
# density approximate
den <- dfptsde1d(out)
den
plot(den)


## Example 2: Stratonovich SDE
## dX(t) = 0.5*X(t)*t *dt + sqrt(1+X(t)^2) o dW(t)
## S(t) = -0.5*sqrt(t) + exp(t^2) (time-dependent boundary)
set.seed(1234)

# SDE 1d
f <- expression( 0.5*x*t )
g <- expression( sqrt(1+x^2) )
mod2 <- snssde1d(drift=f,diffusion=g,x0=2,M=1000,type="srt")
# boundary
St <- expression(-0.5*sqrt(t)+exp(t^2))
# random
out2 <- fptsde1d(mod2,boundary=St)
out2
summary(out2)
# density approximate
plot(dfptsde1d(out2,bw='ucv'))
## Example 3: fptsde1d vs fptdApproximate
## Not run:
f <- expression( -0.5*x+0.5*5 )
g <- expression( 1 )
St <- expression(5+0.25*sin(2*pi*t))
mod <- snssde1d(drift=f,diffusion=g,boundary=St,x0=3,T=10,N=10^4,M =10000)
mod
# random
out3 <- fptsde1d(mod,boundary=St)
out3
summary(out3)
# density approximate:
library("fptdApprox")
# Under `fptdApprox':
# Define the diffusion process and give its transitional density:
OU <- diffproc(c("alpha*x + beta","sigma^2",
                 "dnorm((x-(y*exp(alpha*(t-s)) - beta*(1 - exp(alpha*(t-s)))/alpha))/
(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha))),0,1)/
(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))",
                 "pnorm(x, y*exp(alpha*(t-s)) - beta*(1 - exp(alpha*(t-s)))/alpha,
sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))"))
# Approximate the first passgage time density for OU, starting in X_0 = 3
# passing through 5+0.25*sin(2*pi*t) on the time interval [0,10]:
res <- Approx.fpt.density(OU, 0, 10, 3,"5+0.25*sin(2*pi*t)", list(alpha=-0.5,beta=0.5*5,sigma=1))
##
plot(dfptsde1d(out3,bw='ucv'),main = 'fptsde1d vs fptdApproximate')
lines(res$y~res$x, type = 'l',lwd=2)
legend('topright', lty = c('solid', 'dashed'), col = c(1, 2),
       legend = c('fptdApproximate', 'fptsde1d'), lwd = 2, bty = 'n')
