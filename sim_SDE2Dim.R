# Simulation of 2-D Bridge SDE’s

library(Sim.DiffProc)

## dX(t) = 4*(-1-X(t)) dt + 0.2 dW1(t)
## dY(t) = X(t) dt + 0 dW2(t)
## x01 = 0 , y01 = 0
## x02 = 0, y02 = 0
## W1(t) and W2(t) two correlated Brownian motion 
## with matrix Sigma=matrix(c(1,0.7,0.7,1),nrow=2)
set.seed(1234)
fx <- expression(4*(-1-x) , x)
gx <- expression(0.2 , 0)
Sigma= matrix(c(1,0.7,0.7,1),nrow=2)
res <- bridgesde2d(drift=fx,diffusion=gx,Dt=0.005,M=500,corr=Sigma)
res
summary(res) ## Monte-Carlo statistics at time T/2=2.5
summary(res,at=1) ## Monte-Carlo statistics at time 1
summary(res,at=4) ## Monte-Carlo statistics at time 4

plot(res)
