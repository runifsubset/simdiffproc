library(Sim.DiffProc)

set.seed(1234, kind = "L'Ecuyer-CMRG")
theta = 0.5
f <- expression( (0.5*theta^2*x) )
g <- expression( theta*x )
mod1 <- snssde1d(drift=f,diffusion=g,x0=10,M=1000,type="ito",T=15) # Using Ito
mod2 <- snssde1d(drift=f,diffusion=g,x0=10,M=1000,type="str",T=15) # Using Stratonovich 

mem = MEM.sde(drift = f,diffusion = g)
TEX.sde(mem)
TEX.sde(object=c(drift = f, diffusion = g))


AppendI = dsde1d(mod1, at =1)
AppendS = dsde1d(mod2, at =1)


par(mfrow = c(1, 2))

## Ito
  plot(mod1,ylab=expression(X^mod1))
lines(time(mod1),apply(mod1$X,1,mean),col=2,lwd=2)
lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(mod1),apply(mod1$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of", 95,"% confidence")),inset = .01,col=c(2,4),lwd=2,cex=0.8)
## Stratonovich
  plot(mod2,ylab=expression(X^mod2))
lines(time(mod2),apply(mod2$X,1,mean),col=2,lwd=2)
lines(time(mod2),apply(mod2$X,1,bconfint,level=0.95)[1,],col=4,lwd=2)
lines(time(mod2),apply(mod2$X,1,bconfint,level=0.95)[2,],col=4,lwd=2)
legend("topleft",c("mean path",paste("bound of",95,"% confidence")),col=c(2,4),inset =.01,lwd=2,cex=0.8)

par(mfrow = c(1, 1))

