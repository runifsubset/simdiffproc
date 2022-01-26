library(Sim.DiffProc)

op <- par(mfrow = c(2, 2))

## Brownian motion
title<- 'Brownian motion'
set.seed(1234)
X <- BM(M = 100)
plot(X,plot.type="single",main=title)
lines(as.vector(time(X)),rowMeans(X),col="red")

## Brownian bridge
title <- 'Brownian bridge'
set.seed(1234)
X <- BB(M =100)
plot(X,plot.type="single",main=title)
lines(as.vector(time(X)),rowMeans(X),col="red")

## Geometric Brownian motion
title <- 'Geometric Brownian motion'
set.seed(1234)
X <- GBM(M = 100)
plot(X,plot.type="single",main=title)
lines(as.vector(time(X)),rowMeans(X),col="red")

## Arithmetic Brownian motion
title <- 'Arithmetic Brownian motion'
set.seed(1234)
X <- ABM(M = 100)
plot(X,plot.type="single",main=title)
lines(as.vector(time(X)),rowMeans(X),col="red")
par(op)
