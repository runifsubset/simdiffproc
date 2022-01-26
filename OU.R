# dXt = (θ1 − θ2Xt)dt + θ3dWt .
# Ornstein-Uhlenbeck process

set.seed(123)
N <- 1000
T <- 1
x <-10
theta <- c(0, 5, 3.5)
Dt <- 1/N
Y <- numeric(N+1)
Y[1] <-x
Z  <- rnorm(N)
  for (i in 1:N) {
    Y[i+1] <- Y[i] + (theta[1]-theta[2]*Y[i])*Dt + theta[3]*sqrt(Dt)*Z[i]
  }
Y <- ts(Y, start = 0, deltat = 1/N)
plot(Y, main="Ornstein-Uhlenbeck process")


## Y (k) = Y (1) · (1 − θ2∆t)k−1 + sum (from = k−1, to = Xj =1 of (θ3 ∗ √∆t ∗ Z(j)) ∗ (1 − θ2∆t)(k−j−1)

set.seed(123)
theta <- c(0, 5, 3.5)
N <- 1000
T <- 1
x <-10
Z <- rnorm (N)
Dt <- 1/N
A <- theta [3] * sqrt (Dt)*Z
P <- (1- theta [2] *Dt )^(0:(N -1))
X0 <- x
X <- sapply (2:( N+1) , function (x) X0*(1- theta [2] *Dt )^(x -1) +
               A [1:(x -1)] %*% P[(x -1):1])
Y <- ts(c(X0 ,X), start =0, deltat =1/N)
plot (Y)


## Or : from given values of Theta

OU <- function (a,b, x, N =1000){
   Y <- numeric (N +1)
   Y [1] <- x
   Z <- rnorm (N)
   Dt <- 1/N
   for (i in 1:N){
     Y[i +1] <- Y[i] - a*Y[i]*Dt + b* sqrt (Dt)*Z[i]
   }
     invisible (Y)
}


OU.vec <- function (a, b, x, N =1000)
  {
   Dt <- 1/N
   Z <- rnorm(N)
   A <- b*sqrt(Dt)*Z
   P <- (1-a*Dt)^(0:(N -1))
   X0 <- x
   X <- c(X0 , sapply(2:(N+1),function (x)X0*(1-a*Dt )^(x -1)+
       sum (A[1:(x -1)] * P[(x -1):1])))
   invisible(X)
}

## Using the system.time function, we test the two implementations

set.seed (123)
system.time(OU(10,5,3.5))
system.time(OU.vec(10,5,3.5))
