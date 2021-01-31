# Projection Pursuit Method proposed by Friedman and Tukey in 1973.
# Projection Pursuit paper could be found on Friedman's web page
# Projection Pursuit aims to visualize data in one or two dimension
# Here is an example of one-dim example

# Rosenbrock Algorithm is to compute minimum of a function where derivative is not needed

library(pracma)

Rosenbrock <- function(x, f, epsilon = 10^-2, beta1 = 1.1, beta2 = -0.7){
  # f: the aim function
  # x: starting point
  # epsilon: tolerance
  # beta1, beta2: parameters to adjust step sizes
  n <- length(x)
  # E: directions
  E <- diag(n)
  Z <- matrix(0, nrow = n, ncol = n+1)
  Z[,1] <- x
  delta <- rep(1, n)
  
  repeat{
    for (i in 1:n) {
      if(f(Z[,i] + delta[i] * E[,i]) < f(Z[,i])){
        Z[,i+1] <- Z[,i] + delta[i] * E[,i]
        delta[i] <- delta[i] * beta1
      }
      else{
        Z[,i+1] <- Z[,i]
        delta[i] <- delta[i] * beta2
      }
    }
    
    fn <- f(Z[,n+1])
    f1 <- f(Z[,1])
    f2 <- f(Z[,2])
    
    if(fn < f2){
      Z[,1] <- Z[,n+1]
    }
    else{
      if( fn < f1){
        if(sqrt(sum((Z[,n+1] - Z[,1])^2))<epsilon){
          break
        }
        else{
          Lambda <- solve(E) %*% (Z[,n+1] - Z[,1])
          A <- matrix(0, nrow = n, ncol = n)
          for (j in 1:(n-1)) {
            A[,j] <- (abs(Lambda[j])<10^-3) * E[,j] + (abs(Lambda[j])>=10^-3) * rowSums(E[,j:n])
          }
          A[,n] <- E[,n]
          E <- gramSchmidt(A)$Q
          Z[,1] <- Z[,n+1]
          delta <- rep(1, n)
        }
      }
      else{
        if(max(abs(delta))<epsilon){
          break
        }
        else{
          Z[,1] <- Z[,n+1]
        }
      }
    }
  }
  return(Z[,n+1])
}

# one-dimensional f value
f <- function(r){
  return(R-r)
}

# Projection Pursuit is to minimize S * D

S <- function(X, k, p = 0.05){
  # X: high-dimensional data
  # k: direction
  # p: fraction of p of the points at each of the extremes are omitted
  k <- matrix(k, ncol = 1)
  Xk <- order(X %*% k, decreasing = T)
  n <- nrow(X)
  n1 <- floor(n * p)
  n2 <- floor(n * (1-p))
  
  X_trim <- Xk[n1:n2]
  X_bar <- sum(X_trim)/(n2 - n1)
  s <- sqrt(sum((X_trim - X_bar)^2)/(n2 - n1))
  
  return(s)
}


D <- function(X, k, R){
  # X: Data
  # k: direction
  # R: f values that are larger than R are zet to be zero
  k <- matrix(k, ncol = 1)
  Xk <- X %*% k
  n <- nrow(X)
  A <- matrix(Xk, ncol = 1) %*% matrix(1, nrow = 1, ncol = n)
  A <- abs(A - t(A))
  d <- sum(f(A) * (0+A<R))
  return(d)
}


# I: product of S and D
I <- function(k){
  i <- S(X, k, p=0.05) * D(X, k, R)
  return(i)
}

# In the original paper, transformation is used due to the constraint of \sum_k^2=1
# Solid Angle Transform: E^n-1 (0,1) -> S^n
SAT <- function(eta){
  # check the range
  if(any(eta>=1|eta<=0)){
    stop("eta is not in the range of E^n(0,1)")
  }
  
  n <- length(eta) + 1
  x <- matrix(data = NA, nrow = n)
  
  temp <- 1
  
  if(n%%2==0){
    if(n>2){
      for (i in 1:(n/2 - 1)) {
      x[2*i] <- temp * cos(asin(eta[2*i]^(1/(n - 2*i) ) ) ) * sin(2 * pi * eta[2*i-1])
      x[2*i - 1] <- x[2*i] / tan(2*pi*eta[2*i-1])
      temp <- temp * eta[2 * i]^(1/(n - 2 * i))
    }
    x[n] <- temp * sin(2*pi*eta[n-1])
    x[n-1] <- x[n] / tan(2*pi*eta[n-1])
    }
    else{
      x[2] <- sin(2*pi*eta)
      x[1] <- cos(2*pi*eta)
    }
  }
  else{
    temp <- sqrt(2 * eta[n-1] - eta[n-1]^2)
    if(n>3){
      for (i in 1:((n-3)/2)) {
        x[2*i] <- temp * cos(asin(eta[2*i]^(1/(n - 2*i) ) ) ) * sin(2 * pi * eta[2*i-1])
        x[2*i - 1] <- x[2*i] / tan(2*pi*eta[2*i-1])
        temp <- temp * eta[2 * i]^(1/(n - 2 * i))
    }
    x[n] <- temp * (eta[n-1] - 1)/sqrt(2 * eta[n-1] - eta[n-1]^2)
    x[n-1] <- temp * sin(2*pi*eta[n-2])
    x[n-2] <- x[n-1] / tan(2*pi*eta[n-2])
    }
    else{
      x[3] <- eta[1] - 1
      x[2] <- sqrt(2 * eta[1] - eta[1]^2) * sin(2*pi*eta[2])
      x[1] <-  sqrt(2 * eta[1] - eta[1]^2)* cos(2 * pi * eta[2])
    }
  }
  return(x)
}
# Transform E^n (-inf, inf) -> E^n (0,1)
Tan_transform <- function(z){
  x <- atan(z) /pi + 0.5
  return(x)
}



# Another classic transformation: transform E^n-1 (-inf, inf) -> S^n
Transform_ <- function(eta){
  k <- 4/(4 + sum(eta^2))
  eta_new <- c(eta, 2)
  x <- k * eta_new
  n <- length(x)
  x[n] <- x[n] - 1
  return(x)
}




# The Function that needs to be minimized
# Aim1 uses Transform_ transformation
Aim1 <- function(x){
  k = Transform_(x)
  y = -S(X, k, p=0.05) * D(X, k, R)
  return(y)
}
# Aim2 uses SAT transformation
Aim2 <- function(x){
  k = SAT(Tan_transform(x))
  y = -S(X, k, p=0.05) * D(X, k, R)
  return(y)
}




# SImulation
library(MASS)


m <- 4
N <- 300
X <- mvrnorm(n = N, mu = seq(m), Sigma = diag(m))
R <- max(floor(log(N)), 1)


for (i in 1:10) {
  initial <- rnorm(m-1, mean = 0, sd = 1)
  k1 <- Rosenbrock(initial, f = Aim1)
  k2 <- Rosenbrock(initial, f = Aim2)
  message("Starting point:",initial)
  # Because Aim = - S * D, a minus is added before AIm1 and Aim2
  message("Convergence Direction in Method1: ", Transform_(k1))
  message("Initial function value of Method1: ", -Aim1(initial), " Convergence function value of Method1: ", -Aim1(k1))
  message("Convergence Direction in Method2: ", SAT(Tan_transform(k2)))
  message("Initial function value of Method2: ", -Aim2(initial), " Convergence function value of Method1: ", -Aim2(k2))
}

