## This code contains functions used to generate simulated dataset
## data_gen() is used to generate data in Setting I in the simulation studies
## data_gen_binZ() is used to generate data in Setting II in the simulation studies

logit <- function(prob){ log(prob) - log(1-prob)}

expit <- function(logodds){ 1/(1+exp(-logodds))}

## we simulate Z using a mixture normal; then simulate D given Z and X, and U
## we then simulate Y given Z and X, with additive term of U
## we write U inside the expit() function of P(D=1| X, Z, U), therefore, there is Z-U interaction in predicting D, 
## but we do not have d-U interaction in predicting Y(d)
## Assumption 5d holds in our simulation
data_gen <- function(alpha.true, beta.true, iota.true, theta.true, zeta.true, n, sigma=1, k=0.1, Y.type = "binary",seed){
  set.seed(seed)
  
  ##simulate measured confounder
  intc <- rep(1,n)  #intercept term
  x1   <- rnorm(n, 0, 1)
  x2.ind       = rbinom(n, 1, 0.5)
  # x2           = ((-1)^x2.ind) * runif(n, 0.5, 1)
  x2 <- runif(n, 0, 1)
  x3 <- rbinom(n, 1, 0.4)
  x4 <- rnorm(n, 0, 1)
  X <- cbind(intc, x2, x3)
  
  ## simulate unmeasured confounder
  U <- runif(n, 0, 1)
  kk = 0.8
  B_U <- (1/kk) * log(sin(pi * kk * U)/sin(pi * kk * (1 - U)))
  
  ## E[D|X]=expit(\iota^T X)
  mu.d = expit(X %*% iota.true)
  
  ## E[Z|X] = \zeta^T X
  pi.x = X %*% zeta.true
  
  ## \delta^Z(X) = E[Z|D=1, X] -E[Z|D=0, X] = beta^T X
  delta.z.true = X %*% beta.true
  
  ##E[Z|D=0, X]
  pi0 = pi.x - mu.d * delta.z.true
  
  ##E[Z|D=1, X]
  pi1 = pi.x + (1 - mu.d) * delta.z.true
  
  ## simulate Z using a mixture normal distribution
  D <- rbinom(n, 1, mu.d)
  Z <- rnorm(n, mean=ifelse(D ==1, pi1, pi0), sd = sigma)
  
  mu.d.zx <- expit(X %*% iota.true + delta.z.true/sigma^2 * Z - (pi1^2 - pi0^2)*0.5/sigma^2)
  
  ## E[D|Z, X] derived from bayes rule
  p.d.true.u <- expit((1/kk) * (X %*% iota.true + delta.z.true/sigma^2 * Z - (pi1^2 - pi0^2)*0.5/sigma^2) + B_U)
  
  D <- rbinom(n, 1, p.d.true.u)
  
  if(Y.type == "binary"){
    ## E[Y|X] = expit(theta^T X)
    mu.y <- expit(X %*% theta.true)
    
    atanhrd.delta.true <- X %*% alpha.true
    delta.true <- tanh(atanhrd.delta.true)
    
    ## E[Y|Z, X] = \delta(X)E[D|Z, X] + E[Y|X] - \delta(X)E[D|X]
    mu.y.zx <- delta.true * mu.d.zx + mu.y - delta.true * mu.d
    
    mu.y.true.u <- mu.y.zx + k * (2*U -1)
    
    if (any(!is.finite(mu.y.true.u) | mu.y.true.u < 0 | mu.y.true.u > 1))
      stop("mu.y.true.u out of bounds")
    
    Y <- rbinom(n, 1, mu.y.true.u)
    
  } else{
    mu.y <- X %*% iota.true
    
    delta.true <- X %*% alpha.true
    
    ## E[Y|Z, X] = \delta(X)E[D|Z, X] + E[Y|X] - \delta(X)E[D|X]
    mu.y.zx <- delta.true * mu.d.zx + mu.y - delta.true * mu.d
    
    mu.y.true.u <- mu.y.zx + 2*k * (2*U -1)
    
    Y = rnorm(n, mean = mu.y.true.u, sd=0.5)
  }
  
  dat <- data.frame(cbind(int=intc,X=X[, -1], U=U, Z=Z, D=D, Y=Y))
  colnames(dat) <- c("Intercept","X1", "X2","U", "Z", "D","Y")
  return(list(dat=dat, mu.d.zx=mu.d.zx, mu.y.zx=mu.y.zx, delta.true=delta.true))
}

data_gen_binZ <- function(alpha.true, beta.true, iota.true, theta.true, zeta.true,
                          n, k=0.1, Y.type = "binary", seed) {
  set.seed(seed)
  
  ## measured confounders 
  intc <- rep(1, n)
  x1   <- rnorm(n, 0, 1)
  x2.ind <- rbinom(n, 1, 0.5)
  # x2 <- ((-1)^x2.ind) * runif(n, 0.5, 1)
  x2 <- runif(n, 0, 1)
  x3 <- rbinom(n, 1, 0.4)
  X  <- cbind(intc, x2, x3)
  
  ## unmeasured confounder
  U <- runif(n, 0, 1)
  kk = 0.8
  B_U <- (1/kk) * log(sin(pi * kk * U)/sin(pi * kk * (1 - U)))
  
  ## E[D|X] = expit(iota^T X)
  mu.d <- expit(X %*% iota.true)
  
  ## E[Z|X] = expit(\zeta^T X)
  pi.x <- expit(X %*% zeta.true)
  
  ## delta_Z(X) = tanh(beta^T X)
  delta.z.true <- tanh(X %*% beta.true)
  
  ## Construct pi0, pi1 so that:
  ## E[Z|X] = mu.d*pi1 + (1-mu.d)*pi0 = pi.x, and pi1 - pi0 = delta.z.true
  pi0 <- pi.x - mu.d * delta.z.true
  pi1 <- pi.x + (1 - mu.d) * delta.z.true
  
  S <- rbinom(n, 1, mu.d)
  
  Z <- rbinom(n, 1, pi.x) 
  ## E[D | Z, X] via Bayes for binary Z
  # num1 <- mu.d * pi1
  # den1 <- mu.d * pi1 + (1 - mu.d) * pi0
  # num0 <- mu.d * (1 - pi1)
  # den0 <- mu.d * (1 - pi1) + (1 - mu.d) * (1 - pi0)
  # mu.d.zx <- ifelse(Z == 1, num1 / den1, num0 / den0)
  c0 <- log( (1-pi1)/(1-pi0))
  c1 <- log((pi1 * (1-pi0))/(pi0 * (1-pi1)))
  
  mu.d.zx <- expit(X %*% iota.true + c0 + c1 * Z)
  p.d.true.u <- expit((1/kk) * (X %*% iota.true + c0 + c1 * Z) + B_U)
  
  D <- rbinom(n, 1, p.d.true.u)
  
  if (Y.type == "binary") {
    ## E[Y|X] = expit(theta^T X)
    mu.y <- expit(X %*% theta.true)
    
    ## delta(X) = tanh(alpha^T X)
    atanhrd.delta.true <- X %*% alpha.true
    delta.true <- tanh(atanhrd.delta.true)
    
    ## E[Y|Z, X] = delta(X) E[D|Z, X] + E[Y|X] - delta(X) E[D|X]
    mu.y.zx <- delta.true * mu.d.zx + mu.y - delta.true * mu.d
    
    mu.y.true.u <- mu.y.zx + k * (2 * U - 1)
    if (any(!is.finite(mu.y.true.u) | mu.y.true.u < 0 | mu.y.true.u > 1))
      stop("mu.y.true.u out of bounds")
    
    Y <- rbinom(n, 1, mu.y.true.u)
    
  } else {
    ## continuous Y case 
    mu.y <- X %*% iota.true
    delta.true <- X %*% alpha.true
    mu.y.zx <- delta.true * mu.d.zx + mu.y - delta.true * mu.d
    mu.y.true.u <- mu.y.zx + 2 * k * (2 * U - 1)
    Y <- rnorm(n, mean = mu.y.true.u, sd = 0.5)
  }
  
  dat <- data.frame(cbind(int = intc, X = X[, -1], U = U, Z = Z, D = D, Y = Y))
  colnames(dat) <- c("Intercept", "X1", "X2", "U", "Z", "D", "Y")
  
  return(list(dat = dat,
              mu.d.zx = mu.d.zx,
              mu.y.zx = mu.y.zx,
              delta.true = delta.true,
              delta.z.true = delta.z.true,
              pi0 = pi0, pi1 = pi1, pi.x = pi.x))
}