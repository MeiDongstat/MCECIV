###########################CODE for Estimation###########################

logit <- function(prob){ log(prob) - log(1-prob)}

expit <- function(logodds){ 1/(1+exp(-logodds))}

link_identity <- list(
  name  = "identity",
  fn    = function(eta) eta
)

link_tanh <- list(
  name  = "tanh",
  fn    = function(eta) tanh(eta)
)

## function for estimating \beta_ipw when the IV Z is continuous
beta.est <- function(x.delta.z, z, d, pred.mu.d){
  ### solve estimating equation for \beta
  beta.objective = function(beta){
    obj = t(x.delta.z) %*% ((d - pred.mu.d) /(pred.mu.d * (1-pred.mu.d))*z - x.delta.z %*% beta)
    return(sum(obj^2))
  }
  startpars = rep(0,dim(x.delta.z)[2])
  opt = optim(startpars, beta.objective)
  return(opt)
}

## function for estimating \beta_ipw when the IV Z is binary
beta.bounded.est <- function(x.delta.z, z, d, pred.mu.d){
  beta.objective = function(beta){
    obj = t(x.delta.z) %*% ((d - pred.mu.d) /(pred.mu.d * (1-pred.mu.d))*z - tanh(x.delta.z %*% beta))
    return(sum(obj^2))
  }
  startpars = rep(0,dim(x.delta.z)[2])
  opt = optim(startpars, beta.objective)
  return(opt)
}

## function for estimating \alpha_dr when the outcome Y is continuous
alpha.dr.Est <- function(x.delta, pred.pi, pred.mu.d, pred.mu.y, z, d, y, link){
  ### solve estimating equation for doubly robust \alpha
  alpha.objective = function(alpha){
    obj = t(x.delta) %*% ( (z - pred.pi)*(y - pred.mu.y - link$fn(x.delta %*% alpha) * (d - pred.mu.d)) )
    return(sum(obj^2))
  }
  startpars = rep(0,dim(x.delta)[2])
  opt = optim(startpars, alpha.objective)
  return(opt)
}

## function for estimating \alpha_dr when the outcome Y is binary
alpha.dr.bounded.Est <- function(x.delta, pred.pi, pred.mu.d, pred.mu.y, cov.zd, z, d, y, link){
  alpha.objective = function(alpha){
    obj = t(cbind(1/cov.zd, x.delta[,-1])) %*% ((z - pred.pi)*(y - pred.mu.y - link$fn(x.delta %*% alpha) * (d - pred.mu.d)))
    return(sum(obj^2))
  }
  startpars = rep(0,dim(x.delta)[2])
  opt = optim(startpars, alpha.objective)
  return(opt)
}

## function for estimating alpha_1 for $\Delta_1$ in supplementary Material S.2
alpha1.est <- function(x.delta, pred.mu.d, pred.mu.y, z, d, y, link){
  alpha1.objective = function(alpha){
    obj = t(x.delta) %*% (z * (y - pred.mu.y - link$fn(x.delta %*% alpha) * (d - pred.mu.d)))
    return(sum(obj^2))
  }
  startpars = rep(0,dim(x.delta)[2])
  opt = optim(startpars, alpha1.objective)
  return(opt)
}

## function for estimating alpha_2 for $\Delta_{b-2}$ when the outcome Y is binary in supplementary Material S.2  
alpha2.bouned.est <- function(x.delta, x.delta.z, pred.pi, pred.mu.d, beta.hat, z, y, link){
  
  delta.z <- x.delta.z %*% beta.hat
  alpha2.objective = function(alpha){
    obj = t(x.delta) %*% ( (z - pred.pi)*y/(delta.z * pred.mu.d * (1 - pred.mu.d)) - link$fn(x.delta %*% alpha))
    return(sum(obj^2))
  }
  startpars = rep(0,dim(x.delta)[2])
  opt = optim(startpars, alpha2.objective)
  return(opt)
}

## function for estimating alpha_3 for $\Delta_3$ in supplementary Material S.2
alpha3.est <- function(x.delta, pred.pi, z, d, y, link){
  alpha3.objective = function(alpha){
    obj = t(x.delta) %*% ( (z - pred.pi)*(y - link$fn(x.delta %*% alpha) * d) )
    return(sum(obj^2))
  }
  startpars = rep(0,dim(x.delta)[2])
  opt = optim(startpars, alpha3.objective)
  return(opt)
}

##function for estimating ATE when Z is continuous
estimation <- function(x.delta, x.delta.z, x.mu.d, x.mu.y, x.pi, z, d, y, link_delta = c("identity", "tanh")){
  
  ##this matches the link function for \delta(X; \alpha)
  link_delta <- match.arg(link_delta)
  link <- if (link_delta == "identity") link_identity else link_tanh
  
  ### Estimate E(D \mid X) using mle
  mod.d <- glm(d ~ x.mu.d - 1, family = binomial())
  iota.est <- mod.d$coefficients
  pred.mu.d <- predict(mod.d, type = "response")
  
  ### Estimate E(Y \mid X) using mle
  if(link_delta == "identity"){ ## when Y is continuous, we use identity link
    mod.y <- glm(y ~ x.mu.y - 1, family = gaussian())
    pred.mu.y <- predict(mod.y, type = "response")
  } else { ## when Y is binary, we can use tanh() function or other functions that are bounded in (-1, 1)
    mod.y <- glm(y ~ x.mu.y - 1, family = binomial())
    pred.mu.y <- predict(mod.y, type="response")
  }
  
  ### Estimate E(Z \mid X) using mle  
  mod.z <- glm(z ~ x.pi - 1, family = gaussian())
  zeta.hat = mod.z$coefficients
  pred.pi <- predict(mod.z, type = "response")
  
  ### Estimate beta indexing \delta^Z(X; \beta)
  beta.hat <- beta.est(x.delta.z, z, d, pred.mu.d)$par
  cov.zd <- (x.delta.z %*% beta.hat) * pred.mu.d * (1 - pred.mu.d)
  
  ## estimate \Delta_1
  alpha1 <- alpha1.est(x.delta, pred.mu.d, pred.mu.y, z, d, y, link)$par
  Delta1 <- mean(link$fn(x.delta %*% alpha1))
  
  ## estimate \Delta_2
  Delta2 <- mean( (z - pred.pi)*y/cov.zd)
  alpha2.b <- alpha2.bouned.est(x.delta, x.delta.z, pred.pi, pred.mu.d, beta.hat, z, y, link)$par
  
  ## estimate \Delta_{b-2}  
  Delta2.b <- mean(link$fn(x.delta %*% alpha2.b))
  
  ## estimate \Delta_3
  alpha3 <- alpha3.est(x.delta, pred.pi, z, d, y, link)$par
  Delta3 <- mean(link$fn(as.vector(x.delta %*% alpha3)))

  ## estimate triply robust estimator
  if(link_delta == "identity"){
    
    ### doubly robust estimation for \alpha indexing \delta(X; \alpha)
    alpha.dr.hat <- alpha.dr.Est(x.delta, pred.pi, pred.mu.d, pred.mu.y, z, d, y, link)$par
    
    ## estimate \Delta_{tr}
    Delta.mr = mean( (z-pred.pi) * (y - pred.mu.y - as.vector(x.delta %*% alpha.dr.hat)* (d - pred.mu.d))/ cov.zd + x.delta %*% alpha.dr.hat)
  } else {
    ## estimate \Delta_{b-tr}
    alpha.dr.hat <- alpha.dr.bounded.Est(x.delta, pred.pi, pred.mu.d, pred.mu.y, cov.zd, z, d, y, link)$par
    Delta.mr = mean(link$fn(x.delta %*% alpha.dr.hat))
  }

  return(list(beta.hat=beta.hat,iota.hat=iota.est, zeta.hat = zeta.hat,
              alpha.dr.hat = alpha.dr.hat, Delta1 = Delta1, Delta2 = Delta2, Delta2.b = Delta2.b, 
              Delta3 = Delta3, Delta.mr = Delta.mr))
}

##estimate ATE when Z is binary
estimation_binZ <- function(x.delta, x.delta.z, x.mu.d, x.mu.y, x.pi, z, d, y, link_delta = c("identity", "tanh")){
  
  ##this matches the link function for \delta(X; \alpha)
  link_delta <- match.arg(link_delta)
  link <- if (link_delta == "identity") link_identity else link_tanh
  
  ### Estimate E(D \mid X) using mle
  mod.d <- glm(d ~ x.mu.d - 1, family = binomial())
  iota.est <- mod.d$coefficients
  pred.mu.d <- predict(mod.d, type = "response")
  ### Estimate E(Y \mid X) using mle
  if(link_delta == "identity"){ ## when Y is continuous, we use identity link
    mod.y <- glm(y ~ x.mu.y - 1, family = gaussian())
    pred.mu.y <- predict(mod.y, type = "response")
  } else { ## when Y is binary, we can use tanh() function or other functions that are bounded in (-1, 1)
    mod.y <- glm(y ~ x.mu.y - 1, family = binomial())
    pred.mu.y <- predict(mod.y, type="response")
  }
  
  ### Estimate E(Z \mid X) using mle  
  mod.z <- glm(z ~ x.pi - 1, family = binomial())
  zeta.hat = mod.z$coefficients
  pred.pi <- predict(mod.z, type = "response")
  
  ### Estimate beta indexing \delta^Z(X; \beta)
  beta.hat <- beta.bounded.est(x.delta.z, z, d, pred.mu.d)$par
  cov.zd <- tanh(x.delta.z %*% beta.hat) * pred.mu.d * (1 - pred.mu.d)
  
  ## Estimate \Delta_1
  alpha1 <- alpha1.est(x.delta, pred.mu.d, pred.mu.y, z, d, y, link)$par
  Delta1 <- mean(link$fn(x.delta %*% alpha1))
  
  ## Estimate \Delta_2
  Delta2 <- mean( (z - pred.pi)*y/cov.zd)
  alpha2.b <- alpha2.bouned.est(x.delta, x.delta.z, pred.pi, pred.mu.d, beta.hat, z, y, link)$par
  Delta2.b <- mean(link$fn(x.delta %*% alpha2.b))
  
  ## Estimate \Delta_3
  alpha3 <- alpha3.est(x.delta, pred.pi, z, d, y, link)$par
  Delta3 <- mean(link$fn(as.vector(x.delta %*% alpha3)))
  
  ## Estimate triply robust estimator
  if(link_delta == "identity"){
    
    ### doubly robust estimation for \alpha indexing \delta(X; \alpha)
    alpha.dr.hat <- alpha.dr.Est(x.delta, pred.pi, pred.mu.d, pred.mu.y, z, d, y, link)$par
    
    ## Estimate \Delta_{tr}
    Delta.mr = mean( (z-pred.pi) * (y - pred.mu.y - as.vector(x.delta %*% alpha.dr.hat)* (d - pred.mu.d))/ cov.zd + x.delta %*% alpha.dr.hat)
  } else {
    ## Estimate \Delta_{b-tr}
    alpha.dr.hat <- alpha.dr.bounded.Est(x.delta, pred.pi, pred.mu.d, pred.mu.y, cov.zd, z, d, y, link)$par
    Delta.mr = mean(link$fn(x.delta %*% alpha.dr.hat))
  }
  
  return(list(beta.hat=beta.hat,iota.hat=iota.est, zeta.hat = zeta.hat, 
              alpha.dr.hat = alpha.dr.hat, Delta1 = Delta1, Delta2 = Delta2, Delta2.b = Delta2.b, 
              Delta3 = Delta3, Delta.mr = Delta.mr))
}
