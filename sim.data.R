# Functions to simulate data from a LMM with one longitudinal covariate
# - balanced design: same number of measurements for each individual
# - linear or quadratic functional form
# - random intercept and slope for both outcome and covariate



sim.lin.seq <- function(N, J, beta, sig, D, timerange = c(0,5),...){
  # N: number of individuals
  # J: number of repeated measurements per individual
  # beta: list(beta.y = c(Intercept, x, time), beta.x = c(Intercept, time))
  # sig: c(sig.y, sig.x)
  # D: list(D.y, D.x)
  # timerange: maximum range of time of observations
  
  id <- rep(1:N, each = J)
  time <- unlist(lapply(1:N, function(i){
    sort(runif(n = J, min = timerange[1], max = timerange[2]))
  }))
  
  Z <- cbind(Intercept = 1, time)
  
  b.x <- mvrnorm(N, c(0,0), D[[2]])
  x <- Z %*% beta[[2]] + rowSums(Z * b.x[match(id, unique(id)), ]) + 
    rnorm(N*J, 0, sig[2])
  
  b.y <- mvrnorm(N, c(0,0), D[[1]])
  y <- cbind(1, x,  Z[,2]) %*% beta[[1]] + 
    rowSums(Z * b.y[match(id, unique(id)), ]) + rnorm(N*J, 0, sig[1])
  
  DF <- data.frame(id = id, x=x, y=y, time = time)
  
  return(DF)
}


sim.lin.seq.b <- function(N, J, beta, sig, D, timerange = c(0,5),...){
  # simulate from sequential model, with linear relationship to x, and 
  # correlated random offects
  # N: number of individuals
  # J: number of repeated measurements per individual
  # beta: list(beta.y = c(Intercept, x, time), beta.x = c(Intercept, time))
  # sig: c(sig.y, sig.x)
  # D: covariance matrix of random effects
  # timerange: maximum range of time of observations
  
  id <- rep(1:N, each = J)
  time <- unlist(lapply(1:N, function(i){
    sort(runif(n = J, min = timerange[1], max = timerange[2]))
  }))
  
  Z <- cbind(Intercept = 1, time)
  
  b <- mvrnorm(N, c(0,0,0,0), D)
  x <- Z %*% beta[[2]] + rowSums(Z * b[match(id, unique(id)), 3:4]) +
    rnorm(N*J, 0, sig[2])
  
  y <- cbind(1, x,  Z[,2]) %*% beta[[1]] + 
    rowSums(Z * b[match(id, unique(id)), 1:2]) + rnorm(N*J, 0, sig[1])
  
  DF <- data.frame(id = id, x=x, y=y, time = time)
  
  return(DF)
}


sim.qdr.seq <- function(N, J, beta, sig, D, timerange = c(0,5),...){
  # N: number of individuals
  # J: number of repeated measurements per individual
  # beta: list(beta.y = c(Intercept, x, x2, time), beta.x = c(Intercept, time))
  # sig: c(sig.y, sig.x)
  # D: list(D.y, D.x)
  # timerange: maximum range of time of observations
  
  id <- rep(1:N, each = J)
  time <- unlist(lapply(1:N, function(i){
    sort(runif(n = J, min = timerange[1], max = timerange[2]))
  }))
  
  Z <- cbind(Intercept = 1, time)
  
  b.x <- mvrnorm(N, c(0,0), D[[2]])
  x <- Z %*% beta[[2]] + rowSums(Z * b.x[match(id, unique(id)), ]) +
    rnorm(N*J, 0, sig[2])
  
  b.y <- mvrnorm(N, c(0,0), D[[1]])
  y <- cbind(1, x,  x^2, Z[,2]) %*% beta[[1]] + 
    rowSums(Z * b.y[match(id, unique(id)), ]) + rnorm(N*J, 0, sig[1])
  
  DF <- data.frame(id = id, x=x, y=y, time = time)
  
  return(DF)
}

sim.qdr.seq.b <- function(N, J, beta, sig, D, timerange = c(0,5),...){
  # N: number of individuals
  # J: number of repeated measurements per individual
  # beta: list(beta.y = c(Intercept, x, x2, time), beta.x = c(Intercept, time))
  # sig: c(sig.y, sig.x)
  # D: list(D.y, D.x)
  # timerange: maximum range of time of observations
  
  id <- rep(1:N, each = J)
  time <- unlist(lapply(1:N, function(i){
    sort(runif(n = J, min = timerange[1], max = timerange[2]))
  }))
  
  Z <- cbind(Intercept = 1, time)
  
  b <- mvrnorm(N, c(0,0,0,0), D)
  x <- Z %*% beta[[2]] + rowSums(Z * b[match(id, unique(id)), 3:4]) +
    rnorm(N*J, 0, sig[2])
  
  y <- cbind(1, x,  x^2, Z[,2]) %*% beta[[1]] +
    rowSums(Z * b[match(id, unique(id)), 1:2]) + rnorm(N*J, 0, sig[1])
  
  DF <- data.frame(id = id, x=x, y=y, time = time)
  
  return(DF)
}


sim.lin.joint.b <- function(N, J, beta, sig, D, timerange = c(0,5), ...){
  # N: number of individuals
  # J: number of repeated measurements per individual
  # beta: list(beta.y = c(Intercept, time), beta.x = c(Intercept, time))
  # sig: c(sig.y, sig.x)
  # D: covariance matrix of the random effects
  # timerange: maximum range of time of observations
  
  id <- rep(1:N, each = J)
  time <- unlist(lapply(1:N, function(i){
    sort(runif(n = J, min = timerange[1], max = timerange[2]))
  }))
  
  Z <- cbind(Intercept = 1, time)
  
  b <- mvrnorm(N, c(0,0,0,0), D)
  
  x <- Z %*% beta[[2]] + rowSums(Z * b[match(id, unique(id)), 3:4]) +
    rnorm(N*J, 0, sig[2])
  y <- Z %*% beta[[1]] + rowSums(Z * b[match(id, unique(id)), 1:2]) + 
    rnorm(N*J, 0, sig[1])
  
  
  DF <- data.frame(id = id, x=x, y=y, time = time)
  
  return(DF)
}


sim.lin.joint.be <- function(N, J, beta, D, Sigma, timerange = c(0,5), ...){
  # N: number of individuals
  # J: number of repeated measurements per individual
  # beta: list(beta.y = c(Intercept, time), beta.x = c(Intercept, time))
  # sig: c(sig.y, sig.x)
  # D: covariance matrix of the random effects, 4x4
  # Sigma: covariance matrix of the longitudinal variables, 2x2
  # timerange: maximum range of time of observations
  
  id <- rep(1:N, each = J)
  time <- unlist(lapply(1:N, function(i){
    sort(runif(n = J, min = timerange[1], max = timerange[2]))
  }))
  
  Z <- cbind(Intercept = 1, time)
  
  b <- mvrnorm(N, c(0,0,0,0), D)
  
  mu.x <- Z %*% beta[[2]] + rowSums(Z * b[match(id, unique(id)), 3:4])
  mu.y <- Z %*% beta[[1]] + rowSums(Z * b[match(id, unique(id)), 1:2])
  
  y <- t(sapply(1:length(mu.y), function(i){
    mvrnorm(1, mu = cbind(mu.y, mu.x)[i, ], Sigma = Sigma)
  }))
  
  DF <- data.frame(id = id, x=y[,2], y=y[,1], time = time)
  
  return(DF)
}


sim.data <- list("sim.lin.seq" = sim.lin.seq,
                 "sim.lin.seq.b" = sim.lin.seq.b,
                 "sim.qdr.seq" = sim.qdr.seq,
                 "sim.qdr.seq.b" = sim.qdr.seq.b,
                 "sim.lin.joint.b" = sim.lin.joint.b,
                 "sim.lin.joint.be" = sim.lin.joint.be)
