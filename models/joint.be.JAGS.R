model {
  for(j in 1:T){
    # Linear mixed effects model for y (and x)
    y[j, 1] ~ dnorm(mu.y[j], tau.y)
    y[j, 2] ~ dnorm(mu.x[j], tau.x)
    mu.y[j] <- inprod(b[subj[j], 1:2], Z[j, ]) + eps[j,1]
    mu.x[j] <- inprod(b[subj[j], 3:4], Z[j, ]) + eps[j,2]

    eps[j, 1:2] ~ dmnorm(mu.eps[], inv.S[,])
  }
  mu.eps[1] <- 0
  mu.eps[2] <- 0

  # Priors for the covariance of y (and x)
  for(k in 1:2){
    priorR.invS[k,k] ~ dgamma(0.1,0.01)
  }
  inv.S[1:2, 1:2] ~ dwish(priorR.invS[1:2, 1:2], 2)
  S[1:2, 1:2] <- inverse(inv.S[, ])

  tau.y  <- 100000
  tau.x <- 100000

  sig.y <- 1/sqrt(tau.y)
  sig.x <- 1/sqrt(tau.x)

  for(i in 1:N){
    # random effects of y
    b[i, 1:4] ~ dmnorm(mu.b[i, 1:4], inv.D[ , ])
    mu.b[i, 1] <- beta[1,1]
    mu.b[i, 2] <- beta[2,1]

    # random effects of x
    mu.b[i, 3] <- beta[1,2]
    mu.b[i, 4] <- beta[2,2]
  }

  # Priors for random effects covariance
  for(k in 1:4){
    priorR.invD[k,k] ~ dgamma(0.1,0.01)
  }
  inv.D[1:4, 1:4] ~ dwish(priorR.invD[1:4, 1:4], 4)
  D[1:4, 1:4] <- inverse(inv.D[, ])

  # Priors for regression coefficients
  for(k in 1:2){
    beta[k,1] ~ dnorm(0, 0.001)
    beta[k,2] ~ dnorm(0, 0.001)
  }
}
