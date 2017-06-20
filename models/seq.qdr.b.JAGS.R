model {
  for(j in 1:T){
    # Linear mixed effects model for y (and x)
    y[j, 1] ~ dnorm(mu.y[j], tau.y)
    y[j, 2] ~ dnorm(mu.x[j], tau.x)
    mu.y[j] <- inprod(b[subj[j], 1:2], Z[j, ]) + beta[3]*y[j, 2] + beta[4] * pow(y[j, 2], 2)
    mu.x[j] <- inprod(b[subj[j], 3:4], Z[j, ])
  }

  # Priors for the covariance of y (and x)
  tau.y ~ dgamma(0.01, 0.01)
  sig.y <- sqrt(1/tau.y)

  tau.x ~ dgamma(0.01, 0.01)
  sig.x <- sqrt(1/tau.x)

  for(i in 1:N){
    # random effects of y
    b[i, 1:4] ~ dmnorm(mu.b[i, 1:4], inv.D[ , ])
    mu.b[i, 1] <- beta[1]
    mu.b[i, 2] <- beta[2]

    # random effects of x
    mu.b[i, 3] <- alpha[1]
    mu.b[i, 4] <- alpha[2]
  }

  # Priors for the regression coefficients
  for(k in 1:4){
    beta[k] ~ dnorm(0, 0.001)
  }
  for(k in 1:2){
    alpha[k] ~ dnorm(0, 0.001)
  }


  # Priors for random effects part
  #################################
  for(k in 1:4){
    priorR.invD[k,k] ~ dgamma(0.1,0.01)
  }
  inv.D[1:4, 1:4] ~ dwish(priorR.invD[1:4, 1:4], 4)
  D[1:4, 1:4] <- inverse(inv.D[, ])
}
