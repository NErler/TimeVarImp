model{
  # main model
  for (j in 1:T1) {
    bmi[j] ~ dnorm(mu.bmi[j], tau.bmi)
    mu.bmi[j] <- inprod(b.bmi[subj1[j], ], Z.bmi[j, ])
  }

  for (i in 1:N) {
    b.bmi[i, 1:4] ~ dmnorm(mu.b.bmi[i, ], inv.D.bmi[,])
    mu.b.bmi[i, 1] <- inprod(beta[4:11], Xc[i, ]) + beta[12]*BMI_M.c[i] +
      beta[13] * d.c[i,2] + beta[14] * d.c[i,3] + beta[15] * d.c[i,4]

    mu.b.bmi[i, 2] <- beta[1]
    mu.b.bmi[i, 3] <- beta[2]
    mu.b.bmi[i, 4] <- beta[3]
  }


  # model for GWEIGHT
  for (j in 1:T2) {
    GWEIGHT[j] ~ dnorm(mu.GWEIGHT[j], tau.GWEIGHT)
    mu.GWEIGHT[j] <- inprod(b.GWEIGHT[subj2[j], ], Z.GWEIGHT[j, ])
  }

  for (i in 1:N) {
    b.GWEIGHT[i, 1:3] ~ dmnorm(mu.b.GWEIGHT[i, ], inv.D.GWEIGHT[,])
    mu.b.GWEIGHT[i, 1] <- inprod(alpha[3:10], Xc[i, ])
    mu.b.GWEIGHT[i, 2] <- alpha[1] * Xc[i, 1]
    mu.b.GWEIGHT[i, 3] <- alpha[2] * Xc[i, 1]


    d[i, 2] <- (GWEIGHT[pos[i,2]] - GWEIGHT[pos[i,1]])/14
    d[i, 3] <- (GWEIGHT[pos[i,3]] - GWEIGHT[pos[i,2]])/13
    d[i, 4] <- (GWEIGHT[pos[i,4]] - GWEIGHT[pos[i,3]])/(GESTBIR[i] - 27)

    d.c[i, 2] <- (d[i,2] - center.d1)/scale.d1
    d.c[i, 3] <- (d[i,3] - center.d2)/scale.d2
    d.c[i, 4] <- (d[i,4] - center.d3)/scale.d3

    BMI_M[i] <- GWEIGHT[pos[i,1]]/pow(HEIGHT_M[i]/100, 2)
    BMI_M.c[i] <- (BMI_M[i] - center.BMI_M)/scale.BMI_M


    # Imputation of PARITY
    ############################
    Xc[i, 3] ~ dbern(p.PARITY[i])
    logit(p.PARITY[i]) <- inprod(alpha[11:12], Xc[i, 1:2])


    # Imputation of ETHN
    ############################
    Xc[i, 4] ~ dbern(p.ETHN[i])
    logit(p.ETHN[i]) <- inprod(alpha[13:15], Xc[i, 1:3])


    # Imputation of EDUC
    #################################
    cat.data[i, 1] ~ dcat(p.EDUC[i, 1:3])
    eta.EDUC[i] <- inprod(alpha[16:18], Xc[i, 2:4])

    p.EDUC[i, 1] <- max(1e-6, min(0.999999, f.EDUC[i,1]))
    p.EDUC[i, 2] <- max(1e-6, min(0.999999, f.EDUC[i, 2] - f.EDUC[i, 1]))
    p.EDUC[i, 3] <- 1-max(1e-6, min(0.999999, sum(p.EDUC	[i, 1:2])))

    logit(f.EDUC[i, 1])  <- gamma.EDUC[1] + eta.EDUC[i]
    logit(f.EDUC[i, 2])  <- gamma.EDUC[2] + eta.EDUC[i]

    Xc[i, 5] <- ifelse(cat.data[i, 1] == 2, 1, 0)
    Xc[i, 6] <- ifelse(cat.data[i, 1] == 3, 1, 0)


    # Imputation of SMOKE
    #################################
    cat.data[i, 2] ~ dcat(p.SMOKE[i, 1:3])
    eta.SMOKE[i] <- inprod(alpha[19:23], Xc[i, 2:6])

    p.SMOKE[i, 1] <- max(1e-6, min(0.999999, f.SMOKE[i,1]))
    p.SMOKE[i, 2] <- max(1e-6, min(0.999999, f.SMOKE[i, 2] - f.SMOKE[i, 1]))
    p.SMOKE[i, 3] <- 1-max(1e-6, min(0.999999, sum(p.SMOKE	[i, 1:2])))

    logit(f.SMOKE[i, 1])  <- gamma.SMOKE[1] + eta.SMOKE[i]
    logit(f.SMOKE[i, 2])  <- gamma.SMOKE[2] + eta.SMOKE[i]

    Xc[i, 7] <- ifelse(cat.data[i, 2] == 2, 1, 0)
    Xc[i, 8] <- ifelse(cat.data[i, 2] == 3, 1, 0)

  }

  # Priors for longitudinal part
  ################################
  for (k in 1:12) {
    beta[k] ~ dnorm(0, 0.01)
  }

  for (k in 13:15) {
    beta[k] ~ dnorm(0, tau.bmi*(1/tau.beta[k-12] + lambda[2]))
    tau.beta[k-12] ~ dexp(lambda[1]/2)
  }
  lambda[1] ~ dgamma(1, 0.001)
  lambda[2] ~ dgamma(1, 0.001)


  tau.bmi ~ dgamma(0.01, 0.01)
  sigma.bmi <- sqrt(1/tau.bmi)


  # Priors for imputation part
  ################################
  for (k in 1:10) {
    alpha[k] ~ dnorm(0, 0.001)
  }

  # Prior for PARITY
  for (k in 11:12) {
    alpha[k] ~ dnorm(0, 4/9)
  }

  # Prior for ETHN
  for (k in 13:15) {
    alpha[k] ~ dnorm(0, 4/9)
  }

  # Prior for EDUC
  for (k in 16:18) {
    alpha[k] ~ dnorm(0, 4/9)
  }
  delta.EDUC[1] ~ dnorm(0, 1.0e-3)
  gamma.EDUC[1] <- delta.EDUC[1]
  delta.EDUC[2] ~ dnorm(0, 1.0e-3)
  gamma.EDUC[2] <- gamma.EDUC[1] + exp(delta.EDUC[2])

  # Prior for SMOKE
  for (k in 19:23) {
    alpha[k] ~ dnorm(0, 4/9)
  }
  delta.SMOKE[1] ~ dnorm(0, 1.0e-3)
  gamma.SMOKE[1] <- delta.SMOKE[1]
  delta.SMOKE[2] ~ dnorm(0, 1.0e-3)
  gamma.SMOKE[2] <- gamma.SMOKE[1] + exp(delta.SMOKE[2])


  tau.GWEIGHT ~ dgamma(0.01, 0.01)
  sigma.GWEIGHT <- sqrt(1/tau.GWEIGHT)


  # Priors for random effects part
  #################################
  for (k in 1:4) {
    priorR.invD.bmi[k,k] ~ dgamma(1, 0.001)
  }
  inv.D.bmi[1:4, 1:4] ~ dwish(priorR.invD.bmi[, ], 4)
  D.bmi[1:4, 1:4] <- inverse(inv.D.bmi[, ])


  for (k in 1:3) {
    priorR.invD.GWEIGHT[k,k] ~ dgamma(1,0.001)
  }
  inv.D.GWEIGHT[1:3, 1:3] ~ dwish(priorR.invD.GWEIGHT[, ], 3)
  D.GWEIGHT[1:3, 1:3] <- inverse(inv.D.GWEIGHT[, ])
}
