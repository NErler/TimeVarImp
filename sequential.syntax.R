#*******************#
# Sequential models #
#*******************#


# complete data ----------------------------------------------------------------
run.seq.lin.all <- function(DF, n.iter, n.chain, n.adapt = 100, n.burnin = 0,
                            n.thin = 1, showplot = F){
  # data list
  l <- list()
  l$y <- cbind(DF$y, DF$x)
  l$Z <- cbind(1,DF$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(DF$id))
  l$subj <- DF$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- NA
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$tau.y <- rgamma(1,1,1)
    ll$tau.x <- rgamma(1,1,1)
    ll$beta <- rnorm(3)
    return(ll)
  }
  
  # parameters
  params <- c("beta", "D.y", "sig.y")
  
  # run JAGS
  adapt.seq.all <- jags.model(file = "models/seq.lin.JAGS.R", l,
                              n.chains = n.chain, n.adapt = n.adapt,
                              inits = get.inits)
  
  gr.beta <- gr.other <- 10
  burnin <- 0
  while((gr.beta > 1.1 | gr.other > 1.25) & burnin < 10000){
    mcmc.seq.all <- coda.samples(adapt.seq.all, params, n.iter = n.iter)
    
    gr.beta <- max(gelman.diag(
      mcmc.seq.all[, grep("beta", colnames(mcmc.seq.all[[1]]))],
      autoburnin = F)[[1]][,2])
    
    gr.other <- max(gelman.diag(
      mcmc.seq.all[, !colnames(mcmc.seq.all[[1]]) %in%
                     grep("beta", colnames(mcmc.seq.all[[1]]), value = T)],
      autoburnin = F, multivariate = F)[[1]][,2])
    burnin <- start(mcmc.seq.all) - 1
    print(gr.beta)
    print(gr.other)
  }
  
  
  mcmc.seq.all.combi <- mcmc.seq.all[, grep("beta", colnames(mcmc.seq.all[[1]]))]
  MCErr <- mcse.mat(do.call(rbind,mcmc.seq.all.combi)
  )[,2]/summary(mcmc.seq.all.combi)$stat[,2]
  print(MCErr)
  while(max(MCErr) > 0.05){
    mcmc.seq.all.add <- coda.samples(adapt.seq.all, "beta", n.iter = n.iter)
    for(i in 1:length(mcmc.seq.all)){
      mcmc.seq.all.combi[[i]] <- as.mcmc(rbind(mcmc.seq.all.combi[[i]],
                                               mcmc.seq.all.add[[i]]))
    }
    MCErr <- mcse.mat(do.call(rbind, mcmc.seq.all.combi)
    )[,2]/summary(mcmc.seq.all.combi)$stat[,2]
    print(MCErr)
  }
  
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.seq.all.combi[, 1:min(30, ncol(mcmc.seq.all.combi[[1]]))])
  }
  
  return(list(burnin = burnin, niter = niter(mcmc.seq.all.combi), params = summary(mcmc.seq.all.combi)))
}


run.seq.b.lin.all <- function(DF, n.iter, n.chain, n.adapt = 100, n.burnin = 0,
                              n.thin = 1, showplot = F){
  # data list
  l <- list()
  l$y <- cbind(DF$y, DF$x)
  l$Z <- cbind(1,DF$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(DF$id))
  l$subj <- DF$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- NA
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$tau.y <- rgamma(1,1,1)
    ll$tau.x <- rgamma(1,1,1)
    ll$beta <- rnorm(3)
    return(ll)
  }
  
  # parameters
  params <- c("beta", "D", "sig.y", "sig.x", "alpha")
  
  # run JAGS
  adapt.seq.b.all <- jags.model(file = "models/seq.lin.b.JAGS.R", l,
                                n.chains = n.chain, n.adapt = n.adapt,
                                inits = get.inits)
  
  gr.beta <- gr.other <- 10
  burnin <- 0
  while((gr.beta > 1.1 | gr.other > 1.25) & burnin < 10000){
    mcmc.seq.b.all <- coda.samples(adapt.seq.b.all, params, n.iter = n.iter)
    gr.beta <- max(gelman.diag(
      mcmc.seq.b.all[, grep("beta", colnames(mcmc.seq.b.all[[1]]))],
      autoburnin = F, multivariate = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.seq.b.all[, !colnames(mcmc.seq.b.all[[1]]) %in% 
                       grep("beta", colnames(mcmc.seq.b.all[[1]]), value = T)],
      autoburnin = F, multivariate = F)[[1]][,2])
    burnin <- start(mcmc.seq.b.all) - 1
    print(gr.beta)
    print(gr.other)
  }
  
  
  mcmc.seq.b.all.combi <- mcmc.seq.b.all[, grep("beta", 
                                                colnames(mcmc.seq.b.all[[1]]))]
  MCErr <- mcse.mat(do.call(rbind, mcmc.seq.b.all.combi)
  )[,2]/summary(mcmc.seq.b.all.combi)$stat[,2]
  print(MCErr)
  while(max(MCErr) > 0.05){
    mcmc.seq.b.all.add <- coda.samples(adapt.seq.b.all, "beta", n.iter = n.iter)
    for(i in 1:length(mcmc.seq.b.all)){
      mcmc.seq.b.all.combi[[i]] <- as.mcmc(rbind(mcmc.seq.b.all.combi[[i]], 
                                                 mcmc.seq.b.all.add[[i]]))
    }
    MCErr <- mcse.mat(do.call(rbind, mcmc.seq.b.all.combi)
    )[,2]/summary(mcmc.seq.b.all.combi)$stat[,2]
    print(MCErr)
  }
  
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.seq.b.all.combi[, 1:min(30, ncol(mcmc.seq.b.all.combi[[1]]))])
  }
  
  return(list(burnin = burnin, niter = niter(mcmc.seq.b.all.combi), 
              params = summary(mcmc.seq.b.all.combi)))
}




run.seq.qdr.all <- function(DF, n.iter, n.chain, n.adapt = 100, n.burnin = 0, 
                            n.thin = 1, showplot = F){
  # data list
  l <- list()
  l$y <- cbind(DF$y, DF$x)
  l$Z <- cbind(1,DF$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(DF$id))
  l$subj <- DF$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- NA
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$tau.y <- rgamma(1,1,1)
    ll$tau.x <- rgamma(1,1,1)
    ll$beta <- rnorm(4)
    return(ll)
  }
  
  # parameters
  params <- c("beta", "D.y", "sig.y")
  
  # run JAGS
  adapt.seq.qdr.all <- jags.model(file = "models/seq.qdr.JAGS.R", l,
                                  n.chains = n.chain, n.adapt = n.adapt,
                                  inits = get.inits)
  
  gr.beta <- gr.other <- 10
  burnin <- 0
  while((gr.beta > 1.1 | gr.other > 1.25) & burnin < 10000){
    mcmc.seq.qdr.all <- coda.samples(adapt.seq.qdr.all, params, n.iter = n.iter)
    gr.beta <- max(gelman.diag(
      mcmc.seq.qdr.all[, grep("beta", colnames(mcmc.seq.qdr.all[[1]]))],
      autoburnin = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.seq.qdr.all[, !colnames(mcmc.seq.qdr.all[[1]]) %in% 
                         grep("beta", colnames(mcmc.seq.qdr.all[[1]]), value = T)],
      autoburnin = F, multivariate = F)[[1]][,2])
    burnin <- start(mcmc.seq.qdr.all) - 1
    print(gr.beta)
    print(gr.other)
  }
  
  
  mcmc.seq.qdr.all.combi <- mcmc.seq.qdr.all[, grep("beta",
                                                    colnames(mcmc.seq.qdr.all[[1]]))]
  MCErr <- mcse.mat(do.call(rbind, mcmc.seq.qdr.all.combi)
  )[,2]/summary(mcmc.seq.qdr.all.combi)$stat[,2]
  print(MCErr)
  while(max(MCErr) > 0.05){
    mcmc.seq.qdr.all.add <- coda.samples(adapt.seq.qdr.all, "beta", n.iter = n.iter)
    for(i in 1:length(mcmc.seq.qdr.all)){
      mcmc.seq.qdr.all.combi[[i]] <- as.mcmc(rbind(mcmc.seq.qdr.all.combi[[i]],
                                                   mcmc.seq.qdr.all.add[[i]]))
    }
    MCErr <- mcse.mat(do.call(rbind, mcmc.seq.qdr.all.combi)
    )[,2]/summary(mcmc.seq.qdr.all.combi)$stat[,2]
    print(MCErr)
  }
  
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.seq.qdr.all.combi[, 1:min(30, ncol(mcmc.seq.qdr.all.combi[[1]]))])
  }
  
  return(list(burnin = burnin, niter = niter(mcmc.seq.qdr.all.combi),
              params = summary(mcmc.seq.qdr.all.combi)))
}



run.seq.b.qdr.all <- function(DF, n.iter, n.chain, n.adapt = 100, n.burnin = 0,
                              n.thin = 1, showplot = F){
  # data list
  l <- list()
  l$y <- cbind(DF$y, DF$x)
  l$Z <- cbind(1,DF$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(DF$id))
  l$subj <- DF$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- NA
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$tau.y <- rgamma(1,1,1)
    ll$tau.x <- rgamma(1,1,1)
    ll$beta <- rnorm(4)
    return(ll)
  }
  
  # parameters
  params <- c("beta", "D", "sig.y", "sig.x", "alpha")
  
  # run JAGS
  adapt.seq.b.qdr.all <- jags.model(file = "models/seq.qdr.b.JAGS.R", l, 
                                    n.chains = n.chain, n.adapt = n.adapt,
                                    inits = get.inits)
  
  gr.beta <- gr.other <- 10
  burnin <- 0
  while((gr.beta > 1.1 | gr.other > 1.25) & burnin < 10000){
    mcmc.seq.b.qdr.all <- coda.samples(adapt.seq.b.qdr.all, params, n.iter = n.iter)
    gr.beta <- max(gelman.diag(
      mcmc.seq.b.qdr.all[, grep("beta", colnames(mcmc.seq.b.qdr.all[[1]]))],
      autoburnin = F, multivariate = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.seq.b.qdr.all[, !colnames(mcmc.seq.b.qdr.all[[1]]) %in% 
                           grep("beta", colnames(mcmc.seq.b.qdr.all[[1]]), value = T)],
      autoburnin = F, multivariate = F)[[1]][,2])
    burnin <- start(mcmc.seq.b.qdr.all) - 1
    print(gr.beta)
    print(gr.other)
  }
  
  
  mcmc.seq.b.qdr.all.combi <- mcmc.seq.b.qdr.all[, grep("beta", colnames(mcmc.seq.b.qdr.all[[1]]))]
  MCErr <- mcse.mat(do.call(rbind, mcmc.seq.b.qdr.all.combi)
  )[,2]/summary(mcmc.seq.b.qdr.all.combi)$stat[,2]
  print(MCErr)
  while(max(MCErr) > 0.05){
    mcmc.seq.b.qdr.all.add <- coda.samples(adapt.seq.b.qdr.all, "beta",
                                           n.iter = n.iter)
    for(i in 1:length(mcmc.seq.b.qdr.all)){
      mcmc.seq.b.qdr.all.combi[[i]] <- as.mcmc(rbind(mcmc.seq.b.qdr.all.combi[[i]],
                                                     mcmc.seq.b.qdr.all.add[[i]]))
    }
    MCErr <- mcse.mat(do.call(rbind, mcmc.seq.b.qdr.all.combi)
    )[,2]/summary(mcmc.seq.b.qdr.all.combi)$stat[,2]
    print(MCErr)
  }
  
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.seq.b.qdr.all.combi[, 1:min(30, ncol(mcmc.seq.b.qdr.all.combi[[1]]))])
  }
  
  return(list(burnin = burnin, niter = niter(mcmc.seq.b.qdr.all.combi),
              params = summary(mcmc.seq.b.qdr.all.combi)))
}



# imputation  ##################################################################
run.seq.lin <- function(DF, n.iter, n.chain, n.adapt = 100, n.burnin = 0, 
                        n.thin = 1, showplot = F, n.burn = 0){
  # data list
  l <- list()
  l$y <- cbind(DF$y, DF$x)
  l$Z <- cbind(1,DF$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(DF$id))
  l$subj <- DF$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- NA
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$tau.y <- rgamma(1,1,1)
    ll$tau.x <- rgamma(1,1,1)
    ll$beta <- rnorm(3)
    return(ll)
  }
  
  # parameters
  params <- c("beta", "D.y", "sig.y", "D.x", "sig.x", "alpha")
  
  # run JAGS
  adapt.seq.lin <- jags.model(file = "models/seq.lin.JAGS.R", l,
                              n.chains = n.chain, n.adapt = n.adapt,
                              inits = get.inits)
  
  if(n.burn>0){
    mcmc.seq.lin <- coda.samples(adapt.seq.lin, params, n.iter = n.burn)
  }
  gr.beta <- gr.other <- 10
  burnin <- n.burn
  while((gr.beta > 1.1 | gr.other > 1.25) & burnin < 10000){
    mcmc.seq.lin <- coda.samples(adapt.seq.lin, params, n.iter = n.iter)
    gr.beta <- max(gelman.diag(
      mcmc.seq.lin[, grep("beta", colnames(mcmc.seq.lin[[1]]))], 
      autoburnin = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.seq.lin[, !colnames(mcmc.seq.lin[[1]]) %in%
                     grep("beta", colnames(mcmc.seq.lin[[1]]), value = T)],
      autoburnin = F, multivariate = F)[[1]][,2])
    burnin <- start(mcmc.seq.lin) - 1
    print(gr.beta)
    print(gr.other)
  }
  
  mcmc.seq.lin.combi <- mcmc.seq.lin[, grep("beta", colnames(mcmc.seq.lin[[1]]))]
  MCErr <- mcse.mat(do.call(rbind, mcmc.seq.lin.combi)
  )[,2]/summary(mcmc.seq.lin.combi)$stat[,2]
  while(max(MCErr) > 0.05){
    mcmc.seq.lin.add <- coda.samples(adapt.seq.lin, "beta", n.iter = n.iter)
    for(i in 1:length(mcmc.seq.lin)){
      mcmc.seq.lin.combi[[i]] <- as.mcmc(rbind(mcmc.seq.lin.combi[[i]],
                                               mcmc.seq.lin.add[[i]]))
    }
    MCErr <- mcse.mat(do.call(rbind, mcmc.seq.lin.combi)
    )[,2]/summary(mcmc.seq.lin.combi)$stat[,2]
    print(MCErr)
  }
  
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.seq.lin.combi[, 1:min(30, ncol(mcmc.seq.lin.combi[[1]]))])
  }
  
  return(list(burnin = burnin, niter = niter(mcmc.seq.lin.combi),
              params = summary(mcmc.seq.lin.combi)))
}


run.seq.qdr <- function(DF, n.iter, n.chain, n.adapt = 100, n.burnin = 0,
                        n.thin = 1, showplot = F){
  # data list
  l <- list()
  l$y <- cbind(DF$y, DF$x)
  l$Z <- cbind(1,DF$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(DF$id))
  l$subj <- DF$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- NA
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$tau.y <- rgamma(1,1,1)
    ll$tau.x <- rgamma(1,1,1)
    ll$beta <- rnorm(4)
    return(ll)
  }
  
  # parameters
  params <- c("beta", "D.y", "sig.y", "D.x", "sig.x", "alpha")
  
  # run JAGS
  adapt.seq.qdr <- jags.model(file = "models/seq.qdr.JAGS.R", l,
                              n.chains = n.chain, n.adapt = n.adapt,
                              inits = get.inits)
  
  gr.beta <- gr.other <- 10
  burnin <- 0
  while((gr.beta > 1.1 | gr.other > 1.25) & burnin < 10000){
    mcmc.seq.qdr <- coda.samples(adapt.seq.qdr, params, n.iter = n.iter)
    gr.beta <- max(gelman.diag(
      mcmc.seq.qdr[, grep("beta", colnames(mcmc.seq.qdr[[1]]))], 
      autoburnin = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.seq.qdr[, !colnames(mcmc.seq.qdr[[1]]) %in% 
                     grep("beta", colnames(mcmc.seq.qdr[[1]]), value = T)],
      autoburnin = F, multivariate = F)[[1]][,2])
    burnin <- start(mcmc.seq.qdr) - 1
    print(gr.beta)
    print(gr.other)
  }
  
  
  mcmc.seq.qdr.combi <- mcmc.seq.qdr[, grep("beta", colnames(mcmc.seq.qdr[[1]]))]
  MCErr <- mcse.mat(do.call(rbind, mcmc.seq.qdr.combi)
  )[,2]/summary(mcmc.seq.qdr.combi)$stat[,2]
  while(max(MCErr) > 0.05){
    mcmc.seq.qdr.add <- coda.samples(adapt.seq.qdr, "beta", n.iter = n.iter)
    for(i in 1:length(mcmc.seq.qdr)){
      mcmc.seq.qdr.combi[[i]] <- as.mcmc(rbind(mcmc.seq.qdr.combi[[i]],
                                               mcmc.seq.qdr.add[[i]]))
    }
    MCErr <- mcse.mat(do.call(rbind, mcmc.seq.qdr.combi)
    )[,2]/summary(mcmc.seq.qdr.combi)$stat[,2]
    print(MCErr)
  }
  
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.seq.qdr.combi[, 1:min(30, ncol(mcmc.seq.qdr.combi[[1]]))])
  }
  
  return(list(burnin = burnin, niter = niter(mcmc.seq.qdr.combi), 
              params = summary(mcmc.seq.qdr.combi)))
}



run.seq.b.lin <- function(DF, n.iter, n.chain, n.adapt = 100, n.burnin = 0,
                          n.thin = 1, showplot = F, fixit = F,
                          return.chain = F){
  # data list
  l <- list()
  l$y <- cbind(DF$y, DF$x)
  l$Z <- cbind(1,DF$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(DF$id))
  l$subj <- DF$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- NA
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$tau.y <- rgamma(1,1,1)
    ll$tau.x <- rgamma(1,1,1)
    ll$beta <- rnorm(3)
    return(ll)
  }
  
  # parameters
  params <- c("beta", "D", "sig.y", "sig.x", "alpha")
  
  # run JAGS
  adapt.seq.b.lin <- jags.model(file = "models/seq.lin.b.JAGS.R", l, 
                                n.chains = n.chain, n.adapt = n.adapt,
                                inits = get.inits)
  
  gr.beta <- gr.other <- 10
  burnin <- 0
  while((gr.beta > 1.1 | gr.other > 1.25) & burnin < 10000){
    mcmc.seq.b.lin <- coda.samples(adapt.seq.b.lin, params, n.iter = n.iter)
    gr.beta <- max(gelman.diag(
      mcmc.seq.b.lin[, grep("beta", colnames(mcmc.seq.b.lin[[1]]))],
      autoburnin = F, multivariate = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.seq.b.lin[, !colnames(mcmc.seq.b.lin[[1]]) %in% 
                       grep("beta", colnames(mcmc.seq.b.lin[[1]]), value = T)],
      autoburnin = F, multivariate = F)[[1]][,2])
    burnin <- start(mcmc.seq.b.lin) - 1
    print(gr.beta)
    print(gr.other)
  }
  
  
  mcmc.seq.b.lin.combi <- mcmc.seq.b.lin[, grep("beta", colnames(mcmc.seq.b.lin[[1]]))]
  MCErr <- mcse.mat(do.call(rbind, mcmc.seq.b.lin.combi)
  )[,2]/summary(mcmc.seq.b.lin.combi)$stat[,2]
  if(fixit == F){
    while(max(MCErr) > 0.05){
      mcmc.seq.b.lin.add <- coda.samples(adapt.seq.b.lin, "beta", n.iter = n.iter)
      for(i in 1:length(mcmc.seq.b.lin)){
        mcmc.seq.b.lin.combi[[i]] <- as.mcmc(rbind(mcmc.seq.b.lin.combi[[i]],
                                                   mcmc.seq.b.lin.add[[i]]))
      }
      MCErr <- mcse.mat(do.call(rbind, mcmc.seq.b.lin.combi)
      )[,2]/summary(mcmc.seq.b.lin.combi)$stat[,2]
      print(MCErr)
    }
  }
  print(MCErr)
  
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.seq.b.lin.combi[, 1:min(30, ncol(mcmc.seq.b.lin.combi[[1]]))])
  }
  
  chain <- if(return.chain == T){
    mcmc.seq.b.lin.combi
  }else{
    NULL
  }
  
  return(list(burnin = burnin, niter = niter(mcmc.seq.b.lin.combi), 
              params = summary(mcmc.seq.b.lin.combi), chain = chain))
}




run.seq.b.qdr <- function(DF, n.iter, n.chain, n.adapt = 100, n.burnin = 0,
                          n.thin = 1, showplot = F, fixit = F,
                          return.chain = F){
  # data list
  l <- list()
  l$y <- cbind(DF$y, DF$x)
  l$Z <- cbind(1,DF$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(DF$id))
  l$subj <- DF$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- NA
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$tau.y <- rgamma(1,1,1)
    ll$tau.x <- rgamma(1,1,1)
    ll$beta <- rnorm(4)
    return(ll)
  }
  
  # parameters
  params <- c("beta", "D", "sig.y", "sig.x", "alpha")
  
  # run JAGS
  adapt.seq.b.qdr <- jags.model(file = "models/seq.qdr.b.JAGS.R", l, 
                                n.chains = n.chain, n.adapt = n.adapt,
                                inits = get.inits)
  
  gr.beta <- gr.other <- 10
  burnin <- 0
  while((gr.beta > 1.1 | gr.other > 1.25) & burnin < 10000){
    mcmc.seq.b.qdr <- coda.samples(adapt.seq.b.qdr, params, n.iter = n.iter)
    gr.beta <- max(gelman.diag(
      mcmc.seq.b.qdr[, grep("beta", colnames(mcmc.seq.b.qdr[[1]]))],
      autoburnin = F, multivariate = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.seq.b.qdr[, !colnames(mcmc.seq.b.qdr[[1]]) %in% 
                       grep("beta", colnames(mcmc.seq.b.qdr[[1]]), value = T)],
      autoburnin = F, multivariate = F)[[1]][,2])
    burnin <- start(mcmc.seq.b.qdr) - 1
    print(gr.beta)
    print(gr.other)
  }
  
  mcmc.seq.b.qdr.combi <- mcmc.seq.b.qdr[, grep("beta", colnames(mcmc.seq.b.qdr[[1]]))]
  MCErr <- mcse.mat(do.call(rbind, mcmc.seq.b.qdr.combi)
  )[,2]/summary(mcmc.seq.b.qdr.combi)$stat[,2]
  if(fixit == F){
    while(max(MCErr) > 0.05){
      mcmc.seq.b.qdr.add <- coda.samples(adapt.seq.b.qdr, "beta", n.iter = n.iter)
      for(i in 1:length(mcmc.seq.b.qdr)){
        mcmc.seq.b.qdr.combi[[i]] <- as.mcmc(rbind(mcmc.seq.b.qdr.combi[[i]],
                                                   mcmc.seq.b.qdr.add[[i]]))
      }
      MCErr <- mcse.mat(do.call(rbind, mcmc.seq.b.qdr.combi)
      )[,2]/summary(mcmc.seq.b.qdr.combi)$stat[,2]
      print(MCErr)
    }
  }
  print(MCErr)
  
  chain <- if(return.chain == T){
    mcmc.seq.b.qdr.combi
  }else{
    NULL
  }
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.seq.b.qdr.combi[, 1:min(30, ncol(mcmc.seq.b.qdr.combi[[1]]))])
  }
  
  return(list(burnin = burnin, niter = niter(mcmc.seq.b.qdr.combi), 
              params = summary(mcmc.seq.b.qdr.combi), chain = chain))
}
