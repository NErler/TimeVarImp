#********************************#
# Multivariate normal imputation #
#********************************#


run.joint.b <- function(df.mis, n.iter, n.chain, ipc, n.adapt = 100, 
                        n.burnin = 0, n.thin = 1, showplot = F){
  # data list
  l <- list()
  l$y <- cbind(df.mis$y, df.mis$x)
  l$Z <- cbind(1,df.mis$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(df.mis$id))
  l$subj <- df.mis$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- rep(NA, 4)
  
  
  # parameters
  ys <- paste0("y[", which(is.na(l$y[,2])), ",2]")
  
  params <- c("beta", "D", "sig.y", "sig.x", ys)
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$beta <- matrix(rnorm(4), nrow = 2, ncol = 2)
    ll$tau.x <- rgamma(1,1,1)
    ll$tau.y <- rgamma(1,1,1)
    return(ll)
  }
  
  # run JAGS
  adapt.joint.b <- jags.model(file = "models/joint.b.JAGS.R", l,
                              n.chains = n.chain, n.adapt = n.adapt,
                              inits = get.inits)
  
  
  gr.y <- gr.other <- 10
  burnin <- 0
  while((gr.y > 1.1 | gr.other > 1.25) & burnin < 5000){
    mcmc.joint.b <- coda.samples(adapt.joint.b, params, 
                                 n.iter = n.iter, n.thin = 1)
    gr.y <- max(gelman.diag(
      mcmc.joint.b[, colnames(mcmc.joint.b[[1]]) %in% sample(ys, 50)],
      autoburnin = F, multivariate = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.joint.b[, !colnames(mcmc.joint.b[[1]]) %in% ys], autoburnin = F,
      multivariate = F)[[1]][,2])
    burnin <- start(mcmc.joint.b) - 1
    print(gr.y)
    print(gr.other)
  }
  
  # get imputed values
  ximp.joint.b <- mcmc.joint.b[c(1:ipc) * n.iter/ipc, ys]
  ximp.joint.b <- do.call(rbind, ximp.joint.b)
  
  DF.list <- list()
  for(i in 1:nrow(ximp.joint.b)){
    DF.list[[i]] <- df.mis
    DF.list[[i]]$x[is.na(df.mis$x)] <- ximp.joint.b[i, ]
  }
  MI.list.joint.b <- datalist2mids(DF.list)
  
  if(showplot==T){
    windows()
    par(mfrow = c(5,6), mar = c(3,3,2,1))
    traceplot(mcmc.joint.b[, 1:min(30, ncol(mcmc.joint.b[[1]]))])
  }
  
  
  return(list(burnin = start(mcmc.joint.b) - 1, niter = niter(mcmc.joint.b),
              imp = MI.list.joint.b))
}





run.joint.be <- function(df.mis, n.iter, n.chain, ipc, n.adapt = 100, 
                         n.burnin = 0, n.thin = 1, showplot = F){
  # data list
  l <- list()
  l$y <- cbind(df.mis$y, df.mis$x)
  l$Z <- cbind(1,df.mis$time)
  l$T <- nrow(l$y)
  l$N <- length(unique(df.mis$id))
  l$subj <- df.mis$id
  l$priorR.invD <- matrix(nrow = 4, ncol = 4, data = 0)
  diag(l$priorR.invD) <- rep(NA, 4)
  l$priorR.invS <- matrix(nrow = 2, ncol = 2, data = 0)
  diag(l$priorR.invS) <- rep(NA, 2)
  
  
  # parameters
  ys <- paste0("y[", which(is.na(l$y[,2])), ",2]")
  
  params <- c("beta", "D", "S", ys)
  
  # initial values
  get.inits <- function(){
    ll <- list()
    ll$beta <- matrix(rnorm(4), nrow = 2, ncol = 2)
    return(ll)
  }
  
  # run JAGS
  adapt.joint.be <- jags.model(file = "models/joint.be.JAGS.R", l, 
                               n.chains = n.chain, n.adapt = n.adapt,
                               inits = get.inits)
  
  
  gr.y <- gr.other <- 10
  burnin <- 0
  while((gr.y > 1.1 | gr.other > 1.25) & burnin < 5000){
    mcmc.joint.be <- coda.samples(adapt.joint.be, params, n.iter = n.iter)
    gr.y <- max(gelman.diag(
      mcmc.joint.be[, colnames(mcmc.joint.be[[1]]) %in% sample(ys, 50)],
      autoburnin = F, multivariate = F)[[1]][,2])
    gr.other <- max(gelman.diag(
      mcmc.joint.be[, !colnames(mcmc.joint.be[[1]]) %in% ys], autoburnin = F,
      multivariate = F)[[1]][,2])
    burnin <- start(mcmc.joint.be) - 1
    print(gr.y)
    print(gr.other)
  }
  
  # get imputed values
  ximp.joint.be <- mcmc.joint.be[c(1:ipc) * n.iter/ipc, ys]
  ximp.joint.be <- do.call(rbind, ximp.joint.be)
  
  DF.list <- list()
  for(i in 1:nrow(ximp.joint.be)){
    DF.list[[i]] <- df.mis
    DF.list[[i]]$x[is.na(df.mis$x)] <- ximp.joint.be[i, ]
  }
  MI.list.joint.be <- datalist2mids(DF.list)
  
  if(showplot==T){
    windows()
    par(mfrow = c(6,6), mar = c(3,3,2,1))
    traceplot(mcmc.joint.be[, 1:min(36, ncol(mcmc.joint.be[[1]]))])
  }
  
  
  return(list(burnin = start(mcmc.joint.be) - 1, niter = niter(mcmc.joint.be),
              imp = MI.list.joint.be))
}
