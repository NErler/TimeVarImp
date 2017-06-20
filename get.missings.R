######################################
# functions to create missing values #
######################################

get.MAR.1 <- function(DF, misvar, depvar, const = -3.5, perc, plot = F, ...){
  # DF: complete dataset
  # misvar: name of variable that will be incomplete
  # depvar: name of the variable on which the missingness is depending
  # perc: (approximate) proportion of missing values
  
  df.mis <- DF
  ismis <- rbinom(n = nrow(df.mis), size = 1, 
                  prob = plogis(df.mis[, depvar] + const[1]))
  if(plot == T){
    par(mfrow=c(1,2), mar = c(3,3,2,1), mgp = c(2, 0.6, 0))
    plot(df.mis[, depvar], plogis(df.mis[, depvar] + const[1]), 
         xlab = "y", ylab = "prob")
    hist(plogis(df.mis[, depvar] + const[1]))
  }
  df.mis[ismis == 1, misvar] <- NA
  return(df.mis)
}


get.MAR.2 <- function(DF, misvar, depvar, const = c(-4, 2), N, J, cutoff = 0.2, 
                      perc, plot = F,...){
  # DF: complete dataset, must be balanced (same number of observations per
  #     individual) and sorted by individuals
  # misvar: name of variable that will be incomplete
  # depvar: name of the variable on which the missingness is depending
  # N: number of individuals
  # J: number of repeated values per individual
  # cutoff: value of misvar for which the value of the next observation is more
  #         likely to be missing
  
  df.mis <- DF
  mismat <- matrix(model.matrix(~ get(misvar), df.mis, na.action = NULL)[, -1],
                   nrow = J, ncol = N)
  probmat <- mismat * 0
  for(i in 1:ncol(mismat)){
    for(j in J:2){
      if(mismat[j-1, i] < cutoff){probmat[j, i] <- 1}
    }
  }
  ismis <- rbinom(n = nrow(df.mis), size = 1, 
                  prob = plogis(const[1] + df.mis[, depvar] + const[2] * probmat))
  if(plot == T){
    par(mfrow = c(2,2), mar = c(3,3,2,1), mgp = c(2, 0.6, 0))
    plot(df.mis[, depvar], 
         plogis(const[1] + df.mis[, depvar] + const[2] * probmat),
         xlab = "y", ylab = "prob")
    plot(df.mis[, misvar], 
         plogis(const[1] + df.mis[, depvar] + const[2] * probmat),
         xlab = "x", ylab = "prob")
    hist(plogis(const[1] + df.mis[, depvar] + const[2] * probmat))}
  df.mis[ismis == 1, misvar] <- NA
  return(df.mis)
}


get.mis <- list(MAR.1 = get.MAR.1,
                MAR.2 = get.MAR.2
)
