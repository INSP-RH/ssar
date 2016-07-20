#Tests for checking that time = 0 and time = tmax are included int he
#results. This tests also check that the number of saved iterations are adecquate

context("Checking simulation results")

#Function for getting mean
getMean <- function(simulation, times){
  
  #INPUT:
  #simulation .- Data frame resulting from ssa.R
  #times      .- Vector of times to interpolate
  #
  #OUTPUT: 
  # Mean of the variables "Var" in simulation approximated at time t
  
  #Get number of simulations
  .nsim  <- max(simulation$Simulation)
  
  #Empty list for returning
  .rlist <- list()
  
  #Loop through all simulations
  for (.col in 1:( ncol(simulation) - 3)){
    
    #Matrix for including interpolation points of simulation
    .allsims    <- matrix(NA, ncol = .nsim, nrow = length(times))
    
    #Get data at interpolation points
    for (.i in 1:.nsim){
      
      #Get subse of data to use for aprox
      .subsimulation <- subset(simulation, simulation$Simulation == .i)
      
      #Use aprox
      .allsims[,.i]   <- approx(.subsimulation$Time, .subsimulation[,paste0("Var",.col)], times, method = "constant")$y

    }
    
    #Create list name
    .rlist[[.col]] <- .allsims
    
  }
  
  return(.rlist)
  
}


#Testing
test_that("Simulation results",{
  
  #Evaluating if mean of exponential increase matches theoretical mean
  #by considering that confidence intervals of confidence conf % most
  #include the real value conf% of times. 
  expect_true({
    
    #Seed
    set.seed(6248)
    
    #Simulation values
    params     <- c(b=0.1)
    X          <- matrix(c(N=1.2), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] * X[,1]) }
    v          <- matrix( c(+1), ncol = 1)
    tmin       <- 0
    tmax       <- 10
    nsim       <- 100
    conf       <- 95 #Confidence interval
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, plot.sim = F)
    
    #Matrix for saving the results at specific times
    times      <- seq(tmin, tmax, length.out = 100)
    
    #Function for theoretical mean
    tmean      <- X[,1]*exp(params[1]*times)
      
    #Estimate mean and sd
    summarymat <- getMean(simulation, times)[[1]]
    simmean    <- apply(summarymat, 1, mean)       #Estimate mean
    simsd      <- apply(summarymat, 1, sd)         #Estimate sd 
    
    #Do simean - tmean
    res <- abs(simmean - tmean)
    
    #Check the confidence interval matches conf times
    Z <- qnorm(1-(1-conf/100)/2)
    rejected <- length(which(res > Z*simsd/sqrt(nsim)))
    
    #Check that we are rejecting less than conf%
    (rejected/length(res)) < (100 - conf)/100
    
  })
  
  #Evaluating if mean of exponential decrease matches theoretical mean
  expect_true({
    
    #Seed
    set.seed(6248)
    
    #Simulation values
    params     <- c(b=10)
    X          <- matrix(c(N=10), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] * X[,1]) }
    v          <- matrix( c(-1), ncol = 1)
    tmin       <- 0
    tmax       <- 100
    nsim       <- 100
    conf       <- 95 #Confidence interval
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, plot.sim = T)
    
    #Matrix for saving the results at specific times
    times      <- seq(tmin, tmax, length.out = 1000)
    
    #Function for theoretical mean
    tmean      <- X[,1]*exp(-params[1]*times)
    
    #Estimate mean and sd
    summarymat <- getMean(simulation, times)[[1]]
    simmean    <- apply(summarymat, 1, mean)       #Estimate mean
    simsd      <- apply(summarymat, 1, sd)         #Estimate sd 
    
    #Plot
    lines(times, tmean, col ="black", lwd = 2)
    lines(times, simmean, col = "green", lwd = 2)
    
    #Do simean - tmean
    res <- abs(simmean - tmean)
    
    #Check the confidence interval matches conf times
    Z <- qnorm(1-(1-conf/100)/2)
    rejected <- length(which( abs(res - Z*simsd/sqrt(nsim)) > 1.e-2 ))
    
    #Check that we are rejecting less than conf%
    (rejected/length(res)) < (100 - conf)/100
    
  })
  
  #Evaluating if mean of exponential increase with random matrix
  #using credible intervals
  expect_true({
    
    #Seed
    set.seed(234)
    
    #Simulation values
    params     <- c(mu= 0.001, sd = 0.2)
    X          <- matrix(c(N=10), nrow = 1)
    pfun       <- function(t,X,params){ cbind( rlnorm(nrow(X),params[1],params[2])* X[,1]) }
    v          <- matrix( c(+1), ncol=1)
    tmin       <- 0
    tmax       <- 1
    nsim       <- 1000
    conf       <- 95 #Credibility (%)
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, plot.sim = F)
    
    #Matrix for saving the results at specific times
    times      <- seq(tmin, tmax, length.out = 100)
    
    #Mean and sd of random parameter in pfun
    pmedian    <- exp(params[1])
    pvar       <- exp(2*params[1] + params[2]^2)*(exp(params[2]^2) - 1)
    
    #Function for approximate theoretical mean
    tmean      <- X[,1]*exp(pmedian*times)
      
    #Estimate mean and sd
    summarymat <- getMean(simulation, times)[[1]]
    simmean    <- apply(summarymat, 1, mean)       #Estimate mean
    qvars      <- apply(summarymat, 1, function(x) quantile(x,c(0.025,0.975)))         #Estimate variance
    
    #Plot
    #lines(times, tmean,  lwd = 2 )
    #lines(times, simmean, col ="green", lwd = 2)
    #lines(times, qvars[1,], col = "blue", lwd = 2)
    #lines(times, qvars[2,], col = "blue", lwd = 2)
    
    #Check the confidence interval matches conf times
    rejected <- length(which(tmean > qvars[2,] || tmean < qvars[1,]))
    
    #Check that we are rejecting less than conf%
    (rejected/length(res)) < (100 - conf)/100
    
  })
  
  
})