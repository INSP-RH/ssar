#Tests for checking that time = 0 and time = tmax are included int he
#results. This tests also check that the number of saved iterations are adecquate

context("Checking simulation results using deSolve")

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
test_that("Simulation results using deSolve",{
  
  #Evaluating if mean of exponential increase with time-dependent function
  expect_true({
    
    #Seed
    set.seed(6248)

    #Simulation values
    params     <- c(b=0.15)
    X          <- matrix(10, nrow = 1)
    pfun       <- function(t,X,params){ as.matrix( (params[1] + 0.1*cos(t*pi)) * X[,1]) }
    v          <- matrix( c(+1), ncol=1)
    tmin       <- 0
    tmax       <- 10
    nsim       <- 1000
    conf       <- 95
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, plot.sim = F)

    #Matrix for saving the results at specific times
    times      <- seq(tmin, tmax, length.out = 100)

    
    #Solve with ode (this throws "unlock_solver" not resolved from current namespace (deSolve) when using devtools)
    #http://stackoverflow.com/questions/18193556/what-is-the-cause-of-error-in-cunlock-solver-error-in-desolve-package
    tryCatch({
      
      #Call deSolve
      library(deSolve)
      
      #Function for approximate theoretical mean 
      meanfun    <- function(Time, State, Pars){
        with(as.list(c(State, Pars)),{
          
          #Get cosine constant
          coscons <- cos(Time*pi)
          
          #Get equation
          dN <- (b + 0.1*coscons) * N
          
          return(list(c(dN)))
          
        })
      }
      
      vals <- as.data.frame(ode(c(N = X[,1]), times, meanfun, params))
    },
    error = function(e){
      
      #Read from already processed data
      load("Test_deSolve1.Rda")
      
      vals <<- vals
      return("Loading file...")
      
    })
    
    #Estimate mean and sd
    summarymat <- getMean(simulation, times)[[1]]
    simmean    <- apply(summarymat, 1, mean)       #Estimate mean
    qvars      <- apply(summarymat, 1, function(x) quantile(x,c(0.025,0.975)))         #Estimate variance

    #Do simean - tmean
    res <- abs(simmean - vals$N)

    #Confidence Z
    Z <- qnorm(1-(1-conf/100)/2)

    #Plot
    #lines(vals$time,vals$N, lwd = 2, type = "l", col = "green")
    #lines(times,simmean, col ="black", lwd = 2)
    #lines(times, qvars[1,], col = "blue", lwd = 2)
    #lines(times, qvars[2,], col = "blue", lwd = 2)
    
    
    #Check the confidence interval matches conf times
    rejected <- length(which(vals$N > qvars[2,] || vals$N < qvars[1,]))

    #Check that we are rejecting none
    (rejected/length(res)) <= 0

  })
  
  #Evaluating Lotka-Volterra model with time-dependent function
  expect_true({
    
    #Seed
    set.seed(937299)

    #Get initial parameters
    params     <- c(a = 3, b = 0.01, c = 2)

    #Notice that X must be inputed as matrix
    X          <- matrix(c(100, 100), ncol = 2)

    #Notice that pfun returns a matrix
    pfun       <- function(t, X, params){ cbind(params[1]*t*X[,1] + 1,
                                                 params[2]*X[,1]*X[,2],
                                                   params[3]*X[,2]) }

    #Propensity matrix
    v          <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)

    #Additional variables
    tmin       <- 0
    nsim       <- 100
    tmax       <- 3

    #Run the program
    simulation <- ssa(X, pfun, v, params, tmin, tmax = tmax, nsim = nsim, plot.sim = F)
    
    #Matrix for saving the results at specific times
    times      <- seq(tmin, tmax, length.out = 100)
    

    
    #Solve with ode (this throws "unlock_solver" not resolved from current namespace (deSolve) when using devtools)
    #http://stackoverflow.com/questions/18193556/what-is-the-cause-of-error-in-cunlock-solver-error-in-desolve-package
    tryCatch({
      
      #Call deSolve
      library(deSolve)
      
      #Function for approximate theoretical mean 
      meanfun    <- function(Time, State, Pars){
        with(as.list(c(State, Pars)),{
          
          #Get cosine constant
          consa <- a*Time
          
          #Get equation
          dX1 <- consa*X1 - b*X1*X2
          dX2 <- b*X1*X2  - c*X2
          
          return(list(c(dX1, dX2)))
          
        })
      }
      
      vals <- as.data.frame(ode(c(X1 = X[,1], X2 = X[,2]), times, meanfun, params))
    },
    error = function(e){
      
      #Read from already processed data
      load("Test_deSolve2.Rda")
      
      vals <<- vals
      return("Loading file...")
      
    })
    
    #Estimate mean and sd
    summarymat <- getMean(simulation, times)
    simX1      <- apply(summarymat[[1]], 1, mean)       #Estimate mean
    simX2      <- apply(summarymat[[2]], 1, mean)       #Estimate mean
    qvarsX1    <- apply(summarymat[[1]], 1, function(x) quantile(x,c(0.025,0.975)))         #Estimate variance
    qvarsX2    <- apply(summarymat[[2]], 1, function(x) quantile(x,c(0.025,0.975)))         #Estimate variance
    
    #Plot
    #lines(vals$time,vals$X1, lwd = 2, type = "l", col = "black")
    #lines(vals$time,vals$X2, lwd = 2, type = "l", col = "black")
    #lines(times,simX1, col ="blue", lwd = 2)
    #lines(times,simX2, col ="red", lwd = 2)
    #lines(times, qvarsX1[1,], col = "blue", lwd = 2)
    #lines(times, qvarsX1[2,], col = "blue", lwd = 2)
    #lines(times, qvarsX2[1,], col = "red", lwd = 2)
    #lines(times, qvarsX2[2,], col = "red", lwd = 2)
    
    #Check the confidence interval matches conf times
    rejectedX1 <- length(c(which(vals$X1 > qvarsX1[2,]), which(vals$X1 < qvarsX1[1,])))
    rejectedX2 <- length(c(which(vals$X2 > qvarsX2[2,]), which(vals$X2 < qvarsX2[1,])))
    rejected   <- rejectedX1 + rejectedX2
    
    #Check that we are rejecting none
    (rejected/length(res)) <= 0
    
  })

  #Lotka Volterra with propensity function that goes negative around time = 0.4
  expect_true({

    #Seed
    set.seed(937299)

    #Get initial parameters
    params     <- c(a = 3, b = 0.01, c = 2)

    #Notice that X must be inputed as matrix
    X          <- matrix(c(100, 100), ncol = 2)

    #Notice that pfun returns a matrix
    pfun       <- function(t, X, params){ cbind(rlnorm(nrow(X),log(params[1]),0.1)*t*X[,1] + 1,
                                                rlnorm(nrow(X),log(params[2]),0.1)*cos(t*2*pi)*X[,1]*X[,2],
                                                rlnorm(nrow(X),log(params[3]),0.1)*sin(t*2*pi)*X[,2]) }

    #Propensity matrix
    v          <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)

    #Additional variables
    tmin       <- 0
    nsim       <- 1000
    tmax       <- 1

    #Run the program
    simulation <- ssa(X, pfun, v, params, tmin, tmax = tmax, nsim = nsim, plot.sim = F, print.time = F)
    
    #Check that program ended before tmax
    if (min(simulation$Time) < tmax){
      TRUE
    } else {
      FALSE
    }

  })
  
})