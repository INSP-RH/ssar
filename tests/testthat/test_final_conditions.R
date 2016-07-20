#Tests for checking that time = 0 and time = tmax are included int he
#results. This tests also check that the number of saved iterations are adecquate

context("Check final format")

test_that("Errors for inadecquate outputs",{
  
  #Check that initial and final points of simulation are included
  expect_true({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                        (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    maxiter    <- sample(1:1000,1)
    nsim       <- 5
    simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, maxiter = maxiter, plot.sim = F, kthsave = 2)
    
    #Check that 0 and maxiter are included
    all(c(0, maxiter) %in% simulation$Iteration)
    
  })
  
  #Check that initial and final points of simulation are included
  expect_true({

    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    tmax       <- runif(1,0.1,0.9)
    nsim       <- 10
    simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, tmax = tmax, plot.sim = F, kthsave = 2)

    #Check that 0 and maxiter are included
    all(simulation$Time[which(simulation$Iteration == max(simulation$Iteration))] >= tmax)

  })
  
  #Check that data frame is simulation
  expect_true({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    maxiter    <- 10
    nsim       <- 10
    simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, tmax = tmax, plot.sim = F, kthsave = 2)
    
    #Check that 0 and maxiter are included
    is.data.frame(simulation)
    
  })
  
  #Check that colnames match
  expect_equal({
    
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
    nsim       <- 10
    maxiter    <- 300

    #Run the program
    simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, maxiter = maxiter,
                       print.time = FALSE, plot.sim = FALSE)
    
    #Check that length matches number of variables
    ncol(simulation)

    
  }, 5)
  
})