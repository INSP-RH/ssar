#Test conditions for inputs

context("Initial conditions")

test_that("Errors for inadecquate inputs",{
  
  #Check that nsim is positive
  expect_warning({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                        (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    tmax       <- 1
    nsim       <- runif(1,-10000000,1)
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, maxiter = 2, plot.sim = F,
                       title = "Time dependent Logistic Growth example",
                       xlab = "Time", ylab = "N")

  }, paste("nsim is not a strictly positive integer.",
           "Defaulting to closest positive integer"))
  
  #Check that kthsave is positive
  expect_warning({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    tmax       <- 1
    nsim       <- 10
    kthsave     <- runif(1,-10000000,0.99)
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, maxiter = 2, kthsave = kthsave, plot.sim = F)
    
  }, paste("kthsave is not a strictly positive integer.",
           "Defaulting to closest positive integer"))
  
  #Check that maxiter is integer
  expect_warning({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    tmax       <- 1
    nsim       <- 10
    maxiter    <- pi
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, maxiter = maxiter, plot.sim = F)
    
  }, paste("maxiter is not an integer.",
           "Defaulting to maxiter =", ceiling(maxiter)))
  
  #Check that tmin > tmax
  expect_error({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 1
    tmax       <- 0
    nsim       <- 10
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, maxiter = 2, plot.sim = F)
    
  },"tmin >= tmax")
  
  #Check that maxiter > 0
  expect_error({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    tmax       <- 1
    nsim       <- 10
    maxiter    <- runif(1,-1000,0.99)
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, maxiter = maxiter, plot.sim = F)
    
  },"Maximum iterations < 0")
  
  #Check that maxtime > 0
  expect_error({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    tmax       <- 1
    nsim       <- 10
    maxtime    <- runif(1,-1000,0.99)
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, maxtime = maxtime, plot.sim = F)
    
  },"Maximum time       < 0")
  
  #Check that maxiter < Inf, maxtime < Inf and tmax < Inf (at same time)
  expect_warning({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    nsim       <- 10
    maxtime    <- runif(1,-1000,0.99)
    simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, plot.sim = F)
    
  },paste("Please specify either maxiter or maxtime or tmax", 
          "to avoid infinite loops. Defaulting tmax =",tmin + 1))
  
  #Check that pfun is matrix
  expect_error({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- matrix(c(N=500), nrow = 1)
    pfun       <- function(t,X,params){ c(params[1] *(1 + sin(t))* X[,1],
                                              (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    tmax       <- 1
    nsim       <- 10
    maxiter    <- 2
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, maxiter = maxiter, plot.sim = F)
    
  },"pfun needs to be a matrix valued function")
  
  #Check that xinit is matrix
  expect_error({
    
    params     <- c(b=2, d=1, k=1000)
    X          <- c(N=500)
    pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[1],
                                          (params[2] + (params[1]-params[2])*X[1]/params[3])*X[1]) }
    v          <- matrix( c(+1, -1),ncol=2)
    tmin       <- 0
    tmax       <- 1
    nsim       <- 10
    maxiter    <- 2
    simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, maxiter = maxiter, plot.sim = F)
    
  },"xinit needs to be a matrix")
  
  
})