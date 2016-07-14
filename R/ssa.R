#' @title Gillespie Stochastic Simulation Algorithm
#'
#' @description Implementation of Gillespie's Stochastic Simulation Algorithm for in-homogeneous
#' continuous-time Markov Chains.
#'
#' @param tmin       ( double  ) Initial .time for the simulation
#' @param nsim       ( integer ) Number of simulations
#' @param v          ( matrix  ) Propensity Matrix
#' @param pfun       ( vector  ) Character vector of propensity functions where state variables
#'                                correspond to the names of the elements in xinit.
#' @param xinit      ( matrix  ) Numeric named vector with initial variables.
#' @param params     (  list   ) List of parameters including functions used by the problem
#'                                that are not contained in base R
#'
#' \strong{ Optional }
#' @param tmax       ( double  ) Final .time for the simulation
#' @param print.time ( boolean ) Value indicating whether to print current estimation .time
#' @param maxiter    ( numeric ) Maximum number of iterations for model
#' @param maxtime    ( numeric ) Maximum .time (in seconds) before breaking the proces
#' @param kthsave    ( numeric ) Integer indicating iteration-lapse to save iterations (saves every kthsave iteration)
#' @param keep.file  ( boolean ) Instruction whether to keep temporary file created by the process
#' @param file.only  ( boolean ) Instruction whether to return the information as a txt file when \code{TRUE}.
#'                                this makes the package run faster. 
#' @param fname      ( string  ) Name of file where information will be saved if \code{keep.file == T}
#' @param plot.sim   ( boolean ) Variable indicating whether to plot the results
#' @param title      ( string  ) Title for plot
#' @param xlab       ( string  ) X-axis label for plot
#' @param ylab       ( string  ) Y-axis label for plot
#'
#'
#' @author Rodrigo Zepeda
#' @author Dalia Camacho
#'
#' @import compiler
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices rainbow
#' @importFrom graphics lines plot
#' @importFrom utils read.table
#' @note To include .time dependent functions please express them as matrix functions of t. For example
#' \code{pfun <- function(t){cbind(sin(t),2*t)}}
#'
#' @references INCLUDE REFERENCES HERE TO GILLESPIE SSA WITH TIME-DEPENDENT PARAMETERS
#'
#' @examples
#'
#' #EXAMPLE 1
#' #------------------------
#' #Logistic Growth
#' set.seed(123)
#' params     <- c(b=2, d=1, k=1000)
#' X          <- matrix(c(N=500), nrow = 1)
#' pfun       <- function(t,X,params){ cbind(params[1] * X[,1], 
#'                        (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
#' v          <- matrix( c(+1, -1),ncol=2)
#' tmin       <- 0
#' tmax       <- 1
#' nsim       <- 5
#' simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, 
#'                    title = "Logistic Growth example", xlab = "Time", ylab = "N")
#' 
#' #EXAMPLE 2
#' #------------------------
#' #Time dependent logistic growth
#' set.seed(123)
#' params     <- c(b=2, d=1, k=1000)
#' X          <- matrix(c(N=500), nrow = 1)
#' pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1], 
#'                     (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
#' v          <- matrix( c(+1, -1),ncol=2)
#' tmin       <- 0
#' tmax       <- 1
#' nsim       <- 5
#' simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, 
#'                    title = "Time dependent Logistic Growth example", 
#'                    xlab = "Time", ylab = "N")
#'
#' #EXAMPLE 3
#' #------------------------
#' #Lotka Volterra
#' #Set seed
#' set.seed(123)
#'
#' #Get initial parameters
#' params     <- c(a = 3, b = 0.01, c = 2)
#'
#' #Notice that X must be inputed as matrix
#' X          <- matrix(c(100, 100), ncol = 2)
#'
#' #Notice that pfun returns a matrix
#' pfun       <- function(t, X, params){ cbind(params[1]*t*X[,1] + 1, 
#'                                              params[2]*X[,1]*X[,2], 
#'                                                params[3]*X[,2]) }
#'
#' #Propensity matrix
#' v          <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)
#'
#' #Additional variables
#' tmin       <- 0
#' nsim       <- 10
#' maxiter    <- 300
#'
#' #Run the program
#' simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, maxiter = maxiter, 
#'                    print.time = TRUE, title = "Lotka-Volterra example", xlab = "Time",
#'                     ylab = "Individuals")
#'
#'
#' #EXAMPLE 4
#' #------------------------
#' #Time dependent SIS model
#' set.seed(123)
#'
#' #Initial parameters
#' k          <-  24576.5529836797
#' delta      <-  0.0591113454895868 + 0.208953907151055
#' gamma_ct   <-  0.391237630231631
#' params     <- c(k = k, delta = delta, gamma_ct = gamma_ct)
#' X          <- matrix(c(S = 1000, I = 40), ncol = 2)
#' pfun       <- function(t, X, params){
#'
#'   #Value to return
#'   matreturn  <- matrix(NA, nrow = length(t), ncol = 6)
#'
#'   #Create birth function
#'   lambda     <- function(t){ return(4.328e-4 - (2.538e-7)*t - 
#'                              (3.189e-7)*sin(2 * t * pi/52) - 
#'                              (3.812e-7)*cos(2 * t * pi/52) ) }
#'
#'   #Create death function
#'   mu         <- function(t){ return(9.683e-5 + (1.828e-8)*t + 
#'                              (2.095e-6)*sin(2 * t * pi/52) - 
#'                              (8.749e-6)*cos(2 * t * pi/52))}
#'
#'   #Create infectives function
#'   beta_fun   <- function(t){ return( 0.479120824267286 + 
#'                              0.423263042762498*sin(-2.82494252560096 + 2*t*pi/52) )}
#'
#'   #Estimate values
#'   matreturn[,1] <- lambda(t)*(X[,1] + X[,2])
#'   matreturn[,2] <- mu(t)*X[,1]
#'   matreturn[,3] <- beta_fun(t)*X[,1]*X[,2]/(1 + params[1]*X[,2])
#'   matreturn[,4] <- mu(t)*X[,2]
#'   matreturn[,5] <- params[2]*X[,2]
#'   matreturn[,6] <- params[3]*X[,2]
#'
#'   #Return
#'   return(matreturn)
#'
#' }
#' v          <- matrix(c(1,-1, -1, 0, 0, 1, 0, 0, 1, -1, -1, -1), nrow = 2, byrow = TRUE)
#' tmin       <- 0
#' tmax       <- 0.5
#' nsim       <- 100
#'
#' #Simulate the values
#' simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim = nsim, print.time = FALSE, 
#'                    plot.sim = FALSE, kthsave = 1)
#'
#' #Plot using ggplot2 library
#' \dontrun{
#' library(ggplot2)
#' ggplot(data = simulation, aes(x = Time, y = Var2, group=as.factor(Simulation))) +
#'    geom_line(aes(color = as.factor(Simulation))) + theme_bw() + 
#'    theme(legend.position="none") +
#'    ggtitle(paste0("SIS example; Infected cases ", nsim, " simulations")) + 
#'    xlab("Time") + ylab("Individuals") +
#'    geom_vline(xintercept = tmax, linetype = 2)
#'}
#'
#' #EXAMPLE 5
#' #------------------------------------
#' set.seed(123)
#'
#' #Start the parameters
#' params     <- c(NA)
#' X          <- matrix(c(N = 10), nrow = 1)
#' 
#' #Notice the cbind forces to return matrix
#' pfun       <- function(t, X, params){ cbind(1.1 + sin(pi*t/0.01))*X[,1] } 
#' v          <- matrix( c(+1), ncol=1)
#' tmin       <- 0
#' nsim       <- 25
#' maxtime    <- 1   #Maximum .time for model: 1 seconds and stop iterating
#' simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, maxtime = maxtime, 
#'                    print.time = FALSE, 
#'                    title = "Start of an epidemic with time-dependent R0", 
#'                    xlab = "Time", ylab = "Individuals")
#'
#' @useDynLib ssar
#'
#' @export
#'

ssa <- function(xinit, pfun, v, params = c(), tmin = 0, tmax = Inf, nsim = 10,
                 maxiter = Inf, maxtime = Inf, print.time = FALSE,
                 plot.sim = TRUE, title = "", xlab = "", ylab = "",
                 kthsave = 1, keep.file = FALSE, file.only = FALSE,
                 fname = "Temporary_File_ssa.txt"){

  #-----------------------------------------
  #PREPARATION: HERE WE CHECK PARAMETERS MAKE SENSE
  #----------------------------------------

  #Print to user
  #cat("Starting the process...\n")

  #Check that nsim is integer and positive
  if ( (nsim != ceiling(nsim) & nsim != floor(nsim)) || nsim <= 1){
    warning(paste("nsim is not a strictly positive integer.",
                  "Defaulting to closest positive integer"))
    nsim   <- max( ceiling(nsim), 2)
  }

  #Check that kthsave is integer and positive
  if ( (kthsave != ceiling(kthsave) & kthsave != floor(kthsave))
       || kthsave < 1){
    warning(paste("kthsave is not a strictly positive integer.",
                    "Defaulting to closest positive integer"))
    kthsave   <- max( ceiling(kthsave), 1)
  }

  #Check that tmin < tmax and maxiter, maxtime are > 0
  if (tmin >= tmax){
    stop("tmin >= tmax")
  }
  if (maxiter <= 0){
    stop("Maximum iterations < 0")
  }
  if (maxtime <= 0){
    stop("Maximum time       < 0")
  }

  #Value indicating whether to check maxiter or maxtime
  .opts_check  <- FALSE

  #Compile pfun to make it faster
  .pfun     <- cmpfun(pfun)
  
  #Check that pfun is matrix
  .evalpfun <- pfun(tmin,xinit,params)
  if (!identical(.evalpfun, as.matrix(.evalpfun))){
    stop("pfun needs to be a matrix valued function")
  }
  
  #Check that xinit is matrix
  if (!identical(xinit, as.matrix(xinit))){
    stop("xinit needs to be a matrix")
  }
  
  #Length of pfun vector (as it is in string form)
  .lfun <- length(.evalpfun)

  #Check whether additional options were specified
  if (maxiter < Inf || maxtime < Inf || print.time){
    .opts_check <- TRUE
  }
  
  #Check that params are given
  if (length(params) == 0){
    params   <- rep(1,2)
  }
  
  #Don't allow for infinite loops
  if (maxiter == Inf & maxtime == Inf & tmax == Inf){
    tmax     <- tmin + 1
  }
  
  #Keep file if file.only option is on
  if (file.only){
    keep.file <- TRUE
    plot.sim  <- FALSE
  }

  #-----------------------------------------
  #SIMULATION USING C++ PROGRAM (RCPP)
  #-----------------------------------------

  #Run C program to estimate the loop
  #cat("Running simulation...\n")

  ssa_loop(xinit, .pfun, v, params, tmin, nsim, .lfun,
            tmax, maxiter, maxtime,  print.time, .opts_check, kthsave)

  #-----------------------------------------
  #R AFTERMATH
  #-----------------------------------------

  #Read temporary file and delete
  #cat("Getting information...\n")

  if (!file.only){  
  
  #Try to open file which is in current directory
  .filename  <- "Temporary_File_ssa.txt"
  .datamodel <- as.data.frame(read.table(.filename, header = T))
  
  }
  
  #Delete file from memory if indication is given
  if (!keep.file){
    
    try({
      file.remove("Temporary_File_ssa.txt")
    })
    
  } else {
    
    try({
      file.rename("Temporary_File_ssa.txt", fname)
    })
  }

  #-----------------------------------------
  #PLOT
  #-----------------------------------------

  #Plot if plot indication is given
  if (plot.sim){
    
    #Get plot limits
    .mvals   <- length(xinit)
    .lowval  <- min(.datamodel[, 4:(.mvals + 3)])
    .highval <- max(.datamodel[, 4:(.mvals + 3)])
    .ylims   <- c( .lowval, .highval)
    
    #Check finite time
    .posvals <- which(.datamodel$Time < Inf)
    .xlims   <- c( min(.datamodel$Time[.posvals]), max(.datamodel$Time[.posvals]))

    #Colors will be assigned to each variable
    colors <- rainbow(.mvals)

    #Create the plot
    for (i in 1:nsim){

      #Get one simulation at a time 
      .subsimulation <- subset(.datamodel, Simulation == i & Time < Inf)

      if (i == 1){

        plot(.subsimulation$Time, .subsimulation[, 4], type = "l",
             col = colors[1], xlim = .xlims, ylim = .ylims, main = title,
             xlab = xlab, ylab = ylab)

      } else{

        lines(.subsimulation$Time, .subsimulation[, 4], col = colors[1])

      }

      if (.mvals > 1){
        for (j in 2:.mvals){
          lines(.subsimulation$Time, .subsimulation[, 3 + j], col = colors[j])
        }
      }
    }

    }

  #Rerutn as data frame
  if (file.only){
    return(TRUE)
  } else {
    return(.datamodel)
  }
    
  

}
#' @title Gillespie Stochastic Simulation Algorithm
#'
#' @description Implementation of Gillespie's Stochastic Simulation Algorithm for in-homogeneous
#' continuous-time Markov Chains.
#'
#' @param tmin       ( double  ) Initial .time for the simulation
#' @param nsim       ( integer ) Number of simulations
#' @param v          ( matrix  ) Propensity Matrix
#' @param pfun       ( vector  ) Character vector of propensity functions where state variables
#'                                correspond to the names of the elements in xinit.
#' @param xinit      ( matrix  ) Numeric named vector with initial variables.
#' @param params     (  list   ) List of parameters including functions used by the problem
#'                                that are not contained in base R
#'
#' \strong{ Optional }
#' @param tmax       ( double  ) Final .time for the simulation
#' @param print.time ( boolean ) Value indicating whether to print current estimation .time
#' @param maxiter    ( numeric ) Maximum number of iterations for model
#' @param maxtime    ( numeric ) Maximum .time (in seconds) before breaking the proces
#' @param kthsave    ( numeric ) Integer indicating iteration-lapse to save iterations (saves every kthsave iteration)
#' @param keep.file  ( boolean ) Instruction whether to keep temporary file created by the process
#' @param file.only  ( boolean ) Instruction whether to return the information as a txt file when \code{TRUE}.
#'                                this makes the package run faster. 
#' @param fname      ( string  ) Name of file where information will be saved if \code{keep.file == T}
#' @param plot.sim   ( boolean ) Variable indicating whether to plot the results
#' @param title      ( string  ) Title for plot
#' @param xlab       ( string  ) X-axis label for plot
#' @param ylab       ( string  ) Y-axis label for plot
#'
#'
#' @author Rodrigo Zepeda
#' @author Dalia Camacho
#'
#' @import compiler
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices rainbow
#' @importFrom graphics lines plot
#' @importFrom utils read.table
#' @note To include .time dependent functions please express them as matrix functions of t. For example
#' \code{pfun <- function(t){cbind(sin(t),2*t)}}
#'
#' @references INCLUDE REFERENCES HERE TO GILLESPIE SSA WITH TIME-DEPENDENT PARAMETERS
#'
#' @examples
#'
#' #EXAMPLE 1
#' #------------------------
#' #Logistic Growth
#' set.seed(123)
#' params     <- c(b=2, d=1, k=1000)
#' X          <- matrix(c(N=500), nrow = 1)
#' pfun       <- function(t,X,params){ cbind(params[1] * X[,1], 
#'                        (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
#' v          <- matrix( c(+1, -1),ncol=2)
#' tmin       <- 0
#' tmax       <- 1
#' nsim       <- 5
#' simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, 
#'                    title = "Logistic Growth example", xlab = "Time", ylab = "N")
#' 
#' #EXAMPLE 2
#' #------------------------
#' #Time dependent logistic growth
#' set.seed(123)
#' params     <- c(b=2, d=1, k=1000)
#' X          <- matrix(c(N=500), nrow = 1)
#' pfun       <- function(t,X,params){ cbind(params[1] *(1 + sin(t))* X[,1], 
#'                     (params[2] + (params[1]-params[2])*X[,1]/params[3])*X[,1]) }
#' v          <- matrix( c(+1, -1),ncol=2)
#' tmin       <- 0
#' tmax       <- 1
#' nsim       <- 5
#' simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim, 
#'                    title = "Time dependent Logistic Growth example", 
#'                    xlab = "Time", ylab = "N")
#'
#' #EXAMPLE 3
#' #------------------------
#' #Lotka Volterra
#' #Set seed
#' set.seed(123)
#'
#' #Get initial parameters
#' params     <- c(a = 3, b = 0.01, c = 2)
#'
#' #Notice that X must be inputed as matrix
#' X          <- matrix(c(100, 100), ncol = 2)
#'
#' #Notice that pfun returns a matrix
#' pfun       <- function(t, X, params){ cbind(params[1]*t*X[,1] + 1, 
#'                                              params[2]*X[,1]*X[,2], 
#'                                                params[3]*X[,2]) }
#'
#' #Propensity matrix
#' v          <- matrix(c(+1,-1,0,0,+1,-1),nrow=2,byrow=TRUE)
#'
#' #Additional variables
#' tmin       <- 0
#' nsim       <- 10
#' maxiter    <- 300
#'
#' #Run the program
#' simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, maxiter = maxiter, 
#'                    print.time = TRUE, title = "Lotka-Volterra example", xlab = "Time",
#'                     ylab = "Individuals")
#'
#'
#' #EXAMPLE 4
#' #------------------------
#' #Time dependent SIS model
#' set.seed(123)
#'
#' #Initial parameters
#' k          <-  24576.5529836797
#' delta      <-  0.0591113454895868 + 0.208953907151055
#' gamma_ct   <-  0.391237630231631
#' params     <- c(k = k, delta = delta, gamma_ct = gamma_ct)
#' X          <- matrix(c(S = 1000, I = 40), ncol = 2)
#' pfun       <- function(t, X, params){
#'
#'   #Value to return
#'   matreturn  <- matrix(NA, nrow = length(t), ncol = 6)
#'
#'   #Create birth function
#'   lambda     <- function(t){ return(4.328e-4 - (2.538e-7)*t - 
#'                              (3.189e-7)*sin(2 * t * pi/52) - 
#'                              (3.812e-7)*cos(2 * t * pi/52) ) }
#'
#'   #Create death function
#'   mu         <- function(t){ return(9.683e-5 + (1.828e-8)*t + 
#'                              (2.095e-6)*sin(2 * t * pi/52) - 
#'                              (8.749e-6)*cos(2 * t * pi/52))}
#'
#'   #Create infectives function
#'   beta_fun   <- function(t){ return( 0.479120824267286 + 
#'                              0.423263042762498*sin(-2.82494252560096 + 2*t*pi/52) )}
#'
#'   #Estimate values
#'   matreturn[,1] <- lambda(t)*(X[,1] + X[,2])
#'   matreturn[,2] <- mu(t)*X[,1]
#'   matreturn[,3] <- beta_fun(t)*X[,1]*X[,2]/(1 + params[1]*X[,2])
#'   matreturn[,4] <- mu(t)*X[,2]
#'   matreturn[,5] <- params[2]*X[,2]
#'   matreturn[,6] <- params[3]*X[,2]
#'
#'   #Return
#'   return(matreturn)
#'
#' }
#' v          <- matrix(c(1,-1, -1, 0, 0, 1, 0, 0, 1, -1, -1, -1), nrow = 2, byrow = TRUE)
#' tmin       <- 0
#' tmax       <- 0.5
#' nsim       <- 100
#'
#' #Simulate the values
#' simulation <- ssa(X, pfun, v, params, tmin, tmax, nsim = nsim, print.time = FALSE, 
#'                    plot.sim = FALSE, kthsave = 1)
#'
#' #Plot using ggplot2 library
#' \dontrun{
#' library(ggplot2)
#' ggplot(data = simulation, aes(x = Time, y = Var2, group=as.factor(Simulation))) +
#'    geom_line(aes(color = as.factor(Simulation))) + theme_bw() + 
#'    theme(legend.position="none") +
#'    ggtitle(paste0("SIS example; Infected cases ", nsim, " simulations")) + 
#'    xlab("Time") + ylab("Individuals") +
#'    geom_vline(xintercept = tmax, linetype = 2)
#'}
#'
#' #EXAMPLE 5
#' #------------------------------------
#' set.seed(123)
#'
#' #Start the parameters
#' params     <- c(NA)
#' X          <- matrix(c(N = 10), nrow = 1)
#' 
#' #Notice the cbind forces to return matrix
#' pfun       <- function(t, X, params){ cbind(1.1 + sin(pi*t/0.01))*X[,1] } 
#' v          <- matrix( c(+1), ncol=1)
#' tmin       <- 0
#' nsim       <- 25
#' maxtime    <- 1   #Maximum .time for model: 1 seconds and stop iterating
#' simulation <- ssa(X, pfun, v, params, tmin, nsim = nsim, maxtime = maxtime, 
#'                    print.time = FALSE, 
#'                    title = "Start of an epidemic with time-dependent R0", 
#'                    xlab = "Time", ylab = "Individuals")
#'
#' @useDynLib ssar
#'
#' @export
#'

ssa <- function(xinit, pfun, v, params = c(), tmin = 0, tmax = Inf, nsim = 10,
                 maxiter = Inf, maxtime = Inf, print.time = FALSE,
                 plot.sim = TRUE, title = "", xlab = "", ylab = "",
                 kthsave = 1, keep.file = FALSE, file.only = FALSE,
                 fname = "Temporary_File_ssa.txt"){

  #-----------------------------------------
  #PREPARATION: HERE WE CHECK PARAMETERS MAKE SENSE
  #----------------------------------------

  #Print to user
  #cat("Starting the process...\n")

  #Check that nsim is integer and positive
  if ( (nsim != ceiling(nsim) & nsim != floor(nsim)) || nsim <= 1){
    warning(paste("nsim is not a strictly positive integer.",
                  "Defaulting to closest positive integer"))
    nsim   <- max( ceiling(nsim), 2)
  }

  #Check that kthsave is integer and positive
  if ( (kthsave != ceiling(kthsave) & kthsave != floor(kthsave))
       || kthsave < 1){
    warning(paste("kthsave is not a strictly positive integer.",
                    "Defaulting to closest positive integer"))
    kthsave   <- max( ceiling(kthsave), 1)
  }

  #Check that tmin < tmax and maxiter, maxtime are > 0
  if (tmin >= tmax){
    stop("tmin >= tmax")
  }
  if (maxiter <= 0){
    stop("Maximum iterations < 0")
  }
  if (maxtime <= 0){
    stop("Maximum time       < 0")
  }

  #Value indicating whether to check maxiter or maxtime
  .opts_check  <- FALSE

  #Compile pfun to make it faster
  .pfun     <- cmpfun(pfun)
  
  #Check that pfun is matrix
  .evalpfun <- pfun(tmin,xinit,params)
  if (!identical(.evalpfun, as.matrix(.evalpfun))){
    stop("pfun needs to be a matrix valued function")
  }
  
  #Check that xinit is matrix
  if (!identical(xinit, as.matrix(xinit))){
    stop("xinit needs to be a matrix")
  }
  
  #Length of pfun vector (as it is in string form)
  .lfun <- length(.evalpfun)

  #Check whether additional options were specified
  if (maxiter < Inf || maxtime < Inf || print.time){
    .opts_check <- TRUE
  }
  
  #Check that params are given
  if (length(params) == 0){
    params   <- rep(1,2)
  }
  
  #Don't allow for infinite loops
  if (maxiter == Inf & maxtime == Inf & tmax == Inf){
    tmax     <- tmin + 1
  }
  
  #Keep file if file.only option is on
  if (file.only){
    keep.file <- TRUE
    plot.sim  <- FALSE
  }

  #-----------------------------------------
  #SIMULATION USING C++ PROGRAM (RCPP)
  #-----------------------------------------

  #Run C program to estimate the loop
  #cat("Running simulation...\n")

  ssa_loop(xinit, .pfun, v, params, tmin, nsim, .lfun,
            tmax, maxiter, maxtime,  print.time, .opts_check, kthsave)

  #-----------------------------------------
  #R AFTERMATH
  #-----------------------------------------

  #Read temporary file and delete
  #cat("Getting information...\n")

  if (!file.only){  
  
  #Try to open file which is in current directory
  .filename  <- "Temporary_File_ssa.txt"
  .datamodel <- as.data.frame(read.table(.filename, header = T))
  
  }
  
  #Delete file from memory if indication is given
  if (!keep.file){
    
    try({
      file.remove("Temporary_File_ssa.txt")
    })
    
  } else {
    
    try({
      file.rename("Temporary_File_ssa.txt", fname)
    })
  }

  #-----------------------------------------
  #PLOT
  #-----------------------------------------

  #Plot if plot indication is given
  if (plot.sim){
    
    #Get plot limits
    .mvals   <- length(xinit)
    .lowval  <- min(.datamodel[, 4:(.mvals + 3)])
    .highval <- max(.datamodel[, 4:(.mvals + 3)])
    .ylims   <- c( .lowval, .highval)
    
    #Check finite time
    .posvals <- which(.datamodel$Time < Inf)
    .xlims   <- c( min(.datamodel$Time[.posvals]), max(.datamodel$Time[.posvals]))

    #Colors will be assigned to each variable
    colors <- rainbow(.mvals)

    #Create the plot
    for (i in 1:nsim){

      #Get one simulation at a time 
      .subsimulation <- subset(.datamodel, Simulation == i & Time < Inf)

      if (i == 1){

        plot(.subsimulation$Time, .subsimulation[, 4], type = "l",
             col = colors[1], xlim = .xlims, ylim = .ylims, main = title,
             xlab = xlab, ylab = ylab)

      } else{

        lines(.subsimulation$Time, .subsimulation[, 4], col = colors[1])

      }

      if (.mvals > 1){
        for (j in 2:.mvals){
          lines(.subsimulation$Time, .subsimulation[, 3 + j], col = colors[j])
        }
      }
    }

    }

  #Rerutn as data frame
  if (file.only){
    return(TRUE)
  } else {
    return(.datamodel)
  }
    
  

}
