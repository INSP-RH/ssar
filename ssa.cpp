#include <Rcpp.h>
#include <ctime>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>  //include this to use string streams
#include <string> 
using namespace Rcpp;

//Function to print time variable in console
void showTime(double x) {
  Rcout << "Time = " << x << std::endl;
}

//Warning after passing maxiter iterations
void warnMaxiter(double maxiter) {
  Rcerr << "Program reached " << maxiter << " iterations and stopped" << std::endl;
}

//Warning after passing maxiter seconds
void warnTime(double maxtime) {
  Rcerr << "Program reached " << maxtime << " seconds and stopped" << std::endl;
}

//Warning if a0 < 0
void warnNegative(int iter) {
  Rcerr << "Stopped at iter = " << iter << " due to negative values in pfun" << std::endl;
}

//Method for estimating colSums (http://adv-r.had.co.nz/Rcpp.html)
std::pair<NumericVector, bool>  colSumsinv (NumericMatrix x,  int nrow,  int ncol){
  
  //Boolean variable checking all colSums >= 0
  bool flag = false;
  
  //Pair to return
  std::pair <NumericVector, bool> returnvals;
  
  //Create vector of size nrow
  NumericVector ColSum(nrow);
  
  //Loop through elements in row
  for( int row = 0; row < nrow; row++ )
  {
    //Start the total value
    double sumres = 0;
    
    //Loop through elements in column
    for( int column = 0; column < ncol; column++ )
    {
      
      //Add to total
      sumres += x(row, column);
      
    }
    
    //Change 0s to inf
    if (sumres == 0){
      
      //Update vector
      ColSum[row] = std::numeric_limits<double>::infinity();
      
    } else {
      
      //Update vector
      ColSum[row] = 1/sumres;  
      
    }
    
    
    //Check if all values are positive
    if (sumres < 0){
      flag = true;
    }
    
  }
  
  //Create pair to return
  returnvals = std::make_pair (ColSum,flag);
  
  return(returnvals);
}

//Method for getting the column where change occurred from a matrix
NumericMatrix getMatx (NumericMatrix matx, NumericMatrix v, NumericVector vec, int lfun, int nsim, int mvals){
  
  //Create variables
  int counter = 0;
  int k       = 0;
  bool flag   = true;
  
  //Matrix of changes
  NumericMatrix col(nsim, mvals);
  
  //Loop through all columns of matrix
  while (k < lfun && flag){
    
    //Matrix as cummulative (sum previous values)
    if (k > 0){
      matx( _,k) = matx( _,k-1) + matx( _,k);
    }
    
    //Inner loop
    int j       = 0;
    
    //Loop through values of vec to check if is smaller
    while (flag && j < nsim){
      
      //Check if value is smaller
      if (vec(j) < matx(j,k)){
        
        //Assign to opt vector
        col(j,_) = v(_,k);
        vec(j)   = std::numeric_limits<double>::infinity(); //So that they are not considered in next loop
        
        //Counter to break from loop
        counter += 1;
        
      } //Close if
      
      if (counter == nsim){
        flag = false;
      }
      
      //Update k
      j += 1;
      
    } //Close while
    
    //Update k
    k += 1;
    
  } //Close while

  return col;
  
} //Close function

//Method for matrix vector  product (http://stackoverflow.com/questions/24933290/elementwise-matrix-multiplication-r-versus-rcpp-how-to-speed-this-code-up)
NumericMatrix MatVecHadamard(NumericMatrix X, NumericVector y,  int nrowX,  int ncolX){
  int counter = 0;
  for ( int j=0; j < ncolX; j++) {
    for ( int i=0; i < nrowX; i++)  {
      X[counter++] *= y[i];
    }
  }
  return X;
}

// [[Rcpp::export]]
void ssa_loop (NumericVector X, Function pfun, NumericMatrix v, NumericVector params, double tmin,
                         int nsim, int lfun, double tmax, double maxiter, double maxtime, bool printtime, bool optscheck, int kthsave){
  
  //----------------------------------------------------------------
  //INICIALIZE VARIABLES
  //----------------------------------------------------------------
  
  
  //One dimensional
  //----------------------------------------------------------------
  int mvals           = X.size(); //Amount of variables
  double tf           = tmin;    //Initial time for loop
  
  //Vectors
  //----------------------------------------------------------------
  NumericVector t(nsim, tmin);   //For time
  NumericVector tau(nsim);       //For time jumps
  NumericVector a0(nsim);        //For rate of jumps
  NumericVector vec(nsim);       //For random numbers
  IntegerVector nsimRange(nsim); //For indexing simulation value
  IntegerVector indicator(nsim); // For saving index range
  
  //Matrices
  //----------------------------------------------------------------
  NumericMatrix a(nsim, lfun);              //For probability of change
  NumericMatrix matx(nsim, lfun);           //Matrix for probability of change
  NumericMatrix col(nsim, mvals);           //For changes from propensity function
  NumericMatrix Xvals(nsim, mvals);         //For saving the results
  NumericMatrix Printable(nsim, mvals + 3); //For saving the results
  
  
  //Other
  //----------------------------------------------------------------
  time_t  loopinit    = time(0); //Time counter
  time_t  current;
  std::pair<NumericVector, bool>   a0vals;
  
  //Write file
  //----------------------------------------------------------------
  std::ofstream myfile;
  myfile.open ("Temporary_File_ssa.txt");
  
  //Create header of file
  std::string headerfile = "Simulation Iteration Time";
  
  
  //http://stackoverflow.com/questions/5290089/how-to-convert-a-number-to-string-and-vice-versa-in-c
  for (int j = 0; j < mvals; j++){
    std::ostringstream ostr; 
    ostr << j + 1;
    headerfile +=  " Var"  + ostr.str();
  }
  
  myfile << headerfile + "\n";
  
  //Loop to assign initial valies to Xvals and nsimRange
  //----------------------------------------------------------------
  
  
  for (int j = 0; j < nsim; j++){
    
    //Value to Nsimrange
    nsimRange(j) = j + 1;
    
    //Value to Xvals and valsRange
    for (int k = 0; k < mvals; k++){
      
      //Create xvals
      Xvals(j,k) = X(k); 
      
    }
  }
  
  
  //Assign simulation number that won't change in loop
  Printable(_,0) = nsimRange; //Simulation number
  
  //Write initial data to simulation
  Printable(_,1) = indicator; //Iteration number
  Printable(_,2) = t;         //Time
  
  for (int j = 0; j < mvals; j++){
    Printable(_,j + 3) = Xvals(_,j);  //Result from simulation
  }
  
  //Print to file
  myfile << Printable;
  

  //SIMULATION LOOP: NEEDS TO BE FASTER
  //----------------------------------------------------------------
  //Start simulation loop
  while (tf < tmax){
    
    //Check if additional options are required
    if (optscheck){
      
      //Check current time and print
      if (printtime){
        showTime(tf);
      }
      
      //Check maxiter
      if (indicator(1) > maxiter){
        warnMaxiter(maxiter);
        break;
      }
      
      //Check maxtime
      current       = time(0);
      if(current - loopinit > maxtime){
        warnTime( current - loopinit );
        break;
      }
      
    }
    
    //Update .time counter
    tf     = *std::min_element(t.begin(), t.end());

    //Get the ".a" parameter
    a      = pfun(t, Xvals, params); //pfun is the R matrix-valued function
    
    //Calculate .a0
    a0vals = colSumsinv(a, nsim, lfun);
    a0     = a0vals.first;
    
    //Break if some value gets negative
    if (a0vals.second){
      warnNegative(indicator(1));
      break;
    }
    
    //Simulate time change
    tau    = a0*log(1/runif(nsim));
  
    //Update time
    t      =  t + tau;
    
    //Choose where change occured
    matx   = MatVecHadamard(a, a0, nsim, lfun);
    
    //Random number vector to check when change occurs
    vec    = runif(nsim);
    
    //Get the columns that changed
    col    = getMatx(matx, v, vec, lfun, nsim, mvals);

    //Update X matrix
    Xvals += col;
    
    //To make sure we don't print all results (files might get too memory-intensive) 
    if ( ((indicator(1)-1) % kthsave) == 0){

      Printable(_,1) = indicator; //Iteration number
      Printable(_,2) = t;         //Time
      
      for (int j = 0; j < mvals; j++){
        Printable(_,j + 3) = Xvals(_,j);  //Result from simulation
      }
      
      //Print to file
      myfile << Printable;
      
    }
    
    //Update loop variable
    indicator = indicator + 1;
    
  }
  
  //Close file
  myfile.close();
  
  //Report finished
  Rcerr << "*******\n" << "FINISHED\n"  << "*******\n"<< std::endl;
  
}


