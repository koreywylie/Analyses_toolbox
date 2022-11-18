# Distance Covariance fns.
#   Fast calculation of Distance Covariance, using C++ & Rcpp,
#    with Euclidean Distances calculation as detailed by:
#      http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
#    including decomposition of Euclidean distance as suggested by:
#      http://blog.smola.org/post/969195661/in-praise-of-the-second-binomial-formula
#  
#               d_ij^2 = ||a_i||^2 + ||b_j||^2 - 2 * <a_i, b_j>
#  
#      where a_i: column vector from matrix A
#            b_j: column vector from matrix B
#            d_ij:  Euclidean distance between column is from A & j from B
#            ||.||  : L2 norm
#            <a, b> : inner product of vectors a, b
#  
# Author: Korey Wylie, MD
# License: GPL (>= 2)
#
# Example usage:
# > source('dCov_fns.R')
# > x = matrix(rnorm(100), 5, 20)
# > y = matrix(rnorm(50), 5, 10)
# > dCov(x, y)
#

library(Rcpp)


cppFunction(depends='RcppArmadillo',
            '// Fast calculation of Distance Covariance, using C++ & Rcpp
            
            double distCovar(NumericMatrix Xr, NumericMatrix Yr) {
            int n = Xr.nrow(),
            j = Xr.ncol(),
            k = Yr.ncol();
            arma::mat X = arma::mat(Xr.begin(), n, j, false);
            arma::mat Y = arma::mat(Yr.begin(), n, k, false);
            
            arma::colvec Xn = sum(square(X),1);
            arma::mat A = -2 * (X * X.t());
            A.each_col() += Xn;
            A.each_row() += Xn.t();
            A.diag().zeros();
            A = sqrt(A);
            A.replace(arma::datum::nan, 0);
            
            double muA = mean(mean(A));
            arma::colvec muACols = mean(A,1);
            arma::rowvec muARows = mean(A,0);
            A.each_col() -= muACols;
            A.each_row() -= muARows;
            A += muA;
            
            arma::colvec Yn = sum(square(Y),1);
            arma::mat B = -2 * (Y * Y.t());
            B.each_col() += Yn;
            B.each_row() += Yn.t();
            B.diag().zeros();
            B = sqrt(B);
            B.replace(arma::datum::nan, 0);
            
            double muB = mean(mean(B));
            arma::colvec muBCols = mean(B,1);
            arma::rowvec muBRows = mean(B,0);
            B.each_col() -= muBCols;
            B.each_row() -= muBRows;
            B += muB;
            
            double dCov = n*mean(mean(A % B))/(muA*muB);
            return dCov;
            } ')



#################################################################
flatten_and_format <- function(x){
  #################################################################
  # Formats data for input into Distance Corr. fns.,
  #   typically by flattening arrays & lists
  #   & aligning vectors as columns.
  # Returns data as matrix,
  #   w/ subjs. in rows
  #
  
  if ((class(x)=='matrix') && (nrow(x) > 1)){  # do nothing if x is already a matrix or column vector
  x = as.matrix(x)              # ...pass, data already correctly formatted
  }else if ((class(x)=='matrix') && (nrow(x) == 1) && (ncol(x) > 1)){ # format row vector as column vector
    x = t(x)
  }else if ((class(x)=='array') && (length(dim(x)) == 1)){
    x = as.matrix(x)                # format 1d array as a column vector
  }else if (class(x)=='array'){ # assume subjs. in 1st, dim., concatenates all other dims
    x = t(apply(x,1,c))
  }else if (class(x)=="data.frame"){ # ensure data frames data are numeric
    x = sapply(x, as.numeric)
    x = as.matrix(x, dim(x)[1], dim(x)[2])
  }else if (class(x)=="list"){  # vectorize & flatten list of arrays into matrix
    if (is.vector(x[[1]])){  #check for agreement of array dims. w/n list
      for (sub in 2:length(x)){
        if (!all(length(x[[1]]) == length(x[[sub]]))){
          stop("dims. of vectors in input list must agree")
        }
      }
    }else{
      for (sub in 2:length(x)){
        if (!all(dim(x[[1]]) == dim(x[[sub]]))){
          stop("dims. of arrays in input list must agree")
        }
      }
    }
    x = matrix(sapply(x,c),length(x),prod(dim(x[[1]])), byrow=T)
  }else if (is.vector(x)){      # ensure x is a column vector, subjs. in rows
    x = matrix(x, length(x),1) 
  }else{
    stop("could not understand input var. when formatting")
  }
  return(x)
} #################################################################

#################################################################
dCov <- function(x, y) {
  #################################################################
  # Original Distance Covariance
  #
  # INPUTS:
  #   x, y : matrices or column vectors of obs. in rows
  # OUTPUTS : distance correlation
  # REQUIRES: C++ fns. below, created using Rcpp from 'dCor_fns.R'
  #   flatten_and_format()
  #   distCovar()
  #
  
  x = flatten_and_format(x)
  m = nrow(x)
  y = flatten_and_format(y)
  n = nrow(y)
  
  if (n != m) stop("sample sizes must agree")
  if (! (all(is.finite(x))))
    stop("data contains missing or infinite values")
  if (! (all(is.finite(y))))
    stop("data contains missing or infinite values")
  
  return(distCovar(x, y))
} #################################################################
