#' @rdname LQG_path_marginals
#' @title Obtain marginal distributions of path measure 
#' @description Function to compute marginal distributions of path measure
#' given initial a Gaussian initial distribution and Gaussian Markov
#' transition kernels of the form: x_t | x_{t-1} is Gaussian with 
#' mean = K_t x_{t-1} + r_t and covariance = S_t.
#' "initial" is a list containing the initial parameters.
#' "transitions" is a list containing Kmats, rvecs and Smats, which are arrays containing 
#' the transition parameters. For each, the range of the last index gives the number of steps.
#' @param initial list with keys:
#' \code{mean}
#' \code{cov}
#' @param transitions list with keys: 
#' \code{Kmats} array of size d x d x T 
#' \code{rvecs} array of size d x T
#' \code{Smats} array of size d x d x T
#' @return list with keys: 
#' \code{mean}
#' \code{covs}
#' \code{inv.covs}
#' @export

LQG_path_marginals = function(initial,transitions){
  
  # Initialize
  mean.current = initial$mean
  cov.current = initial$cov
  
  # Extract transitions
  Kmats = transitions$Kmats
  rvecs = transitions$rvecs
  Smats = transitions$Smats
  
  # Parameters
  d = length(mean.current) #dimension of state variable
  T = length(Kmats[1,1,]) #number of Markov transition kernels
  
  # Allocate space for marginal means and covariances
  store.means = array(rep(0,d),c(d,T))
  store.covs = array(0,c(d,d,T))
  store.inv.covs = array(0,c(d,d,T))
  
  for(t in 1:T){
    # update
    mean.current = Kmats[,,t]%*%mean.current + rvecs[,t]
    cov.current = Smats[,,t] + Kmats[,,t]%*%cov.current%*%t(Kmats[,,t])
    
    # store
    store.means[,t] = mean.current
    store.covs[,,t] = cov.current
    store.inv.covs[,,t] = solve(cov.current)
  }
  return(list(means = store.means, covs = store.covs, inv.covs = store.inv.covs))
}
