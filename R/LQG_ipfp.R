#' @rdname LQG_ipfp
#' @title Update Markov transition kernels with new policy 
#' @description Function to perform IPFP iterations in the exact LQG setting.
#' Outputs the updated Markov transition kernels.
#' @param nsteps number of IPFP iterations
#' @param initial list with keys:
#' \code{mean}
#' \code{cov}
#' @param target list with keys: 
#' \code{mean}
#' \code{inv.cov}
#' @param transitions list with keys: 
#' \code{Kmats} array of size d x d x T 
#' \code{rvecs} array of size d x T
#' \code{Smats} array of size d x d x T
#' @return transitions is a list with keys: 
#' \code{Kmats} array of size d x d x T 
#' \code{rvecs} array of size d x T
#' \code{Smats} array of size d x d x T
#' \code{inv.Smats} array of size d x d x T
#' @export

LQG_ipfp = function(nsteps,initial,target,transitions){
  
  for(i in 1:nsteps){
    
    #Compute the marginals of the path measure
    path_marginals = LQG_path_marginals(initial,transitions)
    
    #Extract the terminal path marginal
    terminal = list()
    terminal$cov = path_marginals$covs[,,T]
    terminal$inv.cov = path_marginals$inv.covs[,,T]
    terminal$mean = path_marginals$means[,T]
    
    #Compute Radon-Nikodym derivative: target/terminal
    RNderiv = LQG_RNderiv(terminal,target)
    
    #Calculate updated policies in backward recursion
    policy_params = LQG_policy_recursion(RNderiv,transitions)
    
    #Update the transitions
    transitions = LQG_update_transitions(policy_params,transitions)
  }
  
  return(transitions)
  
}
