#' @rdname LQG_update_transitions
#' @title Update Markov transition kernels with new policy 
#' @description Function that takes the parameters produced by the policy recursion (LQG_policy_recursion)
#' and the parameters of the previous Gaussian transition kernels, and outputs the parameters
#' for the updated Gaussian transitions
#' @param policy_params list with keys:
#' \code{As}
#' \code{bs}
#' \code{cs}
#' @param transitions list with keys: 
#' \code{Kmats} array of size d x d x T 
#' \code{rvecs} array of size d x T
#' \code{Smats} array of size d x d x T
#' @return list with keys: 
#' \code{Kmats} array of size d x d x T 
#' \code{rvecs} array of size d x T
#' \code{Smats} array of size d x d x T
#' \code{inv.Smats} array of size d x d x T
#' @export

LQG_update_transitions = function(policy_params,transitions){
  
  # Extract the relevant parameters
  policy.bs = policy_params$bs
  policy.As = policy_params$As
  
  Kmats = transitions$Kmats
  rvecs = transitions$rvecs
  Smats = transitions$Smats
  inv.Smats = transitions$inv.Smats
  
  d = length(policy.bs[,1]) #dimension of state variable
  T = length(Kmats[1,1,]) #number of Markov transition kernels
  
  # update the transition parameters
  updated.inv.Smats = 2*policy.As + inv.Smats
  updated.Smats = array(apply(updated.inv.Smats, MARGIN = 3, FUN = solve), c(d,d,T))
  updated.Kmats = array(sapply(1:T, function(t) updated.Smats[,,t]%*%inv.Smats[,,t]%*%Kmats[,,t]), c(d,d,T))
  updated.rvecs = array(sapply(1:T, function(t) updated.Smats[,,t]%*%(inv.Smats[,,t]%*%rvecs[,t] - policy.bs[,t])), c(d,T))
  
  return(list(Kmats = updated.Kmats, rvecs = updated.rvecs, Smats = updated.Smats, inv.Smats = updated.inv.Smats))
}
