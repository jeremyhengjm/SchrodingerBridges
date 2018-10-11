#' @rdname LQG_policy_recursion
#' @title Compute backward recursion to obtain the optimal policy
#' @description Compute backward recursion to obtain the optimal policy, which has
#' the form at time t -log(psi_t) = x_t' A_t x_t + x_t'b_t + c_t.
#' @param RNderiv list with keys:
#' \code{A} 
#' \code{b}
#' \code{c}
#' @param transitions list with keys: 
#' \code{Kmats} array of size d x d x T 
#' \code{rvecs} array of size d x T
#' \code{Smats} array of size d x d x T
#' @return list with keys: 
#' \code{As} 
#' \code{bs}
#' \code{cs}
#' @export

LQG_policy_recursion = function(RNderiv,transitions){
  
  # Initialize
  A.current = RNderiv$A
  b.current = RNderiv$b
  c.current = RNderiv$c
  
  # Extract transitions
  Kmats = transitions$Kmats
  rvecs = transitions$rvecs
  Smats = transitions$Smats
  inv.Smats = transitions$inv.Smats
  
  # Parameters
  d = length(b.current) #dimension of state variable
  T = length(Kmats[1,1,]) #number of Markov transition kernels
  
  # Allocate space for (A, b, c) defining psi*
  store.cs = rep(0,T)
  store.bs = array(0,c(d,T))
  store.As = array(0,c(d,d,T))
  
  # psi_T given by RN derivative
  store.cs[T] = c.current
  store.bs[,T] = b.current
  store.As[,,T] = A.current
  
  for(t in (T-1):1){
    # updates
    c.current = c.current + 0.5*t(rvecs[,t+1])%*%inv.Smats[,,t+1]%*%rvecs[,t+1] -
      0.25*t(inv.Smats[,,t+1]%*%rvecs[,t+1]-b.current)%*%solve(A.current+0.5*inv.Smats[,,t+1],inv.Smats[,,t+1]%*%rvecs[,t+1]-b.current)
    
    b.current = t(Kmats[,,t+1])%*%(inv.Smats[,,t+1]%*%rvecs[,t+1] - 0.5*inv.Smats[,,t+1]%*%solve(A.current+0.5*inv.Smats[,,t+1],inv.Smats[,,t+1]%*%rvecs[,t+1]-b.current))
    
    A.current = 0.5*t(Kmats[,,t+1])%*%(inv.Smats[,,t+1] - 0.5*inv.Smats[,,t+1]%*%solve(A.current+0.5*inv.Smats[,,t+1],inv.Smats[,,t+1]))%*%Kmats[,,t+1]
    
    # store
    store.cs[t] = c.current
    store.bs[,t] = b.current
    store.As[,,t] = A.current
  }
  
  return(list(cs = store.cs, bs = store.bs, As = store.As))
}
