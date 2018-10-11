#' @rdname LQG_RNderiv
#' @title Obtain Radon Nikodym derivative 
#' @description Obtain Radon Nikodym derivative between target and terminal distribution of current path measure. 
#' This has the form -log(RN) = x' A x + x'b + c.
#' @param terminal list with keys:
#' \code{mean}
#' \code{inv.cov}
#' @param target list with keys: 
#' \code{mean}
#' \code{inv.cov}
#' @return list with keys: 
#' \code{A} 
#' \code{b}
#' \code{c}
#' @export
################## LQG_RNderiv ##################
# Function to calculate the RN derivative based on the terminal path marginal
LQG_RNderiv = function(terminal,target){
  
  # Store the parameters in a list
  RNderiv = list()
  
  # A matrix
  RNderiv$A = 0.5*(target$inv.cov - terminal$inv.cov)
  
  # b vector
  RNderiv$b = -as.vector((target$inv.cov%*%target$mean - terminal$inv.cov%*%terminal$mean))
  
  # c constant
  RNderiv$c = 0.5*(t(target$mean)%*%target$inv.cov%*%target$mean 
                   - t(terminal$mean)%*%terminal$inv.cov%*%terminal$mean
                   + determinant(terminal$inv.cov, log=T)$modulus[1] 
                   - determinant(target$inv.cov, log=T)$modulus[1])
  
  return(RNderiv)
}
