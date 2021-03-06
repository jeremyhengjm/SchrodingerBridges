% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LQG_path_marginals.R
\name{LQG_path_marginals}
\alias{LQG_path_marginals}
\title{Obtain marginal distributions of path measure}
\usage{
LQG_path_marginals(initial, transitions)
}
\arguments{
\item{initial}{list with keys:
\code{mean}
\code{cov}}

\item{transitions}{list with keys: 
\code{Kmats} array of size d x d x T 
\code{rvecs} array of size d x T
\code{Smats} array of size d x d x T}
}
\value{
list with keys: 
\code{mean}
\code{covs}
\code{inv.covs}
}
\description{
# Function to compute marginal distributions of path measure
given initial a Gaussian initial distribution and Gaussian Markov
transition kernels of the form: x_t | x_{t-1} is Gaussian with 
mean = K_t x_{t-1} + r_t and covariance = S_t.
"initial" is a list containing the initial parameters.
"transitions" is a list containing Kmats, rvecs and Smats, which are arrays containing 
the transition parameters. For each, the range of the last index gives the number of steps.
}
