% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LQG_ipfp.R
\name{LQG_ipfp}
\alias{LQG_ipfp}
\title{Update Markov transition kernels with new policy}
\usage{
LQG_ipfp(nsteps, initial, target, transitions)
}
\arguments{
\item{nsteps}{number of IPFP iterations}

\item{initial}{list with keys:
\code{mean}
\code{cov}}

\item{target}{list with keys: 
\code{mean}
\code{inv.cov}}

\item{transitions}{list with keys: 
\code{Kmats} array of size d x d x T 
\code{rvecs} array of size d x T
\code{Smats} array of size d x d x T}
}
\value{
transitions is a list with keys: 
\code{Kmats} array of size d x d x T 
\code{rvecs} array of size d x T
\code{Smats} array of size d x d x T
\code{inv.Smats} array of size d x d x T
}
\description{
# Function to perform IPFP iterations in the exact LQG setting.
Outputs the updated Markov transition kernels.
}
