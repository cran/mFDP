\name{get.bound}
\alias{get.bound}
\title{Compute a 50 percent confidence upper bound for the number of false positives}

\description{For a p-value rejection threshold \code{t}, compute the 50 percent confidence upper bound for the number of false positives. The bounds are simultaneous over \code{t}.}

\usage{get.bound(t,c,kappa.max)}

\arguments{
  \item{t}{The p-value threshold}
  \item{c}{The tuning parameter, which influences the intercept and slope of the enveloppe. Should be numeric.}
  \item{kappa.max}{This value needs to be computes based on the p-values. Together with \code{c} it defines the bounds.}
}

\value{A non-negative integer, which is a median unbiased (or upward biased) estimate of the number false positives.}

\examples{
#Suppose the envelope that has been computed is defined by c=0.002 and kappa.max=0.001.
#We can then evaluate the envelope at several thresholds t as below.
#This is equivalent to simply entering the formula floor((t+c)/kappa.max).

#50 percent confidence upper bound for nr of false positives, if p-value threshold of 0.01 is used:
get.bound(t=0.01,c=0.002,kappa.max=0.001)   #12

#50 percent confidence upper bound for nr of false positives, if p-value threshold of 0.02 is used:
get.bound(t=0.02,c=0.002,kappa.max=0.001)   #22

#50 percent confidence upper bound for nr of false positives, if p-value threshold of 0.03 is used:
get.bound(t=0.03,c=0.002,kappa.max=0.001)   #32
}

