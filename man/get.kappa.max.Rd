\name{get.kappa.max}
\alias{get.kappa.max}
\title{Based on a vector of raw p-values, compute \code{kappa.max}, which  defines the mFDP envelope}

\description{Based on a vector of unadjusted(!) p-values, compute \code{kappa.max}, which together with \code{c} defines the mFDP envelope}

\usage{get.kappa.max(P, c="1/(2m)", s1=0, s2=0.1)}

\arguments{
  \item{P}{A vector of p-values.}
  \item{c}{The tuning parameter, which influences the intercept and slope of the enveloppe. Should either be numeric or \code{"1/(2m)"} or \code{"1/m"}.}
  \item{s1}{The smallest p-value threshold of interest. Non-negative.}
  \item{s2}{The largest p-value threshold of interest. Should be larger than \code{s1} and at most 1.}
}

\value{\code{kappa.max}, which together with \code{c} defines the mFDP envelope.}

\examples{
set.seed(5193)

### Make some p-values
m=500        #the nr of hypotheses
nrfalse=100  #the nr of false hypotheses

tstats = rnorm(n=m)   #m test statistics
tstats[1:nrfalse] = tstats[1:nrfalse] + 3    #add some signal
P = 1 - pnorm(tstats)                        #compute p-values

### Compute kappa.max. (Taking c to be the default value 1/(2m).)
kappa.max = get.kappa.max(P=P)
kappa.max
}

