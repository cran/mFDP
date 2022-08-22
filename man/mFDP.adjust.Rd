\name{mFDP.adjust}
\alias{mFDP.adjust}
\title{compute mFDP-adjusted p-values}

\description{Provides mFDP-adjusted p-values, given a vector of p-values.}

\usage{mFDP.adjust(P, c="1/(2m)", s1=0, s2=0.1)}

\arguments{
  \item{P}{A vector of (raw, i.e. unadjusted) p-values.}
  \item{c}{The tuning parameter, which influences the intercept and slope of the envelope. Should either be numeric or \code{"1/(2m)"} or \code{"1/m"}.}
  \item{s1}{The smallest p-value threshold of interest.}
  \item{s2}{The largest p-value threshold of interest.}
}

\value{A vector of mFDP-adjusted p-values. Some can be infinity - which can be interpreted as 1.}

\examples{
set.seed(5193)

### make some p-values
m=500        #the nr of hypotheses
nrfalse=100  #the nr of false hypotheses

tstats = rnorm(n=m)   #m test statistics
tstats[1:nrfalse] = tstats[1:nrfalse] + 3    #add some signal
P = 1 - pnorm(tstats)                        #compute p-values

P.adjusted = mFDP.adjust(P=P)     #mFDP-adjusted p-values. Be careful with interpretation.

min(P.adjusted) #0.0208    
}

