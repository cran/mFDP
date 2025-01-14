\name{mFDP.direc}
\alias{mFDP.direc}
\title{mFDP control for directional testing}

\description{Controls the median of the FDP, given a vector of test statistics.}

\usage{mFDP.direc(Ts, delta=0, gamma=0.05)}

\arguments{
  \item{Ts}{A vector containing the test statistics T_1,...,T_m}
  \item{delta}{Determines the null hypotheses H_1,...,H_m. 
  For every j, H_j is the hypothesis that T_j has mean at most delta.}
  \item{gamma}{Desired upper bound for the mFDP}
}

\value{The indices of the rejected hypotheses}

\examples{
set.seed(123)

tstats = rnorm(n=200)+1.5   #make some test statistics

mFDP.direc(Ts = tstats)  #indices of rejected hypotheses

}

