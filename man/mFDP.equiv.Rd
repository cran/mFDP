\name{mFDP.equiv}
\alias{mFDP.equiv}
\title{mFDP control for equivalence testing}

\description{Controls the median of the FDP, given a vector of test statistics.}

\usage{mFDP.equiv(Ts, delta, gamma=0.05)}

\arguments{
  \item{Ts}{A vector containing the test statistics T_1,...,T_m}
  \item{delta}{Determines the null hypotheses H_1,...,H_m. Should be >0.
  For every j, H_j is the null hypothesis that |T_j| has mean at least delta.
  The alternative hypothesis is that |T_j| has mean smaller than delta}
  \item{gamma}{Desired upper bound for the mFDP}
}

\value{The indices of the rejected hypotheses}

\examples{
set.seed(123)

tstats = rnorm(n=200)  #make some test statistics

mFDP.equiv(Ts = tstats,delta=2)  #indices of rejected hypotheses

}

