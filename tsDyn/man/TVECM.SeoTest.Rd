\name{TVECM.SeoTest}
\alias{TVECM.SeoTest}
\title{No cointegration vs threshold cointegration test}
\description{
Test the null of no cointegration against threshold cointegration with bootstrap distribution of Seo (2006)
}
\usage{
TVECM.SeoTest(data,lag, beta, trim=0.1,nboot, plot=FALSE, hpc=c("none", "foreach"),check=FALSE)
}
\value{A list cointaining diverse informations:

Estimated threshold parameters and usual slope parameters.

Value of the test.

Critical and Pvalue from bootstrap distribution.
}
\arguments{
\item{data}{ time series }
\item{lag}{Number of lags to include in each regime}
\item{beta}{Pre-specified cointegarting value}
\item{trim}{trimming parameter indicating the minimal percentage of observations in each regime}
\item{nboot}{Number of bootstrap replications}
\item{plot}{Whether a grid with the SSR of each threshold should be
  printed}
\item{hpc}{Possibility to run the bootstrap on parallel core. See
  details in \code{\link{TVECM.HStest}}}
\item{check}{Possiblity to check the function by no sampling: the test value should be the sme as in the original data}
}
\details{
For this test, the cointegrating value has to be specified by the user. 

The model used is one where the threshold effect concerns only the cointegrating vector, and only in the outer regimes. 

Due to the presence of parameters unidentified under the null hypothesis, the test employed is a Sup-Wald test, that means that for each combination of the thresholds, a Wald Test is computed and the supremum of all tests is taken. 
For each bootstrap replication, this approach is taken, so that the test is really slow. 
}
\seealso{
\code{\link{TVECM}} for estimating a TVECM, \code{\link{TVECM.sim}} for simulating/bootstrap a TVECM, 
}

\author{Matthieu Stigler}
\keyword{ ts }
\references{
Seo, Myunghwan, 2006. "Bootstrap testing for the null of no cointegration in a threshold vector error correction model," Journal of Econometrics, vol. 127(1), pages 129-150, September.
}
\examples{
# As the function takes long long time to be executed, we show in in don't run environement
\dontrun{
data(zeroyld)

#can be useful to check whether the bootstrap is working: sithout sampling, results of boot should be same as original
#this is indeed not always the case duye to floating point algorithm
TVECM.SeoTest(zeroyld,lag=2, beta=1, trim=0.1,nboot=2, plot=FALSE,check=TRUE)

#then run the function:
TVECM.SeoTest(zeroyld,lag=2, beta=1, trim=0.1,nboot=100, plot=FALSE,check=FALSE)
}
}


