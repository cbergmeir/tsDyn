\name{LSTAR}
\alias{LSTAR}
\alias{lstar}
\title{Logistic Smooth Transition AutoRegressive model}
\description{
Logistic Smooth Transition AutoRegressive model.
}
\synopsis{
lstar(x, m, d=1, steps=d, series, mL, mH, mTh, thDelay, 
          thVar, th, gamma, trace=TRUE, control=list()) 
}
\usage{
lstar(x, m, d=1, steps=d, series, mL, mH, thDelay, 
          th, gamma, trace=TRUE, control=list())

lstar(series, m, d, steps, mL, mH, mTh,
    th, gamma, trace=TRUE, control=list())

lstar(series, m, d, steps, mL=m, mH=m, thVar,
    th, gamma, trace=TRUE, control=list())
}
\value{
  An object of class \code{nlar}, subclass \code{lstar}, i.e. a list
  with fitted model informations.
%  Among others:
%  \describe{
%    \item{mL,mH}{low and high regimes autoregressive orders}
%    \item{externThVar}{is the threshold variable external?}
%    \item{thVar}{threshold variable values}
%  }
%  If the thresholding variable isn't external, there is an additional
%  component \code{mTh}, containing coefficients for lagged time series
%  used for obtaining the treshold variable.
}
\arguments{
\item{x}{ time series }
\item{m, d, steps}{ embedding dimension, time delay, forecasting steps }
\item{series}{ time series name (optional) }
\item{mL}{autoregressive order for 'low' regime (dafult: m). Must be <=m}
\item{mH}{autoregressive order for 'high' regime (default: m). Must be
  <=m}
\item{thDelay}{'time delay' for the threshold variable (as multiple of
  embedding time delay d)}
\item{mTh}{coefficients for the lagged time series, to obtain the
  threshold variable}
\item{thVar}{external threshold variable}
\item{th, gamma}{starting values for coefficients in the LSTAR model. If missing, a grid search is performed}
\item{trace}{should additional infos be printed? (logical)}
\item{control}{further arguments to be passed as \code{control} list to
  \code{\link{optim}}}
}
\details{
  \deqn{ x_{t+s} = ( \phi_{1,0} + \phi_{1,1} x_t + \phi_{1,2} x_{t-d} + \dots +
  \phi_{1,mL} x_{t - (mL-1)d} ) G( z_t, th, \gamma ) +
  ( \phi_{2,0} + \phi_{2,1} x_t + \phi_{2,2} x_{t-d} + \dots + \phi_{2,mH}
  x_{t - (mH-1)d} ) (1 - G( z_t, th, \gamma ) ) + \epsilon_{t+steps}}{
  x[t+steps] = ( phi1[0] + phi1[1] x[t] + phi1[2] x[t-d] + \dots +
  phi1[mL] x[t - (mL-1)d] ) G( z[t], th, gamma )
  + ( phi2[0] + phi2[1] x[t] + phi2[2] x[t-d] + \dots + phi2[mH] x[t -
  (mH-1)d] ) (1 - G( z[t], th, gamma ) ) + eps[t+steps]
  }
with \var{z} the treshold variable, and \eqn{G} the logistic function,
computed as \code{plogis(q, location = th, scale = 1/gamma)}, so see
\code{\link{plogis}} documentation for details on the logistic function
formulation and parameters meanings.  
The threshold variable can alternatively be specified by:
\describe{
\item{mTh}{ \eqn{z[t] = x[t] mTh[1] + x[t-d] mTh[2] + \dots +
    x[t-(m-1)d] mTh[m]} } 
\item{thDelay}{ \eqn{z[t] = x[t - thDelay*d ]} }
\item{thVar}{ \eqn{z[t] = thVar[t]} }
}

Note that if starting values for phi1 and phi2 are provided, isn't
necessary to specify mL and mH. Further, the user has to specify only
one parameter between mTh, thDelay and thVar for indicating the
threshold variable. 

Estimation is done by analytically determining \var{phi1} and \var{phi2}
(through linear regression) and then minimizing residuals sum of squares
with respect to \var{th} and \var{gamma}. These two steps are repeated
until convergence is achieved. For the nonlinear estimation of the
parameters \var{th} and \var{gamma}, the program uses the
\code{\link{optim}} function, with its default optimization method. You 
can pass further arguments directly to the 'control' list argument of
this function. For istance, the option \code{maxit} maybe useful when
there are convergence issues (see examples).
}

\seealso{
\code{\link{plot.lstar}} for details on plots produced for this model
from the \code{plot} generic. 
}
\author{ Antonio, Fabio Di Narzo }
\examples{
#fit a LSTAR model. Note 'maxit': slow convergence
mod.lstar <- lstar(log10(lynx), m=2, mTh=c(0,1), control=list(maxit=3000))
mod.lstar
}
\keyword{ ts }
\references{
Non-linear time series models in empirical finance, Philip Hans Franses
and Dick van Dijk, Cambridge: Cambridge University Press (2000). 

Non-Linear Time Series: A Dynamical Systems Approach, Tong, H., Oxford:
Oxford University Press (1990). 
}
