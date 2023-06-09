% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_cochran_q.R
\name{ts_cochran_q}
\alias{ts_cochran_q}
\title{Cochran Q Test}
\usage{
ts_cochran_q(data, success = 1)
}
\arguments{
\item{data}{dataframe with the binary scores}

\item{success}{indicator for what is considered a success (default is 1)}
}
\value{
A dataframe with:
\item{statistic}{the test statistic}
\item{df}{the degrees of freedom}
\item{pValue}{the significance (p-value)}
}
\description{
Cochran Q Test
}
\details{
The formula used (Cochran, 1950, p. 259):
\deqn{Q = \frac{\left(k - 1\right)\times\sum_{j=1}^k\left(C_j - \bar{C}\right)^2}{k\times\sum_{i=1}^n R_i - \sum_{i=1}^n R_i^2}}
\deqn{sig. = 1 - \chi^2\left(Q, df\right)}
With:
\deqn{df = k - 1}

\emph{Symbols used:}
\itemize{
\item \eqn{C_j} the number of successes in category j
\item \eqn{k} the number of categories (factors)
\item \eqn{R_i} the number of succeses in case i
\item \eqn{n} the number of cases
}

\strong{Alternatives}

\emph{library(nonpar)}

matr = cbind(var1, var2, var3, var4)

cochrans.q(matr)

\emph{library(RVAideMemoire)}

myData.long<-reshape(dFr, varying=c("var1", "var2", "var3", "var4"), v.names="score", timevar="var", times=c("var1", "var2", "var3", "var4"),new.row.names = 1:1000, direction="long")

cochran.qtest(score~var |id, data=myData.long)
}
\examples{
var1 = c(0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1)
var2 = c(0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)
var3 = c(0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0)
var4 = c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1)
dFr = data.frame(var1, var2, var3, var4)
ts_cochran_q(dFr)

}
\references{
Cochran, W. G. (1950). The comparison of percentages in matched samples. \emph{Biometrika, 37}(3/4), 256–266. https://doi.org/10.2307/2332378
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
