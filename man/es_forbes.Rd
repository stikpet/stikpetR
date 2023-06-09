% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eff_size_forbes.R
\name{es_forbes}
\alias{es_forbes}
\title{Forbes Coefficient}
\usage{
es_forbes(var1, var2)
}
\arguments{
\item{var1}{A vector with the binary data from the first variable}

\item{var2}{A vector with the binary data from the second variable}
}
\value{
Forbes Coefficient
}
\description{
Forbes Coefficient
}
\details{
The formula used (Forbes, 1907, p. 279):
\deqn{F = \frac{n\times\min\left(a, d\right)}{C_1\times R_1}}

\emph{Symbols used:}
\itemize{
\item \eqn{a} the count in the top-left cell of the cross table
\item \eqn{b} the count in the top-right cell of the cross table
\item \eqn{c} the count in the bottom-left cell of the cross table
\item \eqn{d} the count in the bottom-right cell of the cross table
\item \eqn{R_1} the sum of counts in the 1st row
\item \eqn{C_1} the sum of counts in the 1st column
}

The coefficient has a value of 1 if there is no association, while it has a value of 0 or 2 when there is a perfect one.
To adjust to the more traditional range of -1 to 1, Cole 1 simply subtracts one from the Forbes coefficient.

Note that Alroy added an adjustment to Forbes F, see \emph{es_alroy_f()} for this adjusted version.
}
\examples{
bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
es_forbes(bin1, bin2)

}
\references{
Forbes, S. A. (1907). On the local distribution of certain Illinois fishes: An essay in statistical ecology. \emph{Illinois Natural History Survey Bulletin, 7}(8), 273–303.
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
