% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/help_kendall_alg.R
\name{he_kendall}
\alias{he_kendall}
\title{Helper Function - Kendall Algorithm}
\usage{
he_kendall(n, c)
}
\arguments{
\item{n}{the sample size}

\item{c}{the number of concordant pairs}
}
\value{
\item{pValue}{upper tail p-value of Kendall tau Distribution}
}
\description{
Helper Function - Kendall Algorithm
}
\details{
An algorithm found at https://github.com/scipy/scipy/blob/v1.10.1/scipy/stats/_mstats_basic.py#L774-L898
was adapted. This refers to Kendall (1970), and uses the helper function \emph{he_kendall(n, C)}.
Where \eqn{C = n_c}, i.e. the number of concordant pairs.
This algorithm already returns a two-tailed result.
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
