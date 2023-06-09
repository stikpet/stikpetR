% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thumb_yule_q.R
\name{th_yule_q}
\alias{th_yule_q}
\title{Rules of thumb for Yule Q}
\usage{
th_yule_q(q, qual = "glen")
}
\arguments{
\item{q}{the Yule Q value}

\item{qual}{optional for which rule-of-thumb to use. Currently only "glen"}
}
\value{
the qualification/classification
}
\description{
Rules of thumb for Yule Q
}
\details{
Glen rule of thumb for Yule Q (2017):\tabular{ll}{
   |Q| \tab Interpretation \cr
   0.00 < 0.30 \tab negligible \cr
   0.30 < 0.50 \tab moderate \cr
   0.50 < 0.70 \tab substantial \cr
   0.70 or more \tab very strong \cr
}
}
\examples{
q = 0.4285714
th_yule_q(q)

}
\references{
Glen, S. (2017, August 16). Gamma Coefficient (Goodman and Kruskal’s Gamma) & Yule’s Q. Statistics How To. https://www.statisticshowto.com/gamma-coefficient-goodman-kruskal/
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
