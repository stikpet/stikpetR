% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thumb_cohen_h.R
\name{th_cohen_h}
\alias{th_cohen_h}
\title{Rule-of-Thumb for Cohen h}
\usage{
th_cohen_h(h, qual = "cohen")
}
\arguments{
\item{h}{the Cohen h value}

\item{qual}{optional setting for which rule of thumb to use. Currently only 'cohen'}
}
\value{
A dataframe with:
\item{classification}{the qualification of the effect size}
\item{reference}{a reference for the rule of thumb used}
}
\description{
Simple function to use a rule-of-thumb for the Cohen h effect size.

This function is shown in this \href{https://youtu.be/sGfFB7Zzeas}{YouTube video} and the effect size is also described at \href{https://peterstatistics.com/Terms/EffectSizes/CohenH.html}{PeterStatistics.com}
}
\details{
Cohen's rule of thumb for Cohen g (1988, p. 198):\tabular{ll}{
   |h| \tab Interpretation \cr
   0.00 < 0.20 \tab negligible \cr
   0.20 < 0.50 \tab small \cr
   0.50 < 0.80 \tab medium \cr
   0.80 or more \tab large \cr
}


Note that Cohen actually just lists small = 0.20, medium = 0.50, and large = 0.80.
}
\section{Before, After and Alternatives}{

\code{\link{es_cohen_h}}, to determine Cohen h
\code{\link{es_cohen_h_os}}, to determine Cohen h', then use \code{\link{es_convert}}
}

\examples{
es = 0.6
th_cohen_d(es)

}
\references{
Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
