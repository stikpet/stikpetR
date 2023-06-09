% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eff_size_glass_delta.R
\name{es_glass_delta}
\alias{es_glass_delta}
\title{Glass Delta}
\usage{
es_glass_delta(scores, groups, control = NULL)
}
\arguments{
\item{scores}{A vector with the scores data}

\item{groups}{A vector with the group data}

\item{control}{Optional to indicate which group is the 'control' group, otherwise the 1st found will be used.}
}
\value{
Hedges g value
}
\description{
Glass Delta
}
\details{
The formula used is (Glass, 1976, p. 7):
\deqn{\delta = \frac{\bar{x}_1 - \bar{x}_2}{s_2}}

With:
\deqn{s_2 = \sqrt{\frac{\sum_{i=1}^{n_2} \left(x_{2,i} - \bar{x}_2\right)^2}{n_2 - 1}}}
\deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}

\emph{Symbols used:}
\itemize{
\item \eqn{x_{i,j}} the j-th score in category i
\item \eqn{n_i} the number of scores in category i
}

Glass actually uses a ‘control group’ and \eqn{s_2} is then the standard deviation of the control group.
}
\examples{
scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
es_glass_delta(scores, groups)
es_glass_delta(scores, groups, control="international")

}
\references{
Glass, G. V. (1976). Primary, secondary, and meta-analysis of research. \emph{Educational Researcher, 5}(10), 3–8. https://doi.org/10.3102/0013189X005010003
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
