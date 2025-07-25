% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_student_t_ps.R
\name{ts_student_t_ps}
\alias{ts_student_t_ps}
\title{Student t Test (Paired Samples)}
\usage{
ts_student_t_ps(field1, field2, dmu = 0)
}
\arguments{
\item{field1}{the scores on the first variable}

\item{field2}{the scores on the second variable}

\item{dmu}{difference according to null hypothesis (default is 0)}
}
\value{
A dataframe with:
\item{n}{the number of scores}
\item{statistic}{the test statistic (t-value)}
\item{df}{the degrees of freedom}
\item{pValue}{the significance (p-value)}
}
\description{
The assumption about the population (null hypothesis) for this test is a pre-defined difference between two means, usually zero (i.e. the difference between the (arithmetic) means is zero, they are the same in the population). If the p-value (significance) is then below a pre-defined threhold (usually 0.05), the assumption is rejected.
}
\details{
The formula used is:
\deqn{t_p = \frac{\bar{d} - d_{H0}}{SE}}
\deqn{sig. = 2\times\left(1 - \text{T}\left(\left|t_p\right|, df\right)\right)}
With:
\deqn{\bar{d} = \bar{x}_1 - \bar{x}_2}
\deqn{SE = \sqrt{\frac{\sigma_{s^2}}{n}}}
\deqn{s_d^2 = \frac{\sum_{i=1}^n \left(d_i - \bar{d}_i\right)^2}{n - 1}}
\deqn{d_i = x_{i,1} - x_{i,2}}
\deqn{\bar{d} = \frac{\sum_{i=1}^n d_i}{n}}

\strong{Symbols used:}
\itemize{
\item \eqn{n} the number of pairs (sample size)
\item \eqn{x_{i,1}} the i-th score of the first variable
\item \eqn{x_{i,2}} the i-th score of the second variable
\item \eqn{d_{H0}} the expected difference in the population
\item \eqn{\text{T}\left(\dots, \dots\right)} the cumulative distribution function of the Student t distribution
}

\strong{Alternatives}

R's \emph{stats} library

t.test(var1, var2, paired=TRUE)

t.test(var1, var2, paired=TRUE, mu=5)
}
\references{
Student. (1908). The probable error of a mean. \emph{Biometrika, 6}(1), 1-25. doi:10.1093/biomet/6.1.1
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
