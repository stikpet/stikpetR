% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_fligner_policello.R
\name{ts_fligner_policello}
\alias{ts_fligner_policello}
\title{Fligner-Policello Test}
\usage{
ts_fligner_policello(dataVar, groupVar, ties = TRUE, cc = FALSE)
}
\arguments{
\item{dataVar}{A vector with the scores data}

\item{groupVar}{A vector with the group data}

\item{ties}{boolean to indicate the use of a ties correction}

\item{cc}{boolean to indicate the use of a continuity correction}
}
\value{
dataframe with the test statistic, p-value, and the test used
}
\description{
Fligner-Policello Test
}
\details{
The formula used is:
\deqn{z = \frac{N_Y - N_X}{2\times\sqrt{SS_X + SS_Y - M_X\times M_Y}}}
With:
\deqn{SS_X = \sum_{x\in X} \left(N_X - M_X\right)^2, SS_Y = \sum_{y\in Y} \left(N_Y - M_Y\right)^2}
\deqn{M_X = \frac{N_X}{n_x}, M_Y = \frac{N_Y}{n_y}}
\deqn{N_X = \sum_{x \in X} N\left(x\right), N_Y = \sum_{y \in Y} N\left(y\right)}
\deqn{N\left(y\right) = \sum_{x\in X} f\left(y, x\right)}
\deqn{N\left(x\right) = \sum_{y\in Y} f\left(x, y\right)}
\deqn{f\left(a, b\right) = \begin{cases} 1 & \text{ if } a> b \\ 0 & \text{ if } a\leq b \end{cases}}

In case of a tie correction (Hollander et al., 2014, p. 146):
\deqn{z = \frac{\left|N_Y - N_X\right| - 0.5}{2\times\sqrt{SS_X + SS_Y - M_X\times M_Y}}}
\deqn{f\left(a, b\right) = \begin{cases} 1 & \text{ if } a> b \\ 0.5 & \text{ if } a = b \\ 0 & \text{ if } a\leq b \end{cases}}

\emph{Symbols used:}
\itemize{
\item \eqn{X} the scores in the first category
\item \eqn{Y} the scores in the second category
\item \eqn{n_i} the number of scores in the i category
}

The test is described by Fligner and Policello (1981), and can also be found in Kloke and McKean (2015, p. 68)
}
\examples{
scores = c(5, 12, 3, 4, 6, 1, 11, 13, NA)
groups = c("A","A","A","B","B","B","B", NA, "C")
ts_fligner_policello(scores, groups)
ts_fligner_policello(scores, groups, ties=FALSE, cc=FALSE)
ts_fligner_policello(scores, groups, ties=FALSE, cc=TRUE)
ts_fligner_policello(scores, groups, ties=TRUE, cc=FALSE)
ts_fligner_policello(scores, groups, ties=TRUE, cc=TRUE)

}
\references{
Fligner, M. A., & Policello, G. E. (1981). Robust rank procedures for the Behrens-Fisher problem. \emph{Journal of the American Statistical Association, 76}(373), 162–168. https://doi.org/10.1080/01621459.1981.10477623

Hollander, M., Wolfe, D. A., & Chicken, E. (2014). \emph{Nonparametric statistical methods} (3rd ed.). John Wiley & Sons, Inc.

Kloke, J., & McKean, J. W. (2015). \emph{Nonparametric statistical methods using R}. CRC Press, Taylor & Francis.
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
