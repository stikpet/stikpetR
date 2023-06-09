% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_fisher_owa.R
\name{ts_fisher_owa}
\alias{ts_fisher_owa}
\title{Fischer/Classic One-Way ANOVA}
\usage{
ts_fisher_owa(scores, groups)
}
\arguments{
\item{scores}{the numeric scores variable}

\item{groups}{the groups variable}
}
\value{
A dataframe with:
\item{statistic}{the F-statistic from the test}
\item{df1}{the first degrees of freedom}
\item{df2}{the second degrees of freedom}
\item{pValue}{the significance (p-value)}
}
\description{
Fischer/Classic One-Way ANOVA
}
\details{
The formula used:
\deqn{F_F = \frac{df_w\times SS_b}{df_b\times SSw}}
\deqn{df_b = k - 1}
\deqn{df_w = n - k}
\deqn{sig. = 1 - F\left(F_F, df_b, df_w\right)}
With:
\deqn{SS_b = \sum_{j=1}^k n_j\times\left(\bar{x}_j - \bar{x}\right)^2}
\deqn{SS_w = SS_t - SS_b}
\deqn{SS_t = \sum_{j=1}^k \sum_{i=1}^{n_j}\left(x_{i,j} - \bar{x}\right)^2}
\deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}
\deqn{\bar{x} = \frac{\sum_{j=1}^k n_j\times\bar{x}_j }{n} = \frac{\sum_{j=1}^k \sum_{i=1}^{n_j} x_{i,j}}{n}}
\deqn{n = \sum_{j=1}^k n_j}

Alternative format of the F-statistic equation (but the same result):
\deqn{F_F = \frac{MS_b}{MSw}}

\emph{Symbols}
\itemize{
\item \eqn{x_{i,j}} the i-th score in category j
\item \eqn{n} the total sample size
\item \eqn{n_j} the number of scores in category j
\item \eqn{k} the number of categories
\item \eqn{\bar{x}_j} the mean of the scores in category j
\item \eqn{SS_i} the sum of squares of i (sum of squared deviation of the mean)
\item \eqn{df_i} the degrees of freedom of i
\item \eqn{b} is between = factor = treatment = model
\item \eqn{w} is within = error (the variability within the groups)
}

Note that the Fisher-Pitman test (Pitman, 1937a, 1937b, 1938) uses a different approach but will lead to the same result.

I'm not fully sure what the original source is for the Fisher test, but likely either of his
sources from 1918, 1921, 1925 or 1935.
}
\examples{
scores = c(20, 50, 80, 15, 40, 85, 30, 45, 70, 60, NA, 90, 25, 40, 70, 65, NA, 70, 98, 40, 65, 60, 35, NA, 50, 40, 75, NA, 65, 70, NA, 20, 80, 35, NA, 68, 70, 60, 70, NA, 80, 98, 10, 40, 63, 75, 80, 40, 90, 100, 33, 36, 65, 78, 50)
groups = c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
ts_fisher_owa(scores, groups)

}
\references{
Fisher, R. A. (1918). The correlation between relatives on the supposition of mendelian inheritance. \emph{Transactions of the Royal Society of Edinburgh, 52}, 399–433.

Fisher, R. A. (1921). On the “probable error” of a coefficient of correlation deduced from a small sample. \emph{Metron, 1}, 3–32.

Fisher, R. A. (1925). \emph{Statistical methods for research workers}. Oliver and Boyd.

Fisher, R. A. (1935). \emph{The design of experiments}. Oliver and Boyd.

Pitman, E. J. G. (1937a). Significance tests which may be applied to samples from any populations. \emph{Supplement to the Journal of the Royal Statistical Society, 4}(1), 119–130. https://doi.org/10.2307/2984124

Pitman, E. J. G. (1937b). Significance tests which may be applied to samples from any populations. II. The correlation coefficient test. \emph{Supplement to the Journal of the Royal Statistical Society, 4}(2), 225–232. https://doi.org/10.2307/2983647

Pitman, E. J. G. (1938). Significance tests which may be applied to samples from any populations: III The analysis of variance test. \emph{Biometrika, 29}(3–4), 322–335. https://doi.org/10.1093/biomet/29.3-4.322
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
