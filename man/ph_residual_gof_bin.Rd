% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poho_residual_gof_bin.R
\name{ph_residual_gof_bin}
\alias{ph_residual_gof_bin}
\title{Post-Hoc Residuals Using Binary Tests for GoF}
\usage{
ph_residual_gof_bin(
  data,
  test = "std-residual",
  expCount = NULL,
  mtc = "bonferroni",
  ...
)
}
\arguments{
\item{data}{dataframe with scores}

\item{test}{{"adj-residual", "std-residual", "binomial", "wald", "score"} optional test to use}

\item{expCount}{optional dataframe with categories and expected counts}

\item{mtc}{optional string. Any of the methods available in p_adjust() to correct for multiple tests}

\item{...}{optional additional parameters to be passed to the test}
}
\value{
a dataframe with:
\item{category}{the label of the first category}
\item{obs. count}{the observed count}
\item{exp. count}{the expected count}
\item{statistic}{the test statistic}
\item{p-value}{the unadjusted significance}
\item{adj. p-value}{the adjusted significance}
\item{test}{description of the test used}
}
\description{
This function will perform a residuals post-hoc test for each of the categories in a nominal field. This could either be a z-test using the standardized residuals, the adjusted residuals, or any of the one-sample binary tests.

The unadjusted p-values and Bonferroni adjusted p-values are both determined.
}
\details{
The formula used is for the adjusted residual test:
\deqn{z = \frac{F_i - E_i}{\sqrt{E_i\times\left(1 - \frac{E_i}{n}\right)}}}
\deqn{sig = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}

The formula used for the standardized residual test:
\deqn{z = \frac{F_i - E_i}{\sqrt{E_i}}}
\deqn{sig = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}

With:
\itemize{
\item \eqn{F_i}, the observed count for category $i$
\item \eqn{E_i}, the expected count for category $i$
\item \eqn{\Phi\left(\dots\right)}, the cumulative distribution function of the standard normal distribution
}

If no expected counts are provide it is assumed they are all equal for each category, i.e. \eqn{E_i = \frac{n}{k}}

The Bonferroni adjustment is calculated using:
\deqn{p_{adj} = \min \left(p \times k, 1\right)}

The other tests use the formula from the one-sample test variant, using the expected count/n as the expected proportion.

The adjusted residuals will gave the same result as using a one-sample score test. Some sources will also call these adjusted residuals as standardized residuals (Agresti, 2007, p. 38), and the standardized residuals used in this function as Pearson residuals (R, n.d.). Haberman (1973, p. 205) and Sharpe (2015, p. 3) are sources for the terminology used in this function.
}
\references{
Agresti, A. (2007). \emph{An introduction to categorical data analysis} (2nd ed.). Wiley-Interscience.
Haberman, S. J. (1973). The analysis of residuals in cross-classified tables. \emph{Biometrics, 29}(1), 205–220. doi:10.2307/2529686
R. (n.d.). Chisq.test \link{Computer software}. https://stat.ethz.ch/R-manual/R-devel/library/stats/html/chisq.test.html
Sharpe, D. (2015). Your chi-square test is statistically significant: Now what? Practical Assessment, \emph{Research & Evaluation, 20}(8), 1–10.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
