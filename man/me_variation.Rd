% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meas_variation.R
\name{me_variation}
\alias{me_variation}
\title{Measures of Quantitative Variation}
\usage{
me_variation(
  data,
  levels = NULL,
  measure = "std",
  ddof = 1,
  center = "mean",
  azs = "square"
)
}
\arguments{
\item{data}{: list or dataframe}

\item{levels}{: dictionary, optional coding to use}

\item{measure}{: {"std", "var", "mad", "madmed", "medad", "stddm", "cv", "cd", "own"}, optional the measure to determine. Default is "std"}

\item{ddof}{: float, optional option to adjust the division in standard deviation or variance with. Default is 1.}

\item{center}{: {"mean", "median", "mode"} or float, optional if measure is "own" the value to use as center. Default is "mean"}

\item{azs}{: {"square", "abs"}, optional if measure is "own" the way to avoid a zero sum. Either by squaring or absolute value}
}
\value{
A dataframe with:
\itemize{
\item \emph{value}, the value of the measure
\item \emph{measure}, description of the measure
}
}
\description{
Probably the most famous measure of dispersion is the standard deviation, but there are more. This function provides a variety of measures and allows the creation of your own version.
}
\details{
\strong{Standard Deviation} (std)

The formula used is:
\deqn{s = \sqrt{\frac{\sum_{i=1}^n \left(x_i - \bar{x}\right)^2}{n - d}}}

Where \eqn{d} is the offset specified at \emph{ddof}. By default this is 1, giving the sample standard deviation.

\strong{Variance} (var)

The formula used is:

\deqn{s^2 = \frac{\sum_{i=1}^n \left(x_i - \bar{x}\right)^2}{n - d}}

Where \eqn{d} is the offset specified at \emph{ddof}. By default this is 1, giving the sample standard deviation.

\strong{Mean Absolute Deviation} (mad)

The formula used is:
\deqn{MAD = \frac{\sum_{i=1}^n \left| x_i - \bar{x}\right|}{n}}

\strong{Mean Absolute Deviation from the Median} (madmed)

The formula used is:

\deqn{MAD = \frac{\sum_{i=1}^n \left| x_i - \tilde{x}\right|}{n}}

Where \eqn{\tilde{x}} is the median

\strong{Median Absolute Deviation} (medad)

The formula used is:
\deqn{MAD = MED\left(\left| x_i - \tilde{x}\right|\right)}

\strong{Decile Standard Deviation}

The formula used is (Siraj-Ud-Doulah, 2018, p. 310):
\deqn{s_{dm} = \sqrt{\frac{\sum_{i=1}^n \left(x_i - DM\right)^2}{n - d}}}

Where DM is the decile mean.

\strong{Coefficient of Variation} (cv)

The formula used is (Pearson, 1896, p. 277):
\deqn{CV = \frac{s}{\bar{x}}}

\strong{Coefficient of Diversity} (cd)

The formula used is (Siraj-Ud-Doulah, 2018, p. 310):
\deqn{CD = \frac{s_{dm}}{DM}}

\strong{Own}
it's possible to create one's own method. Decide on a specific center. Default options are the mean, median and mode. Then on either to sum the squared deviations or the absolute differences.
}
\references{
Pearson, K. (1896). Contributions to the mathematical theory of evolution. III. Regression, Heredity, and Panmixia. \emph{Philosophical Transactions of the Royal Society of London}. (A.), 1896, 253–318.

Siraj-Ud-Doulah, M. (2018). Alternative measures of standard deviation coefficient of variation and standard error. \emph{International Journal of Statistics and Applications, 8}(6), 309–315. https://doi.org/10.5923/j.statistics.20180806.04
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
