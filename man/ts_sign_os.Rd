% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_sign_os.R
\name{ts_sign_os}
\alias{ts_sign_os}
\title{one-sample sign test}
\usage{
ts_sign_os(data, levels = NULL, mu = NULL)
}
\arguments{
\item{data}{A vector or dataframe}

\item{levels}{optional vector with levels in order}

\item{mu}{optional hypothesized median, otherwise the midrange will be used}
}
\value{
Dataframe with:
\item{mu}{the mean tested}
\item{p-value}{he significance (p-value)}
\item{test}{a description of the test used}
}
\description{
This function will perform one-sample sign test.

This function is shown in this \href{https://youtu.be/Hj15elBte44}{YouTube video} and the test is also described at \href{https://peterstatistics.com/Terms/Tests/SignOneSample.html}{PeterStatistics.com}
}
\details{
The test statistic is calculated using (Stewart, 1941, p. 236):
\deqn{p = 2\times B\left(n, \text{min}\left(n_+, n_-\right), \frac{1}{2}\right)}

\emph{Symbols used:}
\itemize{
\item \eqn{B\left(\dots\right)} is the binomial cumulative distribution function
\item \eqn{n} is the number of cases
\item \eqn{n_+} is the number of cases above the hypothesized median
\item \eqn{n_-} is the number of cases below the hypothesized median
\item \eqn{min} is the minimum value of the two values
}

The test is described in Stewart (1941), although there are earlier uses.

The paired version for example was already described by Arbuthnott (1710)
}
\section{Before, After and Alternatives}{

Before this measure you might want an impression using a frequency table or a visualisation:
\code{\link{tab_frequency}}, for a frequency table
\code{\link{vi_bar_stacked_single}}, or Single Stacked Bar-Chart.
\code{\link{vi_bar_dual_axis}}, for Dual-Axis Bar Chart.

After this you might want to determine an effect size measure:
\code{\link{es_common_language_os}}, for the Common Language Effect Size.
\code{\link{es_dominance}}, for the Dominance score.
\code{\link{r_rank_biserial_os}}, for the Rank-Biserial Correlation

Alternative tests:
\code{\link{ts_trinomial_os}}, for One-Sample Trinomial Test.
\code{\link{ts_wilcoxon_os}}, for One-Sample Wilcoxon Signed Rank Test.
}

\examples{
#Example 1: Text dataframe
file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
ex1 = df2[['Teach_Motivate']]
order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
ts_sign_os(ex1, levels=order)

#Example 2: Numeric data
ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
ts_sign_os(ex2) 

}
\references{
Arbuthnott, J. (1710). An argument for divine providence, taken from the constant regularity observ'd in the births of both sexes. \emph{Philosophical Transactions of the Royal Society of London, 27}(328), 186-190. https://doi.org/10.1098/rstl.1710.0011

Stewart, W. M. (1941). A note on the power of the sign test. \emph{The Annals of Mathematical Statistics, 12}(2), 236-239. https://doi.org/10.1214/aoms/1177731755
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
