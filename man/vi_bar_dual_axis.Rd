% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_bar_dual_axis.R
\name{vi_bar_dual_axis}
\alias{vi_bar_dual_axis}
\title{Dual-Axis Bar Chart}
\usage{
vi_bar_dual_axis(data, varname = NULL)
}
\arguments{
\item{data}{the data from which to create a Pareto chart}

\item{varname}{a name for the data, if not provided the name of the data variable is used}
}
\value{
a chart in the plot window
}
\description{
A dual axis bar-chart is a bar-chart with two vertical axis. In this function it will
show both the count and cumulative proportion.

This chart could be used with a single ordinal variable.

This function is shown in this \href{https://youtu.be/J62B1zyTO3U}{YouTube video} and the visualisation is also described at \href{https://peterstatistics.com/Terms/Visualisations/bar-chart.html}{PeterStatistics.com}
}
\section{Before, After and Alternatives}{

Before the visualisation you might first want to get an impression using a frequency table:
\code{\link{tab_frequency}}, for a frequency table

After visualisation you might want some descriptive measures:
\code{\link{me_consensus}}, for the Consensus.
\code{\link{me_hodges_lehmann_os}}, for the Hodges-Lehmann Estimate (One-Sample).
\code{\link{me_median}}, for the Median.
\code{\link{me_quantiles}}, for Quantiles.
\code{\link{me_quartiles}}, for Quartiles / Hinges.
\code{\link{me_quartile_range}}, for Interquartile Range, Semi-Interquartile Range and Mid-Quartile Range.

or perform a test:
\code{\link{ts_sign_os}}, for One-Sample Sign Test.
\code{\link{ts_trinomial_os}}, for One-Sample Trinomial Test.
\code{\link{ts_wilcoxon_os}}, for One-Sample Wilcoxon Signed Rank Test.

Alternatives for this visualisation could be:
\code{\link{vi_bar_stacked_single}}, or Single Stacked Bar-Chart.
}

\examples{
#Example 1: dataframe
dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
ex1 = df1['mar1']
vi_bar_dual_axis(ex1);
vi_bar_dual_axis(ex1, varname="marital status");

}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
