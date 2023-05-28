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
}
\section{Alternatives}{

It is also possible to do this with the library \emph{ggplot2}.
A video on how to use ggplot2 for this can be found \href{https://youtu.be/EhK9K9jwsco}{here}
}

\examples{
ordData <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
vi_bar_stacked_single(ordData)

}
\seealso{
An alternative chart for a single ordinal variable could be a single stacked bar-chart,
see \code{\link{vi_bar_stacked_single}}
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}
}