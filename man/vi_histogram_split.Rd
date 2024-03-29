% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_histogram_split.R
\name{vi_histogram_split}
\alias{vi_histogram_split}
\title{Split Histogram}
\usage{
vi_histogram_split(catField, scaleField, categories = NULL, ...)
}
\arguments{
\item{catField}{list or dataframe with the categories}

\item{scaleField}{list or dataframe with the scores}

\item{categories}{optional list with categories to use}

\item{...}{other parameters for use in geom_histogram function}
}
\value{
The split histogram
}
\description{
Based on a categorical field the scores for each category are plotted in a separate histogram and each of the histograms is placed underneath each other.

See \strong{vi_histogram()} for more details on histograms.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
