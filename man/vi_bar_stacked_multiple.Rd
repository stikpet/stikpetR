% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_bar_stacked_multiple.R
\name{vi_bar_stacked_multiple}
\alias{vi_bar_stacked_multiple}
\title{Multiple Stacked Bar-Chart}
\usage{
vi_bar_stacked_multiple(catField, ordField, levels = NULL, ...)
}
\arguments{
\item{catField}{list or dataframe with the categories}

\item{ordField}{list or dataframe with the scores}

\item{levels}{optional list with the scores in order}

\item{...}{optional, other parameters for use in barplot function}
}
\value{
multiple stacked bar-chart
}
\description{
To visualise an ordinal variable, it often makes sense to stack the results. Stacking the results creates a compound bar chart, or sometimes stacked bar chart (Wilkinson, 2005, p. 157) or component bar chart (Zedeck, 2014, p. 54). It can be defined as: "a bar chart showing multiple bars stacked at each x-axis category, each representing a value of the stacking variable" (Upton & Cook, 2014, p. 88).

Instead of one bar (see \strong{vi_bar_stacked_single()}), we can create two or more (one for each group). This could then be considered a multiple compound bar-chart.
}
\details{
This function is more like a wrapper for the \strong{barplot()} from R \emph{graphics} library.
}
\examples{
file1 = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
df1 <- read.csv(file1, sep=",", na.strings=c("", "NA"))
vi_bar_stacked_multiple(df1[['mar1']], df1[['accntsci']], ylab= "percent", col=1:5)

cats = c(1, 1, 2, 2, 2, 3, 3, 3, 3)
scor = c(1, 2, 1, 1, 2, 1, 1, 1, 2)
vi_bar_stacked_multiple(cats, scor, ylab= "percent", col=1:5)


}
\references{
Upton, G., & Cook, I. (2014). \emph{Oxford: Dictionary of statistics} (3rd ed.). Oxford University Press.

Wilkinson, L. (2005). \emph{The grammar of graphics} (2nd ed). Springer.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
