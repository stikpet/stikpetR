% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eff_size_alt_ratio.R
\name{es_alt_ratio}
\alias{es_alt_ratio}
\title{Alternative Ratio}
\usage{
es_alt_ratio(data, codes = NULL, p0 = 0.5, category = NULL)
}
\arguments{
\item{data}{vector with the data}

\item{codes}{optional vector with the two codes to use}

\item{p0}{optional the hypothesized proportion for the first category (default is 0.5)}

\item{category}{optional category to label as 'success', otherwise the first category found is used.}
}
\value{
Dataframe with:
\item{AR1}{the alternative category for one category}
\item{AR2}{the alternative category for the other category}
}
\description{
The Alternative Ratio is an effect size measure that could be accompanying a one-sample binomial, score or Wald test.It is simply the sample proportion (percentage), divided by the expected population proportion (often set at 0.5)

The Alternative Ratio is only mentioned in the documentation of a program called PASS from NCSS (n.d.), and referred to as Relative Risk by JonB (2015).
}
\details{
If codes and category are not provided the first category will be the first data point.

If codes only are provided the first category in the codes is used.

The formula used is:
\deqn{AR=\frac{p}{\pi}}

\emph{Symbols used}:
\itemize{
\item \eqn{p} is the sample proportion of one of the categories
\item \eqn{\pi} the expected proportion
}
}
\examples{
# Example 1: Numeric list
ex1 = c(1, 1, 2, 1, 2, 1, 2, 1)
es_alt_ratio(ex1)
es_alt_ratio(ex1, p0=0.3)

# Example 2: Text list
ex2 = c("Female", "Male", "Male", "Female", "Male", "Male")
es_alt_ratio(ex2)
es_alt_ratio(ex2, category='Female')
es_alt_ratio(ex2, codes=c('Male', 'Female'))

# Example 3: dataframe
file1 <- "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
df1 <- read.csv(file1, sep=",", na.strings=c("", "NA"))
es_alt_ratio(df1['sex'])
es_alt_ratio(df1['mar1'], codes=c("DIVORCED", "NEVER MARRIED"))

}
\references{
JonB. (2015, October 14). Effect size of a binomial test and its relation to other measures of effect size. StackExchange - Cross Validated. https://stats.stackexchange.com/q/176856

NCSS. (n.d.). Tests for one proportion. In PASS Sample Size Software (pp. 100-1-100–132). Retrieved November 10, 2018, from https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Tests_for_One_Proportion.pdf
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}
}
