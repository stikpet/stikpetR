% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eff_size_odds_ratio.R
\name{es_odds_ratio}
\alias{es_odds_ratio}
\title{Odds Ratio}
\usage{
es_odds_ratio(field1, field2, categories1 = NULL, categories2 = NULL)
}
\arguments{
\item{field1}{: dataframe field with categories for the rows}

\item{field2}{: dataframe field with categories for the columns}

\item{categories1}{: optional list with selection and/or order for categories of field1}

\item{categories2}{: optional list with selection and/or order for categories of field2}
}
\value{
Dataframe with:
\item{OR}{the odds ratio}
\item{n}{the sample size}
\item{statistic}{the test statistic (z-value)}
\item{p-value}{the significance (p-value)}
}
\description{
Determines the odds ratio from a 2x2 table.

Odds can sometimes be reported as 'a one in five odds', but sometimes as 1 : 4. This later notation is less often seen, but means for every one event on the left side, there will be four on the right side.

The Odds is the ratio of that something will happen, over the probability that it will not. For the Odds Ratio, we compare the odds of the first category with the second group.

If the result is 1, it indicates that one variable has no influence on the other. A result higher than 1, indicates the odds are higher for the first category. A result lower than 1, indicates the odds are lower for the first.
}
\details{
The formula used is (Fisher, 1935, p. 50):
\deqn{OR = \frac{a/c}{b/d} = \frac{a\times d}{b\times c}}

\emph{Symbols used:}
\itemize{
\item \eqn{a} the count in the top-left cell
\item \eqn{b} the count in the top-right cell
\item \eqn{c} the count in the bottom-left cell
\item \eqn{d} the count in the bottom-right cell
\item \eqn{\Phi\left(\dots\right)} the cumulative density function of the standard normal distribution
}

As for the test (McHugh, 2009, p. 123):
\deqn{sig. = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}

With:
\deqn{SE = \sqrt{\frac{1}{a} + \frac{1}{b} + \frac{1}{c} + \frac{1}{d}}}
\deqn{z = \frac{\ln{\left(OR\right)}}{SE}}

The p-value is for the null-hypothesis that the population OR is 1.

The term Odds Ratio can for example be found in Cox (1958, p. 222).
}
\section{Alternatives}{


R's \emph{stats} library has a function that also shows an odds ratio: \emph{fisher.test()}
}

\examples{
#Example: dataframe
dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
es_odds_ratio(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))

}
\references{
Cox, D. R. (1958). The regression analysis of binary sequences. \emph{Journal of the Royal Statistical Society: Series B (Methodological), 20}(2), 215-232. https://doi.org/10.1111/j.2517-6161.1958.tb00292.x

Fisher, R. A. (1935). The logic of inductive inference. \emph{Journal of the Royal Statistical Society, 98}(1), 39-82. https://doi.org/10.2307/2342435

McHugh, M. (2009). The odds ratio: Calculation, usage, and interpretation. \emph{Biochemia Medica, 19}(2), 120-126. https://doi.org/10.11613/BM.2009.011
}
\seealso{
\code{\link{th_odds_ratio}}, rules of thumb for odds ratio

\code{\link{es_convert}}, to convert an odds ratio to Yule Q, Yule Y, or Cohen d.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
