% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poho_pairwise_t.R
\name{ph_pairwise_t}
\alias{ph_pairwise_t}
\title{Post-Hoc Pairwise Student T}
\usage{
ph_pairwise_t(nomField, scaleField, categories = NULL)
}
\arguments{
\item{nomField}{the groups variable}

\item{scaleField}{the numeric scores variable}

\item{categories}{vector, optional. the categories to use from catField}
}
\value{
A dataframe with:
\item{category 1}{the first category in the pair}
\item{category 2}{the second category in the pair}
\item{n1}{sample size of first category}
\item{n2}{sample size of second category}
\item{mean 1}{arithmetic mean of scores in first category}
\item{mean 2}{arithmetic mean of scores in second category}
\item{sample diff.}{difference between the two arithmetic means}
\item{hyp diff.}{the hypothesized difference}
\item{statistic}{the test-statistic}
\item{df}{the degrees of freedom}
\item{p-value}{the unadjusted p-value (significance)}
\item{adj. p-value}{the Bonferroni adjusted p-values}
\item{test}{description of test used}
}
\description{
This function performs pairwise independent samples Student t tests, for use after a one-way ANOVA, to determine which categories significantly differ from each other.

It differs slightly in the calculation of the standard error, than the version used by using ph_pairwise_is(nomField, scaleField, isTest = "student"). This version appears to be producing the same results as SPSS shows, when using a Bonferroni correction. SPSS refers to Winer (1962) for their procedures.

A simple Bonferroni correction is also applied.
}
\details{
The formula used:
\deqn{t_{1,2} = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{MS_w \times\left(\frac{1}{n_1}+ \frac{1}{n_2}\right)}}}
\deqn{df_w = n - k}
\deqn{sig. = 2\times\left(1 - T\left(\left|t_{1,2}\right|, df_w\right)\right)}

With:
\deqn{MS_w = \frac{SS_w}{df_w}}
\deqn{SS_w = \sum_{j=1}^k \sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}_j\right)^2}
\deqn{\bar{x}_j = \frac{\sum_{i=1}^{n_j} x_{i,j}}{n_j}}

\emph{Symbols used}
\itemize{
\item \eqn{x_{i,j}}, the i-th score in category j
\item \eqn{n}, the total sample size
\item \eqn{n_j}, the number of scores in category j
\item \eqn{k}, the number of categories
\item \eqn{\bar{x}_j}, the mean of the scores in category j
\item \eqn{MS_w}, the mean square within
\item \eqn{SS_w}, the sum of squares of within (sum of squared deviation of the mean)
\item \eqn{df_w}, the degrees of freedom of within
}

The Bonferroni adjustment is simply:
\deqn{p_{adj} = \min \left(p \times n_{comp}, 1\right)}
\deqn{n_{comp} = \frac{k\times\left(k-1\right)}{2}}

\emph{Symbols used:}

\itemize{
\item \eqn{n_{comp}}, number of comparisons (pairs)
\item \eqn{k}, number of categories
}
}
\references{
Winer, B. J. (1962). \emph{Statistical principles in experimental design}. McGraw Hill.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
