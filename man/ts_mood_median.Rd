% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_mood_median.R
\name{ts_mood_median}
\alias{ts_mood_median}
\title{Mood Median Test}
\usage{
ts_mood_median(
  catField,
  ordField,
  categories = NULL,
  levels = NULL,
  test = "pearson",
  cc = c(NULL, "yates", "pearson", "williams"),
  lambd = 2/3
)
}
\arguments{
\item{catField}{vector with categories}

\item{ordField}{vector with the scores}

\item{categories}{vector, optional. the categories to use from catField}

\item{levels}{vector, optional. the levels or order used in ordField.}

\item{test}{string, optional. the test of independence to use. Default is "pearson". Other options are "pearson", "fisher", "freeman-tukey", "g", "mod-log", "neyman", "power"}

\item{cc}{: string, optional. method for continuity correction. Either NULL (default), "yates", "pearson", "williams"}

\item{lambd}{float or string, optional. either name of test or specific value. Default is "cressie-read" i.e. lambda of 2/3. Only applies to Power Divergence test. Other options include float, "cressie-read", "likelihood-ratio", "mod-log", "pearson", "freeman-tukey", "neyman"}
}
\value{
A dataframe with the results of the specified test.
}
\description{
This test looks if the median from different categories would be the same in the population. If not, at least one is different then at least one other category. A Kruskal-Wallis test (see ts_kruksal_wallis()) is very similar but checks the average ranks instead of median.

The test only looks at the number of scores above the overall median and those that are equal or below. A cross table is made with each category and the numbers below and above the overall median. From this table a test of independence can be used.
}
\details{
The Mood Median test creates a 2xk cross table, with k being the number of categories. The two rows are one for the number of scores in that category that are above the overall median, and the second row the number of scores in that category that are equal or below the overall median.

A chi-square test of independence on this cross table can then be performed. There are quite some different options for this:

\itemize{
\item "pearson", will perform a Pearson chi-square test of independence using the ts_pearson_ind() function.
\item "fisher", will perform a Fisher exact test using the ts_fisher() function, but only if there are 2 categories, if there are more the test will be set to "pearson"
\item "freeman-tukey", will perform a Freeman-Tukey test of independence using the ts_freeman_tukey_ind() function
\item "g", will perform a G test of independence using the ts_g_ind() function
\item "mod-log", will perform a Mod-Log Likelihood test of independence using the ts_mod_log_likelihood_ind() function
\item "neyman", will perform a Neyman test of independence using the ts_neyman_ind() function
\item "power", will perform a Power Divergence test of independence using the ts_powerdivergence_ind() function.
}

The formula using the default Pearson test is:
\deqn{\chi_{M}^2 = \sum_{i=1}^2 \sum_{j=1}^k \frac{\left(F_{i,j}-E_{i,j}\right)^2}{E_{i,j}}}
\deqn{df = k - 1}
\deqn{sig. = 1 - \chi^2\left(\chi_{M}^2, df\right)}

With:
\deqn{E_{i,j} = \frac{R_i \times C_j}{n}}
\deqn{R_i = \sum_{j=1}^k F_{i,j}}
\deqn{C_j = \sum_{i=1}^2 F_{i,j}}
\deqn{n = \sum_{i=1}^2 \sum_{j=1}^k F_{i,j} = \sum_{i=1}^2 R_i = \sum_{j=1}^k C_j}

The original source for the formula is most likely Mood (1950), but the ones shown are based on Brown and Mood (1951).

\emph{Symbols used:}

\itemize{
\item \eqn{k}, the number of categories (columns)
\item \eqn{F_{1,j}}, the number of scores is category j that are above the overall median
\item \eqn{F_{2,j}}, the number of scores is category j that are equal or below the overall median
\item \eqn{E_{i,j}}, the expected count in row i and column j.
\item \eqn{R_i}, the row total of row i
\item \eqn{C_j}, the column total of column j
\item \eqn{n}, the overall total.
\item \eqn{df}, the degrees of freedom
\item \eqn{\chi^2\left(\dots\right)}, the cumulative distribution function of the chi-square distribution.
}
}
\references{
Brown, G. W., & Mood, A. M. (1951). On median tests for linear hypotheses. Proceedings of the Second Berkeley Symposium on Mathematical Statistics and Probability, 2, 159-167.

Mood, A. M. (1950). \emph{Introduction to the theory of statistics}. McGraw-Hill.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
