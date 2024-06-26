% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eff_size_pairwise_bin_ord.R
\name{es_pairwise_bin_ord}
\alias{es_pairwise_bin_ord}
\title{Pairwise Binary-Ordinal Effect Sizes}
\usage{
es_pairwise_bin_ord(
  catField,
  ordField,
  categories = NULL,
  levels = NULL,
  es = "vda"
)
}
\arguments{
\item{catField}{vector with categories}

\item{ordField}{vector with the scores}

\item{categories}{vector, optional. the categories to use from catField}

\item{levels}{vector, optional. the levels or order used in ordField.}

\item{es}{string, optional. the effect size to determine.Either "vda" (default), "rb" or "rosenthal"}
}
\value{
dataframe with
\item{cat1}{label of first category in pair}
\item{cat2}{label of second category in pair}
\item{effect size}{the value of the effect size}
}
\description{
This function determines the effect size for each comparison in a post-hoc analysis of a nominal vs. ordinal variable (e.g. a Kruskal-Wallis test).
}
\details{
The function simply goes over each possible pair of categories from the \emph{catField} (adjusted with \emph{categories} if used). It then runs for only the scores of those two categories the Vargha-Delaney A or (Glass) Rank Biserial (Cliff delta). If the Rosenthal correlation is requested, it will perform the post-hoc Dunn test to obtain the z-statistic.
}
\seealso{
\code{\link{es_vargha_delaney_a}}, Vargha-Delaney A effect size
\code{\link{r_rank_biserial_is}}, rank biserial for independent samples
\code{\link{ph_dunn}}, post-hoc Dunn test, used to obtain z-value for Rosenthal correlation
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
