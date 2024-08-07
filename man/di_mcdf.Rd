% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_multinomial.R
\name{di_mcdf}
\alias{di_mcdf}
\title{Multinomial Cumulative Distribution Function}
\usage{
di_mcdf(F, P, method = "loggamma")
}
\arguments{
\item{F}{list with the observed counts}

\item{P}{list with the probabilities for each category}

\item{method}{optional the calculation method to use. Either "loggamma" (default), "factorial", "gamma", "mprob".}
}
\value{
A float with the requested probability
}
\description{
This is a function for the cumulative multinomial probability. It returns the probability of a distribution as given in F for a sample size of sum of F, where the probability for each category is given as in P, or a distribution even more rare. It is a generalization of the binomial distribution.
}
\details{
The function first determines all possible arrangements over k categories that sum to n, using the \strong{find_combinations()} function. It then uses the \strong{di_pmf()} function to determine the probability for each of these, and sums those that are less or equal to the sample version.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
