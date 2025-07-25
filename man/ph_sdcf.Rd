% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poho_sdcf.R
\name{ph_sdcf}
\alias{ph_sdcf}
\title{Post-Hoc Steel-Dwass-Critchlow-Fligner Test}
\usage{
ph_sdcf(catField, ordField, categories = NULL, levels = NULL)
}
\arguments{
\item{catField}{vector with categories}

\item{ordField}{vector with the scores}

\item{categories}{vector, optional. the categories to use from catField}

\item{levels}{vector, optional. the levels or order used in ordField.}
}
\value{
A dataframe with:
\item{cat. 1}{one of the two categories being compared}
\item{cat. 2}{second of the two categories being compared}
\item{n1}{number of cat. 1. cases in comparison}
\item{n2}{number of cat. 2 cases in comparison}
\item{mean rank 1}{mean rank of cases in cat. 1, based on all cases (incl. categories not in comparison)}
\item{mean rank 2}{mean rank of cases in cat. 2, based on all cases (incl. categories not in comparison)}
\item{statistic}{the q-value of the test}
\item{std. statistic}{the standardized q value}
\item{p-value}{the p-value (significance)}
}
\description{
This can be used as a post-hoc test for a Kruskal-Wallis test (see ts_kruskal_wallis()).

The test compares each possible pair of categories from the catField and their mean rank. The null hypothesis is that these are then equal.

Other post-hoc tests that could be considered are Dunn, Nemenyi, Conover, a pairwise Mann-Whitney U, or pairwise Mood-Median.

Unlike the Dunn, Nemenyi and Conover-Iman test, this test re-calculates the mean ranks for each pair, using only the scores from the two categories.
}
\details{
The formula used (Hollander & Wolfe, 1999, p. 241):
\deqn{q_{1,2} = \frac{\left| R_1 - E_1\right|}{\sqrt{\sigma^2}}}

With:
\deqn{R_1 = \sum_{i=1}^{n_1} r_{i,1}}
\deqn{n_{1,2} = n_1 + n_2}
\deqn{E_1 = \frac{n_1\times\left(n_{1,2} + 1\right)}{2}}
\deqn{\sigma^2 = \frac{n_1\times n_2}{12}\times\left(n_{1,2} + 1 - \frac{T}{n_{1,2}-1}\right)}
\deqn{T = \sum t_j^3 - t_j}

The p-value is then determined using (Critchlow & Fligner, 1991, p. 131):
\deqn{sig. = 1 - Q\left(q_{1,2}, k, df=\infty\right)}

Note that while looking at the R-code for this, posted by Shigenobu (n.d.), who references Nagata and Yoshida (1997), an alternative but same result equation for the variance can be used:
\deqn{\sigma^2 = \frac{n_1\times n_2}{n_{1,2}\times\left(n_{1,2} - 1\right)}\times\left(\sum_{i=1}^{n_1} r_{i,1}^2 + \sum_{i=1}^{n_2} r_{i,2}^2 - \frac{n_{1,2}\times\left(n_{1,2}+1\right)^2}{4}\right)}.

Steel (1960) and Dwass (1960) independently derived the basics for this test. Critchlow and Fligner (1991) added the case for larger samples using the Tukey Range Distribution, and in Hollander and Wolfe (1999) the version used here can be found, which includes a ties correction.

\emph{Symbols used}

\itemize{
\item \eqn{k}, the number of categories
\item \eqn{t_j}, the frequency of the j-th unique rank.
\item \eqn{n_i}, the number of scores in category i
\item \eqn{r_{i,j}}, the rank of the j-th score in category i using only the scores from the two categories in the comparison.
\item \eqn{Q\left(\dots\right)}, the cumulative distribution function of the standardized range distribution.
}
}
\references{
Critchlow, D. E., & Fligner, M. A. (1991). On distribution-free multiple comparisons in the one-way analysis of variance. \emph{Communications in Statistics - Theory and Methods, 20}(1), 127-139. doi:10.1080/03610929108830487

Dwass, M. (1960). \emph{Some k-sample rank-order tests}. In I. Olkin, S. G. Ghurye, W. Hoeffding, W. G. Madow, & H. B. Mann (Eds.), Contributions to probability and statistics; Essays in honor of Harold Hotelling. Stanford University Press.

Hollander, M., & Wolfe, D. A. (1999). \emph{Nonparametric statistical methods} (2nd ed.). Wiley.

Nagata, Y., & Yoshida, M. (1997). \emph{The Basics of Multiple Comparisons in Statistics}. Scientist Co.

Shigenobu. (2004, July 28). Multiple comparisons using the Steel-Dwass method. http://aoki2.si.gunma-u.ac.jp/R/Steel-Dwass.html

Steel, R. G. D. (1960). A rank sum test for comparing all pairs of treatments. Technometrics, 2(2), 197-207. doi:10.1080/00401706.1960.10489894
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
