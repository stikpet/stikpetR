% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eff_size_phi.R
\name{es_phi}
\alias{es_phi}
\title{Pearson/Yule Phi Coefficient}
\usage{
es_phi(var1, var2)
}
\arguments{
\item{var1}{A vector with the binary data from the first variable}

\item{var2}{A vector with the binary data from the second variable}
}
\value{
phi coefficient
}
\description{
Pearson/Yule Phi Coefficient
}
\details{
the Pearson Phi coefficient (Pearson, 1900, p. 12), is the same as Yule's Phi (Yule, 1912, p. 596)
Cole's C2 (Cole, 1949, p. 415) and Cohen's w (Cohen, 1988, p. 216). It is also sometimes referred to as the Mean Square Contingency.

This is also the same result as if values 0 and 1 would be used for both variables, and the regular Pearson Correlation calculated

The formula used is (Pearson, 1900, p. 12):
\deqn{\phi = \frac{a\times d - b\times c}{\sqrt{R_1\times R_2 \times C_1 \times C_2}}}

\emph{Symbols used:}
\itemize{
\item \eqn{a} the count in the top-left cell of the cross table
\item \eqn{b} the count in the top-right cell of the cross table
\item \eqn{c} the count in the bottom-left cell of the cross table
\item \eqn{d} the count in the bottom-right cell of the cross table
\item \eqn{R_i} the sum of counts in the i-th row
\item \eqn{C_i} the sum of counts in the i-th column
}

Note that Cohen w did not limit the size of the table, but uses the same formula.
For classification/qualification the \emph{th_cohen_w(phi)} function could be used.
}
\examples{
bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
es_phi(bin1, bin2)

}
\references{
Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). L. Erlbaum Associates.

Cole, L. C. (1949). The measurement of interspecific associaton. \emph{Ecology, 30}(4), 411–424. https://doi.org/10.2307/1932444

Pearson, K. (1900). Mathematical Contributions to the Theory of Evolution. VII. On the Correlation of Characters not Quantitatively Measurable. \emph{Philosophical Transactions of the Royal Society of London. Series A, Containing Papers of a Mathematical or Physical Character}, 195, 1–405.

Yule, G. U. (1912). On the methods of measuring association between two attributes. \emph{Journal of the Royal Statistical Society, 75}(6), 579–652. https://doi.org/10.2307/2340126
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
