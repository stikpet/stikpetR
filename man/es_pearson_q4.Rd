% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eff_size_pearson_q4.R
\name{es_pearson_q4}
\alias{es_pearson_q4}
\title{Pearson Q4}
\usage{
es_pearson_q4(var1, var2)
}
\arguments{
\item{var1}{A vector with the binary data from the first variable}

\item{var2}{A vector with the binary data from the second variable}
}
\value{
Pearson Q4
}
\description{
Pearson Q4
}
\details{
This is an approximation for a tetrachoric correlation coefficient.

The formula used (Pearson, 1900, p. 16):
\deqn{Q_4 = \sin\left(\frac{\pi}{2} \times \frac{1}{1 + \frac{2\times b\times c\times n}{\left(a\times d -b\times c\right) \times \left(b+c\right)}}\right)}

\emph{Symbols used:}
\itemize{
\item \eqn{a} the count in the top-left cell of the cross table
\item \eqn{b} the count in the top-right cell of the cross table
\item \eqn{c} the count in the bottom-left cell of the cross table
\item \eqn{d} the count in the bottom-right cell of the cross table
\item \eqn{n} the sum of all counts
}

Note that Pearson Q2 is the same as Yule Q and Pearson Q3 is the same as Yule r.
}
\examples{
bin1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
bin2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
es_pearson_q4(bin1, bin2)

}
\references{
Pearson, K. (1900). Mathematical contributions to the theory of evolution. VII. On the correlation of characters not quantitatively measurable. \emph{Philosophical Transactions of the Royal Society of London, 195}, 1–405. https://doi.org/10.1098/rsta.1900.0022
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
