% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_powerdivergence_gof.R
\name{ts_powerdivergence_gof}
\alias{ts_powerdivergence_gof}
\title{Power Divergence Goodness-of Fit Tests}
\usage{
ts_powerdivergence_gof(
  data,
  expCounts = NULL,
  lambd = c("cressie-read", "g", "mod-log", "freeman-tukey", "neyman"),
  cc = c("none", "yates", "pearson", "williams")
)
}
\arguments{
\item{data}{A vector with the data}

\item{expCounts}{Optional counts according to null hypothesis}

\item{lambd}{Optional either name of test or specific value. Default is "cressie-read" i.e. lambda of 2/3}

\item{corr}{Optional correction to be used.}
}
\value{
Dataframe with:
\item{statistic}{the chi-square statistic}
\item{df}{the degrees of freedom}
\item{pValue}{two-sided p-value}
\item{minExp}{the minimum expected count}
\item{propBelow5}{the proportion of expected counts below 5}
\item{testUsed}{a description of the test used}
}
\description{
A test that can be used with a single nominal variable, to test if the probabilities in all the categories are equal (the null hypothesis)

There are quite a few tests that can do this. Perhaps the most commonly used is the Pearson chi-square test (\eqn{\chi^2}), but also an exact multinomial, G-test (\eqn{G^2}), Freeman-Tukey (\eqn{T^2}), Neyman (\eqn{NM^2}), Mod-Log Likelihood (\eqn{GM^2}), and Freeman-Tukey-Read test are possible.

Cressie and Read (1984, p. 463) noticed how the \eqn{\chi^2}, \eqn{G^2}, \eqn{T^2}, \eqn{NM^2} and \eqn{GM^2} can all be captured with one general formula. The additional variable lambda (\eqn{\lambda}) was then investigated, and they settled on a \eqn{\lambda} of 2/3.

By setting \eqn{\lambda} to different values, we get the different tests:
\itemize{
\item{\eqn{\lambda = 1}}{Pearson chi-square}
\item{\eqn{\lambda = 0}}{G/Wilks/Likelihood-Ratio}
\item{\eqn{\lambda = -\frac{1}{2}}}{Freeman-Tukey}
\item{\eqn{\lambda = -1}}{Mod-Log-Likelihood}
\item{\eqn{\lambda = -2}}{Neyman}
\item{\eqn{\lambda = \frac{2}{3}}}{Cressie-Read}
}
}
\details{
The formula used is (Cressie & Read, 1984, p. 442):
\deqn{\chi_{C}^{2} = \begin{cases} 2\times\sum_{i=1}^{r}\sum_{j=1}^c\left(F_{i,j}\times ln\left(\frac{F_{i,j}}{E_{i,j}}\right)\right) & \text{ if } \lambda=0 \\ 2\times\sum_{i=1}^{r}\sum_{j=1}^c\left(E_{i,j}\times ln\left(\frac{E_{i,j}}{F_{i,j}}\right)\right) & \text{ if } \lambda=-1 \\ \frac{2}{\lambda\times\left(\lambda + 1\right)} \times \sum_{i=1}^{r}\sum_{j=1}^{c} F_{i,j}\times\left(\left(\frac{F_{i,j}}{E_{i,j}}\right)^{\lambda} - 1\right) & \text{ else } \end{cases}}
\deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
\deqn{sig. = 1 - \chi^2\left(\chi_{C}^{2},df\right)}

With:
\deqn{n = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
\deqn{E_{i,j} = \frac{R_i\times C_j}{n}}
\deqn{R_i = \sum_{j=1}^c F_{i,j}}
\deqn{C_j = \sum_{i=1}^r F_{i,j}}

\emph{Symbols used:}
\itemize{
\item \eqn{r} the number of categories in the first variable (the number of rows)
\item \eqn{c} the number of categories in the second variable (the number of columns)
\item \eqn{F_{i,j}} the observed count in row i and column j
\item \eqn{E_{i,j}} the expected count in row i and column j
\item \eqn{R_i} the i-th row total
\item \eqn{C_j} the j-th column total
\item \eqn{n} the sum of all counts
\item \eqn{\chi^2\left(\dots\right)}	the chi-square cumulative density function
}

Cressie and Read (1984, p. 463) suggest to use \eqn{\lambda = \frac{2}{3}},  which
is therefor the default in this function.

The \strong{Pearson chi-square statistic} can be obtained by setting \eqn{\lambda = 1}. Pearson's original
formula is (Pearson, 1900, p. 165):
\deqn{\chi_{P}^2 = \sum_{i=1}^r \sum_{j=1}^c \frac{\left(F_{i,j} - E_{i,j}\right)^2}{E_{i,j}}}

The \strong{Freeman-Tukey test} has as a formula (Bishop et al., 2007, p. 513):
\deqn{T^2 = 4\times\sum_{i=1}^r \sum_{j=1}^c \left(\sqrt{F_{i,j}} - \sqrt{E_{i,j}}\right)^2}
This will be same as setting lambda to \eqn{-\frac{1}{2}}. Note that the source for the formula is often quoted to be from Freeman and Tukey (1950)
but couldn't really find it in that article.

\strong{Neyman test} formula was very similar to Pearson's, but the observed and expected counts swapped (Neyman, 1949, p. 250):
\deqn{\chi_{N}^2 = \sum_{i=1}^r \sum_{j=1}^c \frac{\left(E_{i,j} - F_{i,j}\right)^2}{F_{i,j}}}
This will be same as setting lambda to \eqn{-2}.

The Yates correction (yates) is calculated using (Yates, 1934, p. 222):

Use instead of \eqn{F_{i,j}} the adjusted version defined by:
\deqn{F_{i,j}^\ast = \begin{cases} F_{i,j} - 0.5 & \text{ if } F_{i,j}>E_{i,j}  \\ F_{i,j} & \text{ if } F_{i,j}= E_{i,j}\\ F_{i,j} + 0.5 & \text{ if } F_{i,j}<E_{i,j} \end{cases}}

The Pearson correction (pearson) is calculated using (E.S. Pearson, 1947, p. 157):
\deqn{\chi_{PP}^2 = \chi_{P}^{2}\times\frac{n - 1}{n}}

The Williams correction (williams) is calculated using (Williams, 1976, p. 36):
\deqn{\chi_{PW}^2 = \frac{\chi_{P}^2}{q}}
With:
\deqn{q = 1 + \frac{\left(n\times\left(\sum_{i=1}^r \frac{1}{R_i}\right)-1\right) \times \left(n\times\left(\sum_{j=1}^c \frac{1}{C_j}\right)-1\right)}{6\times n\times df}}
}
\examples{
 
nom1 <- c("female", "female","female","female","female","female","female","female", "female","female","female", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male","male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", "male")
nom2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other", "other")
ts_powerdivergence(nom1)
ts_powerdivergence(nom1, nom2)

}
\references{
Bishop, Y. M. M., Fienberg, S. E., & Holland, P. W. (2007). \emph{Discrete multivariate analysis}. Springer.

Cressie, N., & Read, T. R. C. (1984). Multinomial goodness-of-fit tests. \emph{Journal of the Royal Statistical Society: Series B (Methodological), 46}(3), 440–464. https://doi.org/10.1111/j.2517-6161.1984.tb01318.x

Freeman, M. F., & Tukey, J. W. (1950). Transformations related to the angular and the square root. \emph{The Annals of Mathematical Statistics, 21}(4), 607–611. https://doi.org/10.1214/aoms/1177729756

Neyman, J. (1949). Contribution to the theory of the chi-square test. Berkeley Symposium on Math. Stat, and Prob, 239–273. https://doi.org/10.1525/9780520327016-030

Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. \emph{Biometrika, 34}(1/2), 139–167. https://doi.org/10.2307/2332518

Pearson, K. (1900). On the criterion that a given system of deviations from the probable in the case of a correlated system of variables is such that it can be reasonably supposed to have arisen from random sampling. \emph{Philosophical Magazine Series 5, 50}(302), 157–175. https://doi.org/10.1080/14786440009463897

Wilks, S. S. (1938). The large-sample distribution of the likelihood ratio for testing composite hypotheses. \emph{The Annals of Mathematical Statistics, 9}(1), 60–62. https://doi.org/10.1214/aoms/1177732360

Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. \emph{Biometrika, 63}(1), 33–37. https://doi.org/10.2307/2335081

Yates, F. (1934). Contingency tables involving small numbers and the chi square test. \emph{Supplement to the Journal of the Royal Statistical Society, 1}(2), 217–235. https://doi.org/10.2307/2983604
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}
}