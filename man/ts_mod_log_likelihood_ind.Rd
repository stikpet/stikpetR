% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_mod_log_likelihood_ind.R
\name{ts_mod_log_likelihood_ind}
\alias{ts_mod_log_likelihood_ind}
\title{Mod-Log Likelihood Test of Independence}
\usage{
ts_mod_log_likelihood_ind(
  field1,
  field2,
  categories1 = NULL,
  categories2 = NULL,
  cc = NULL
)
}
\arguments{
\item{field1}{list or dataframe with the first categorical field}

\item{field2}{list or dataframe with the second categorical field}

\item{categories1}{optional list with order and/or selection for categories of field1}

\item{categories2}{optional list with order and/or selection for categories of field2}

\item{cc}{optional methdod for continuity correction. Either NULL (default), "yates", "pearson", "williams".}
}
\value{
A dataframe with:
\item{n}{the sample size}
\item{n rows}{number of categories used in first field}
\item{n col.}{number of categories used in second field}
\item{statistic}{the test statistic (chi-square value)}
\item{df}{the degrees of freedom}
\item{p-value}{the significance (p-value)}
\item{min. exp.}{the minimum expected count}
\item{prop. exp. below 5}{proportion of cells with expected count less than 5}
\item{test}{description of the test used}
}
\description{
Mod-Log Likelihood Test of Independence
}
\details{
The formula used (Cressie & Read, 1984, p. 441):
\deqn{MG=2\times\sum_{i=1}^{r}\sum_{j=1}^c\left(E_{i,j}\times ln\left(\frac{E_{i,j}}{F_{i,j}}\right)\right)}
\deqn{df = \left(r - 1\right)\times\left(c - 1\right)}
\deqn{sig. = 1 - \chi^2\left(MG,df\right)}

With:
\deqn{n = \sum_{i=1}^r \sum_{j=1}^c F_{i,j}}
\deqn{E_{i,j} = \frac{R_i\times C_j}{n}}
\deqn{R_i = \sum_{j=1}^c F_{i,j}}
\deqn{C_j = \sum_{i=1}^r F_{i,j}}

Cressie and Read (1984) is not the original source, but the source where I found
the formula.

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
 
nom1 <- c("female", "female","female","female","female","female","female",
"female", "female","female","female", "male", "male", "male", "male", "male", 
"male", "male", "male", "male", "male", "male", "male", "male", "male", "male", 
"male","male", "male", "male", "male", "male", "male", "male", "male", "male", 
"male", "male", "male", "male", "male", "male")
nom2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other",
"nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl",
"other", "other", "other", "other", "other", "other", "other", "other", "other", 
"other", "other", "other", "other", "other", "other")
ts_g_ind(nom1, nom2)
ts_g_ind(nom1, nom2, cc="yates")
ts_g_ind(nom1, nom2, cc="pearson")
ts_g_ind(nom1, nom2, cc="williams")
 
}
\references{
Cressie, N., & Read, T. R. C. (1984). Multinomial goodness-of-fit tests. \emph{Journal of the Royal Statistical Society: Series B (Methodological), 46}(3), 440-464. https://doi.org/10.1111/j.2517-6161.1984.tb01318.x

Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 x 2 table. \emph{Biometrika, 34}(1/2), 139-167. https://doi.org/10.2307/2332518

Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. \emph{Biometrika, 63}(1), 33-37. https://doi.org/10.2307/2335081

Yates, F. (1934). Contingency tables involving small numbers and the chi square test. \emph{Supplement to the Journal of the Royal Statistical Society, 1}(2), 217-235. https://doi.org/10.2307/2983604
}
\author{
P. Stikker

Please visit: https://PeterStatistics.com

YouTube channel: https://www.youtube.com/stikpet
}
