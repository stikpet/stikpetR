% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_kruskal_wallis.R
\name{ts_kruskal_wallis}
\alias{ts_kruskal_wallis}
\title{Kruskal-Wallis H Test}
\usage{
ts_kruskal_wallis(
  catField,
  ordField,
  categories = NULL,
  levels = NULL,
  method = c("chi2", "kw-gamma", "kw-gamma-chi2", "kw-beta", "kw-beta-f",
    "wallace-I-beta", "wallace-II-beta", "wallace-III-beta", "wallace-I-f",
    "wallace-II-f", "wallace-III-f", "iman"),
  tiescorr = c(TRUE, FALSE)
)
}
\arguments{
\item{catField}{vector with categories}

\item{ordField}{vector with the scores}

\item{categories}{vector, optional. the categories to use from catField}

\item{levels}{vector, optional. the levels or order used in ordField.}

\item{method}{string, optional. the specific variation of the test to use. Default is "chi2". Options are: : "chi2", "kw-gamma", "kw-gamma-chi2", "kw-beta", "kw-beta-f", "wallace-f1", "wallace-f2", "wallace-f3", "wallace-beta1", "wallace-beta2", "wallace-beta3", "ids"}

\item{tiescorr}{boolean, optional. use of a ties correction. Default is True.}
}
\value{
Returns a dataframe with:
\item{n}{the sample size}
\item{H}{the H value}

Depending on the test used additional items might be added.
In case of a chi-square approximation
\item{statistic}{if not the same as H, the chi-square value used}
\item{df}{the degrees of freedom}

In case of a gamma, or beta approximation
\item{statistic}{the test statistic used}
\item{alpha}{the alpha value used}
\item{beta}{the beta value used}

In case of a F approximation
\item{statistic}{the test statistic used}
\item{df1}{the first degrees of freedom}
\item{df2}{the second degrees of freedom}
}
\description{
This test is an extension of the Mann-Whitney U test (see ts_mann_whitney()) to more than two categories. It is also seen as the non-parametric version of the one-way ANOVA (see ts_fisher_owa()).

The test can indicate if any of the scores in one or more categories, has a significant different mean rank than one or more of the other categories. More strickly the null hypothesis is that the probability of a randomly selected case having a score greater than a random score from the other category is 50\% (Divine et al., p. 286).

Alternative there is a Mood Median test (see ts_mood_median()).

To pin-point which category or categories differ significantly, a post-hoc analysis could be used.
}
\details{
\strong{The H value}

The formula used is (Kruskal & Wallis, 1952, p. 586):
\deqn{H = \frac{12}{n\times\left(n + 1\right)}\times\sum_{i=1}^k \frac{R_i^2}{n_i}-3\times\left(n+1\right)}

With:
\deqn{R_i = \sum_{j=1}^{n_i} r_{i,j}}

The ties correction (Kruskal & Wallis, 1952, p. 586):
\deqn{H_{adj} = \frac{H}{1-\frac{\sum T}{n^3-n}}}

With:
\deqn{T_j = t_j^3 - t_j}

Or alternatively:
\deqn{H_{adj} = \left(n - 1\right)\times\frac{\sum_{i=1}^k n_i\left(\bar{r}_i - \bar{r}\right)^2}{\sum_{i=1}^k \sum_{j=1}^{n_i}\left(r_{i,j}-\bar{r}\right)^2}}

With:
\deqn{\bar{r}_i = \frac{R_i}{n_i}}
\deqn{\bar{r} = \frac{\sum_{i=1}^k \bar{r}_i}{\sum_{i=1}^k n_i}}

\strong{The Test}

\emph{"chi2", Kruskal-Wallis Chi-Square Approximation}
\deqn{sig. \approx 1 - \chi^2\left(H, df\right)}
\deqn{df = k - 1}

\emph{"kw-gamma", Kruskal-Wallis incomplete gamma approximation} (Kruskal & Wallis, 1952, p. 609)
\deqn{sig. \approx 1 - \gamma\left(H, \alpha, \beta\right)}
\deqn{\alpha = \frac{\mu^2}{\sigma^2}}
\deqn{\beta = \frac{\sigma^2}{\mu}}
\deqn{\mu = k -1}
\deqn{\sigma^2 = 2\times\left(k - 1\right) - \frac{2\times\left(k\times k^2 - 6\times k + n\times\left(2\times k^2-6\times k+1\right)\right)}{5\times n\times\left(n + 1\right)} - \frac{6}{5}\times\sum_{i=1}^k \frac{1}{n_i}}

\emph{"kw-gamma-chi2", Kruskal-Wallis Chi-square approximation of gamma approximation}
\deqn{sig. = 1 - \chi^2\left(\chi_a^2, df\right)}
\deqn{\chi_a^2 = \frac{2\times\mu}{\sigma^2}\times H}
\deqn{df = 2\times\frac{\mu^2}{\sigma^2}}

\emph{"kw-beta", Kruskal-Wallis incomplete Beta distribution approximation} (Kruskal & Wallis, 1952, p. 609)
\deqn{sig. = 1 - \beta\left(\frac{H}{M}, \alpha, \beta\right)}
\deqn{M = \frac{n^3 - \sum_{i=1}^k n_i^3}{n\times\left(n + 1\right)}}
\deqn{\alpha = df_1\times\frac{1}{2}}
\deqn{\beta = df_2\times\frac{1}{2}}
\deqn{df_1 = \mu\times\frac{\mu\times\left(M-\mu\right)-\sigma^2}{\frac{1}{2}\times M\times\sigma^2}}
\deqn{df_2 = df_1\times\frac{M-\mu}{\mu}}

\emph{"kw-beta-f", F-approximation of the Kruskal-Wallis incomplete Beta distribution approximation} (Kruskal & Wallis, 1952, p. 610)
\deqn{sig. \approx 1 - F\left(F_{\alpha}, df_1, df_2\right)}
\deqn{F_{\alpha} = \frac{H\times\left(M - \mu\right)}{\mu\times\left(M - H\right)}}

\emph{Wallace F distribution approximations} (Wallace, 1959, p. 226)
\deqn{sig. = 1 - F\left(F_2, df_1^i, df_2^i\right)}

With:
\deqn{F_2 = \frac{\left(n - k\right)\times H}{\left(k - 1\right)\times\left(n - 1 - H\right)}}
\deqn{df_1^i = \left(k - 1\right)\times d_i}
\deqn{df_2^i = \left(n - k\right)\times d_i}

\emph{Wallace Beta distribution approximations} (Wallace, 1959, p. 226)
\deqn{sig. \approx 1 - \beta\left(B_2, \alpha, \beta\right)}
\deqn{B_2 = \frac{H}{n-1}}
\deqn{\alpha = df_1\times\frac{1}{2}}
\deqn{\beta = df_2\times\frac{1}{2}}

\emph{"wallace-f1"} and \emph{"wallace-b1"}
\deqn{d_i = \frac{\left(n - k\right)\times\left(k - 1\right)-\sigma^2}{\frac{1}{2}\times\left(n - 1\right)\times \sigma^2}}

\emph{"wallace-f2"} and \emph{"wallace-b2"}
\deqn{d_i = 1 - \frac{6\times\left(n + 1\right)}{5\times\left(n - 1\right)\times\left(n+1.2\right)}}

\emph{"wallace-f3"} and \emph{"wallace-b3"}
\deqn{d_i = 1}

\emph{"ids", Iman-Davenport Satterwaite approximation} (Iman & Davenport, 1976, p. 1338)
\deqn{sig. = 1 - F\left(F_2, df_1, df_2\right)}

With:
\deqn{df_1 = k - 1}
\deqn{df_2 = \frac{\left(\sum_{i=1}^k\left(n_i-1\right)\times v_i\right)^2}{\sum_{i=1}^k\frac{\left(\left(n_i-1\right)\times v_i\right)^2}{n_i-1}}}
\deqn{v_i = \frac{\sum_{j=1}^{n_i}\left(r_{i,j}-\bar{r}_i\right)^2}{n_i-1}}
\deqn{\bar{r}_i = \frac{\sum_{j=1}^{n_i} r_{i,j}}{n_i}}

\emph{Symbols used:}
\itemize{
\item \eqn{k}, the number of categories
\item \eqn{t_j}, the frequency of the j-th unique rank.
\item \eqn{n}, the total sample size
\item \eqn{n_i}, the number of scores in category i
\item \eqn{r_{i,j}}, the rank of the j-th score in category i
\item \eqn{R_i}, the sum of the ranks in category i
\item \eqn{\bar{r}_i}, the average of the ranks in category i
\item \eqn{\bar{r}}, the average of all ranks
\item \eqn{\chi^2\left(\dots\right)}, the cumulative distribution function of the chi-square distribution.
\item \eqn{F\left(\dots\right)}, the cumulative distribution function of the F distribution.
\item \eqn{\beta\left(\dots\right)}, the cumulative distribution function of the beta distribution.
}

I have not been able to find an exact distribution for H in R. A good starting point
might be Choi et al. (2003) and let me know if you manage.
}
\references{
Iman, R. L., & Davenport, J. M. (1976). New approximations to the exact distribution of the kruskal-wallis test statistic. \emph{Communications in Statistics - Theory and Methods, 5}(14), 1335-1348. doi:10.1080/03610927608827446

Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis. \emph{Journal of the American Statistical Association, 47}(260), 583-621. doi:10.1080/01621459.1952.10483441

Wallace, D. L. (1959). Simplified beta-approximations to the Kruskal-Wallis H test. \emph{Journal of the American Statistical Association, 54}(285), 225. doi:10.2307/2282148
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
