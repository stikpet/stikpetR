#' Kruskal-Wallis H Test
#'
#' @param scores the variable with the scores
#' @param groups the variable with the groups
#' @param tiescorr boolean to indicate the use of a ties correction
#' @param approx approximation to use (see details)
#'
#' @returns
#' Returns a dataframe with:
#' \item{n}{the sample size}
#' \item{H}{the H value}
#' \item{pValue}{the two-tailed significance, a.k.a. p-value}
#'
#' Depending on the test used additional items might be added.
#' In case of a chi-square approximation
#' \item{statistic}{if not the same as H, the chi-square value used}
#' \item{df}{the degrees of freedom}
#'
#' In case of a gamma, or beta approximation
#' \item{statistic}{the test statistic used}
#' \item{alpha}{the alpha value used}
#' \item{beta}{the beta value used}
#'
#' In case of a F approximation
#' \item{statistic}{the test statistic used}
#' \item{df1}{the first degrees of freedom}
#' \item{df2}{the second degrees of freedom}
#'
#' @details
#' The formula used in case of no ties correction is (Kruskal & Wallis, 1952, p. 586):
#' \deqn{H = \frac{12}{n\times\left(n + 1\right)}\times\sum_{i = 1}^c \frac{R_i^2}{n_i} - 3\times\left(n + 1\right)}
#' With:
#' \deqn{R_i \sum_{j=1}^{n_i} r_{i,j}}
#'
#' *Symbols used*
#' \itemize{
#' \item \eqn{r_{i,j}} is the rank of the i-th score in category j
#' \item \eqn{k} is the number of categories
#' \item \eqn{n_i} is the number of ranks in the i-th category
#' \item \eqn{n} is the total number of scores
#' }
#'
#' Applying a ties correction will adjust the formula to (Kruskal & Wallis, 1952, p. 586):
#' \deqn{H_{adj} = \frac{H}{1 - \frac{\sum T}{n^3 - n}}}
#' With:
#' \deqn{T_i = t_i^3 - t_i}
#'
#' *Symbols added*
#' *Symbols used*
#' \itemize{
#' \item \eqn{t_i} is how often the i-th unique rank occurs
#' }
#'
#' I have not been able to find an exact distribution for H in R. A good starting point
#' might be Choi et al. (2003) and let me know if you manage.
#'
#' As for approximations, there are quite a few....
#'
#' Kruskal and Wallis (1952) had a few of their own.
#' The most frequently used one (Kruskal & Wallis, 1952, p. 586):
#' \deqn{sig. \approx 1 - \chi\left(H, df\right)}
#' With:
#' \deqn{df = k - 1}
#' This is the one used when *approx="chi2"*
#'
#' Kruskal and Wallis (1952) also propose a incomplete gamma distribution approximation (p. 609):
#' \deqn{sig. \approx 1 - \Gamma\left(H, \alpha, \beta\right)}
#' With:
#' \deqn{\alpha = \frac{\mu^2}{\sigma^2}}
#' \deqn{\beta = \frac{\sigma^2}{\mu}}
#' \deqn{\mu = k - 1}
#' \deqn{\sigma^2 = 2\times\left(k - 1\right) - \frac{2\times\left(3\times k^2 - 6\times k + n\times\left(2\times k^2 - 6\times k + 1\right)\right)}{5\times n\left(n + 1\right)} - \frac{6}{5}\times \sum_{i=1}^k \frac{1}{n_i}}
#' This is the one used when *approx="kw-gamma"*
#'
#' And a chi-square approximation of this gamma approximation with:
#' \deqn{sig. \approx 1 - \chi\left(\chi_a^2, df\right)}
#' \deqn{\chi_a^2 = 2\times\frac{\mu}{\sigma^2}\times H}
#' \deqn{df = 2\times\frac{\mu^2}{\sigma^2}\times H}
#' This is the one used when *approx="kw-gamma-chi2"*
#'
#' They also provided an incomplete beta distribution approximation with (Kruskal & Wallis, 1952, pp. 609-610):
#' \deqn{sig. \approx 1 - B\left(\frac{H}{M}, \alpha, \beta \right)}
#' \deqn{M = \frac{n^3 - \sum_{i=1}^k n_i^3}{n\times\left(n + 1\right)}}
#' \deqn{\alpha = \frac{df_1}{2}}
#' \deqn{\beta = \frac{df_2}{2}}
#' \deqn{df_1 = \mu\times\frac{\mu\times\left(M - \mu\right) - \sigma^2}{\frac{1}{2}\times M\times\sigma^2}}
#' \deqn{df_2 = df_1\times\frac{M - \mu}{\mu}}
#' This is the one used when *approx="kw-beta"*
#'
#' And a F-distribution approximation of this using:
#' \deqn{sig. \approx 1 - F\left(F_H, df_1, df_2\right)}
#' \deqn{F_H = \frac{H\times\left(M - \mu\right)}{\mu\times\left(M - H\right)}}
#' This is the one used when *approx="kw-beta-f"*
#'
#' Wallace (1959) used a few different approximations using the F and Beta distributions.
#' First the test values are adjuste (Wallace, 1959, p. 226)
#' \deqn{B_2 = \frac{H}{n - 1}}
#' \deqn{F_2 = \frac{\left(n - k\right)\times H}{\left(k - 1\right)\times\left(n - 1 - H\right)}}
#' Wallace then has three different approximations.
#' First we define d for each (Wallace, 1959, pp. 226-227):
#' \deqn{d_{I} = \frac{\left(n - k\right)\times\left(k - 1\right) - \sigma^2}{\frac{1}{2}\times\left(n - 1\right)\times\sigma^2}}
#' \deqn{d_{II} = 1 - \frac{6\times\left(n + 1\right)}{5\times\left(n - 1\right)\times\left(n + 1.2\right)}}
#' \deqn{d_{III} = 1}
#'
#' Then the degrees of freedom for the F-tests will be:
#' \deqn{df_1^i = \left(k - 1\right)\times d_i}
#' \deqn{df_2^i = \left(n - k\right)\times d_i}
#'
#' And the test itself:
#' \deqn{sig. = 1 - F\left(F_2, df_1^i, df_2^i\right)}
#' This is used when *approx="wallace-I-f"*, *approx="wallace-II-f"*, or *approx="wallace-III-f"*
#'
#' For the beta approximations the same as before:
#' \deqn{\alpha_i = \frac{df_1^i}{2}}
#' \deqn{\beta_i = \frac{df_2^i}{2}}
#' \deqn{sig. \approx 1 - B\left(B_2, \alpha_i, \beta_i \right)}
#' This is used when *approx="wallace-I-beta"*, *approx="wallace-II-beta"*, or *approx="wallace-III-beta"*
#'
#' Finally, Iman and Davenport (1976) describe an approximation from Sattertwaite as:
#' \deqn{sig. = 1 - F\left(F_2, df_1, df_2\right)}
#' With:
#' \deqn{df_1 = k - 1}
#' \deqn{df_2 = \frac{\left(\sum_{i=1}^k \left(n_i -1\right)\times v_i\right)^2}{\sum_{i=1}^k \frac{\left(\left(n_i - 1\right)\times v_i\right)^2}{n_i - 1}}}
#' \deqn{v_i = \frac{\sum_{j=1}^{n_i} \left(r_{i,j} - \bar{r}_i\right)^2}{n_i -1}}
#' \deqn{\bar{r}_i = \frac{\sum_{j = 1}^{n_i} r_{i,j}}{n_i}}
#' This is the one used when *approx="iman"*
#'
#' @references
#' Choi, W., Lee, J. W., Huh, M.-H., & Kang, S.-H. (2003). An algorithm for computing the exact distribution of the Kruskal–Wallis test. *Communications in Statistics - Simulation and Computation, 32*(4), 1029–1040. https://doi.org/10.1081/SAC-120023876
#'
#' Iman, R. L., & Davenport, J. M. (1976). New approximations to the exact distribution of the kruskal-wallis test statistic. *Communications in Statistics - Theory and Methods, 5*(14), 1335–1348. https://doi.org/10.1080/03610927608827446
#'
#' Kruskal, W. H., & Wallis, W. A. (1952). Use of ranks in one-criterion variance analysis. *Journal of the American Statistical Association, 47*(260), 583–621. https://doi.org/10.1080/01621459.1952.10483441
#'
#' Wallace, D. L. (1959). Simplified beta-approximations to the Kruskal-Wallis H test. *Journal of the American Statistical Association, 54*(285), 225. https://doi.org/10.2307/2282148
#'
#' @author
#' P. Stikker
#'
#' Please visit: https://PeterStatistics.com
#'
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @examples
#' scores = c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, NA, 3, 4, 2, 4, 2, 1, 3, 2, 2, 4, 1, 1, 3, 1, 2, 4, 1, 5, 4, 2, 3, 4, 1, 2, 5, 1, 1, 3, 3, 3, 1, 4, 3, 1, 1, 2, 3, 1)
#' groups = c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="chi2")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="kw-gamma")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="kw-gamma-chi2")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="kw-beta")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="kw-beta-f")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="wallace-I-beta")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="wallace-II-beta")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="wallace-III-beta")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="wallace-I-f")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="wallace-II-f")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="wallace-III-f")
#' ts_kruskal_wallis(scores, groups, tiescorr=TRUE, approx="iman")
#'
#' @export
ts_kruskal_wallis <- function(scores, groups,
                              tiescorr=c(TRUE, FALSE),
                              approx=c("chi2", "kw-gamma", "kw-gamma-chi2", "kw-beta", "kw-beta-f", "wallace-I-beta", "wallace-II-beta", "wallace-III-beta", "wallace-I-f", "wallace-II-f", "wallace-III-f", "iman"))
  {

  if (length(tiescorr)>1) {
    tiescorr=tiescorr[1]
  }

  if (length(approx)>1) {
    approx=approx[1]
  }

  dfr = na.omit(data.frame(scores, groups))

  #add ranks to dataframe
  dfr["r"] = rank(dfr$scores)

  #sample size
  n = nrow(dfr)

  #Sum of ranks and count for each category
  Rc = aggregate(dfr$r, by=list(group=dfr$groups), FUN=sum)$x
  nc = aggregate(dfr$r, by=list(group=dfr$groups), FUN=length)$x

  #Unadjusted H
  H = 12/(n*(n + 1)) * sum(Rc^2/nc) - 3*(n + 1)

  if (tiescorr) {
    #Frequencies of each rank
    t = table(dfr["r"])

    #Ties adjustment
    ts = sum(t**3 - t)

    #Adjusted test statistic
    H = H / (1 - ts/(n**3 - n))
  }

  #Test
  k = length(nc)

  if (approx=="chi2") {
    statistic = H
    df = k - 1
    pValue = pchisq(statistic, df, lower.tail = FALSE)
    results <- data.frame(n, H, df, pValue)
  }
  else {
    E = k - 1
    V = 2*(k - 1) - 2*(3*k**2 - 6*k + n*(2*k**2 - 6*k + 1))/(5*n*(n + 1)) - 6/5*sum(1/nc)
    M = (n**3 - sum(nc**3))/(n*(n + 1))
    if(approx=="kw-gamma"){
      #gamma approximation (Kruskal & Wallace, 1952, p. 609)
      alpha = E**2/V
      beta = V/E
      pValue = 1 - pgamma(q=H, shape=alpha, scale=beta)
      results <- data.frame(n, H, alpha, beta, pValue)
    }
    else if(approx=="kw-gamma-chi2"){
      #approximate gamma approximation (Kruskal & Wallace, 1952, p. 609)
      chiA = 2*H*E/V
      dfA = 2*E**2/V
      pValue = pchisq(chiA, dfA, lower.tail = FALSE)

      statistic=chiA
      df = dfA
      results <- data.frame(n, H, statistic, df, pValue)
    }

    else if(approx=="kw-beta"){
      #beta approximation (Kruskal & Wallace, 1952, pp. 609-610)

      f1 = E*((E*(M - E) - V)/(0.5*M*V))
      f2 = ((M - E)/E)*f1
      alpha = 0.5*f1
      beta = 0.5*f2

      statistic = H/M
      pValue = 1 - pbeta(q=statistic, shape1=alpha, shape2=beta)

      results <- data.frame(n, H, statistic, alpha, beta, pValue)

    }
    else if(approx=="kw-beta-f"){
      #approximate with F (Kruskal & Wallace, 1952, pp. 609-610)
      Fvalue = H*(M - E)/(E*(M - H))
      df1 = E*((E*(M - E) - V)/(0.5*M*V))
      df2 = ((M - E)/E)*df1
      pValue = 1 - pf(q=Fvalue, df1, df2, lower.tail = TRUE)

      statistic = Fvalue
      results <- data.frame(n, H, statistic, df1, df2, pValue)
    }
    else {
      d= 1
      if(approx=="wallace-I-beta" || approx=="wallace-I-f"){
        #adjustment to beta approximation I (Wallace, 1959, p. 226)
        d = ((n - k)*(k - 1) - V) / (0.5*(n - 1)*V)}
      else if(approx=="wallace-II-beta" || approx=="wallace-II-f"){
        #adjustment to beta approximation II (Wallace, 1959, p. 227)
        d = 1 - 6/5*(n + 1)/(n - 1)*1/(n + 1.2)}
      else if(approx=="wallace-III-beta" || approx=="wallace-III-f"){
        #adjustment to beta approximation III (Wallace, 1959, p. 227)
        d = 1}

      df1 = (k - 1)*d
      df2 = (n - k)*d

      if(approx %in% c("wallace-I-beta", "wallace-II-beta", "wallace-III-beta")){
        B2 = H/(n - 1)
        alpha = 0.5*df1
        beta = 0.5*df2
        pValue = 1 - pbeta(q=B2, shape1=alpha, shape2=beta)

        statistic = B2
        results <- data.frame(n, H, statistic, alpha, beta, pValue)
      }
      else if (approx %in% c("wallace-I-f", "wallace-II-f", "wallace-III-f", "iman")){
        F2 = (n - k)*H / ((k - 1)*(n - 1 - H))
        pValue = 1 - pf(q=F2, df1, df2, lower.tail = TRUE)

        statistic = F2
        results <- data.frame(n, H, statistic, df1, df2, pValue)
      }

      if (approx=="iman") {
        #Approximation by Iman and Davenport (1976, p. 1338, based on Satterthwaite (1941, 1946))
        RcAvg = Rc/nc
        cats = aggregate(dfr$r, by=list(group=dfr$groups), FUN=sum)$group
        for (i in 1:n) {
          for (j in 1:k) {
            if (cats[j]==dfr[i,2]) {
              dfr[i, 4] = (dfr$r[i] - RcAvg[j])**2
            }
          }
        }

        vi = aggregate(dfr$V4, by=list(group=dfr$groups), FUN=sum)
        vi = vi$x / (nc - 1)
        df1 = k - 1
        df2 = sum((nc - 1)*vi)**2 / (sum(((nc - 1)*vi)**2/(nc - 1)))
        pValue = 1 - pf(q=F2, df1, df2, lower.tail = TRUE)

        statistic = F2
        results <- data.frame(n, H, statistic, df1, df2, pValue)
      }
    }
  }
  return(results)
}
