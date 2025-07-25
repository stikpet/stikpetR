#' Spearman Rho Distribution
#' @description
#' The Spearman Rank Correlation Coefficient Distribution. Will return a two-tailed p-value
#'
#' This function makes use of the *pspearman* library for exact computations.
#'
#' @param n the number of scores (should be equal in both variables)
#' @param rs the Spearman rank correlation coefficient
#' @param method the test to be used
#' @param iters the number of iterations to use, only applicable if Iman-Conover is used
#' @returns
#' A dataframe with:
#' \item{statistic}{the statistic from the test (only if applicable)}
#' \item{df}{the degrees of freedom (only if applicable)}
#' \item{pValue}{the significance (p-value)}
#'
#' @details
#'
#' The exact distribution is calculated using the following steps:
#' \enumerate{
#' \item Determine all possible permutations of the scores in the first variable
#' \item Determine for each permutation the Spearman rho with the second variable
#' \item Count how often the Spearman rho is above the Spearman rho between the original two variables
#' \item Divide the results by \eqn{n!}
#' }
#'
#' This procedure can be used by using the *he_spearman_exact(ord1, ord2)* function.
#'
#' This function however, makes use of the *pspearman* function from the *pspearman* library.
#' It seems this uses van de Wiel and Bucchianico (2001) method for the exact distribution,
#' which can handle larger sample sizes (up to n = 22).
#'
#' The Student t distribution approximation uses (Kendall & Stuart, 1979, p. 503; Iman & Conover, 1978):
#' \deqn{t_s = r_s\times\sqrt{\frac{n - 2}{1 - r_s^2}}}
#' \deqn{df = n - 2}
#' \deqn{sig = 2\times\left(1 - T\left(\left|t_s\right|, df\right)\right)}
#'
#' Iman and Conover refer to Pitman (1937) for the test.
#'
#' The Fieller's standard normal distribution approximation uses (Fieller et al., 1957, p. 472; Choi, 1977, p. 646):
#' \deqn{z_F = \sqrt{\frac{n - 3}{1.06}}\times\text{atanh}\left(r_s\right)}
#' \deqn{sig = 2\times\left(1 - \Phi\left(\left|z_F\right|\right)\right)}
#'
#' The Old's standard normal distribution approximation uses (Olds, 1938, p. 142; Olds, 1949, p. 117):
#' \deqn{z_O = \frac{x}{ASE}}
#' \deqn{sig = 2\times\left(1 - \Phi\left(\left|z_O\right|\right)\right)}
#' With:
#' \deqn{x = \frac{S}{2} - \frac{n^3 - n}{12}}
#' \deqn{ASE = \sqrt{n - 1}\times\frac{n^2 + n}{12}}
#' \deqn{S = \frac{\left(n^3 - n\right)\times\left(1 - r_s\right)}{6}}
#'
#' A combination of the Student t and Normal approximation by Iman and Conover (1978, p. 272) uses:
#' \deqn{J = \frac{r_s}{2}\times\left(\sqrt{n - 1} + \sqrt{n - 2}{1 - r_s^2}\right)}
#' And reject the null hypothesis if:
#' \deqn{J > J_{crit}}
#' With:
#' \deqn{J_{crit} = \frac{Q\left(\Phi\left(1 - \frac{\alpha}{2}\right)\right) + Q\left(T\left(1 - \frac{\alpha}{2}, df\right)\right)}{2}}
#' \deqn{df = n - 2}
#'
#' The function will use a binary search to find \eqn{alpha} such that \eqn{J = J_{crit}}.
#'
#' One of the more popular methods is Algorithm AS 89 (Best & Roberts, 1975).
#' This is available as a separate helper function *he_AS89(n, S)*.
#'
#' Often in publications the test statistic \eqn{S} is mentioned, this can be
#' defined as:
#' \deqn{S = \sum_{i = 1}^n d_i^2 = \sum_{i=1}^n \left(r_{x_i} - r_{y_i}\right)^2}
#' Which if there are no ties is equal to:
#' \deqn{S = \frac{\left(n^3 - n\right)\times\left(1 - r_s\right)}{6}}
#'
#' @references
#' Best, D. J., & Roberts, D. E. (1975). Algorithm AS 89: The upper tail probabilities of Spearman's rho. *Applied Statistics, 24*(3), 377-379. https://doi.org/10.2307/2347111
#'
#' Choi, S. C. (1977). Tests of equality of dependent correlation coefficients. *Biometrika, 64*(3), 645-647. https://doi.org/10.1093/biomet/64.3.645
#'
#' Fieller, E. C., Hartley, H. O., & Pearson, E. S. (1957). Tests for rank correlation coefficients. I. *Biometrika, 44*(3-4), 470-481. https://doi.org/10.1093/biomet/44.3-4.470
#'
#' Iman, R. L., & Conover, W. J. (1978). Approximations of the critical region for spearman's rho with and without ties present. *Communications in Statistics - Simulation and Computation, 7*(3), 269-282. https://doi.org/10.1080/03610917808812076
#'
#' Kendall, M., & Stuart, A. (1979). *The advanced theory of statistics. Volume 2: Inference and relationship* (4th ed.). Griffin.
#'
#' Olds, E. G. (1938). Distributions of sums of squares of rank differences for small numbers of individuals. *The Annals of Mathematical Statistics, 9*(2), 133-148. https://doi.org/10.1214/aoms/1177732332
#'
#' Olds, E. G. (1949). The 5% significance levels for sums of squares of rank differences and a correction. *The Annals of Mathematical Statistics, 20*(1), 117-118. https://doi.org/10.1214/aoms/1177730099
#'
#' Pitman, E. J. G. (1937). Significance tests which may be applied to samples from any populations. II. The correlation coefficient test. Supplement to the *Journal of the Royal Statistical Society, 4*(2), 225-232. https://doi.org/10.2307/2983647
#'
#' van de Wiel, M. A., & Bucchianico, A. D. (2001). Fast computation of the exact null distribution of Spearman's rho and Page's L statistic for samples with and without ties. *Journal of Statistical Planning and Inference, 92*(1-2), 133-145. https://doi.org/10.1016/S0378-3758(00)00166-X
#'
#' @author
#' P. Stikker
#'
#' Please visit: https://PeterStatistics.com
#'
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @examples
#' n = 12
#' rs = 0.8
#' di_scdf(n, rs)
#' di_scdf(n, rs, method="exact")
#' di_scdf(n, rs, method="t")
#' di_scdf(n, rs, method="z-fieller")
#' di_scdf(n, rs, method="z-olds")
#' di_scdf(n, rs, method="iman-conover")
#' di_scdf(n, rs, method="AS89")
#'
#' @export
di_scdf <- function(n, rs,
                        method=c("t", "z-fieller", "z-olds", "iman-conover", "AS89", "exact"),
                        iters=500){
  if (length(method)>1) {
    method="t"
  }
  df = n - 2

  if (method=="t") {
    #(Kendall & Stuart, 1979, p. 503)
    ts = rs*sqrt((n - 2)/(1 - rs**2))
    pValue = 2*(1 - pt(abs(ts), df))

    statistic = ts
    results = data.frame(statistic, df, pValue)
  }

  else if (method=="z-fieller") {
    #(Fieller et al. 1957, p. 472; Choi, 1977, p. 646)
    zs = atanh(rs)/sqrt(1.06/(n - 3))
    pValue = 2*(1 - pnorm(abs(zs)))

    statistic = zs
    results = data.frame(statistic, pValue)
  }

  else if (method=="z-olds"){
    #(Olds, 1938, p. 142; Olds, 1949, p. 117)
    S = (n**3 - n)*(1 - rs)/6
    x = S/2 - (n**3 - n)/12
    ASE = sqrt(n - 1)*(n*(n+1)/12)
    z = x/ASE
    pValue = 2*(1 - pnorm(abs(z)))

    statistic = z
    results = data.frame(statistic, pValue)
  }

  else if (method=="iman-conover"){
    #(Iman & Conover, 1978, p. 272)
    J = abs(rs)/2*(sqrt(n - 1) + sqrt((n - 2)/ (1 - abs(rs)**2)))
    zpval = 2*(1 - pnorm(abs(J)))
    tpval = 2*(1 - pt(abs(J), df))
    pValue = (zpval + tpval) / 2

    pLow <- 0
    pHigh <- 1
    pValue <- 0.05
    nIter = 1
    repeat {
      zCrit = qnorm(1-pValue/2)
      tCrit = qt(1 - pValue/2, df)
      Jcrit = (zCrit + tCrit)/2
      if (Jcrit == J | nIter == iters) {
        break
      }
      if (Jcrit < J) {
        pHigh <- pValue
        pValue <- (pLow + pValue)/2
      }
      else if (Jcrit > J) {
        pLow <- pValue
        pValue <- (pHigh + pValue)/2
      }
      nIter = nIter + 1
    }
    statistic = J
    results = data.frame(statistic, df, pValue)
  }

  else if (method=="AS89") {
    S = (n**3 - n)*(1 - rs)/6
    statistic = S

    if (S < (n**3 - n)/6) {
      S = (n**3 - n)/3 - S
    }
    pValue = he_AS89(n, S)
    if (pValue > 0.5) {
      pValue = 2*(1 - pValue)
    }
    else{
      pValue = 2*pValue
    }
    results = data.frame(statistic, pValue)
  }

  else if (method=="exact") {
    S = (n**3 - n)*(1 - rs)/6
    if (S > (n**3 - n)/6) {
      pValue = pspearman::pspearman(S, n, lower.tail=FALSE, approximation="exact")
    }
    else{
      pValue = pspearman::pspearman(S, n, lower.tail=TRUE, approximation="exact")
    }
    pValue = min(2*pValue, 1)
    results = data.frame(pValue)
  }
  return(results)

}



