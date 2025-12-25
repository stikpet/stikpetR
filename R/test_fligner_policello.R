#' Fligner-Policello Test
#' @description
#' An alternative for the more famous Mann-Whitney U test. The MWU test has as an assumption that the scores in the two categories have the same shape and have unequal variances (Fong & Huang, 2019). The Fligner-Policello test does not, although the distribution around their medians should be symmetric in the population Zaiontz (n.d.).
#' 
#' Roughly put the assumption for this test is that the two categories have the same median in the population.
#' 
#' The function is shown in this [YouTube video](https://youtu.be/Ek36vYbzl9s) and the test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tests/FlignerPolicello.html)
#' 
#' @param catField A vector with the scores data
#' @param ordField A vector with the group data
#' @param categories optional vector with categories to use and order for the categorical field. Otherwise the first two found will be used.
#' @param levels optional vector with the labels of the ordinal field in order.
#' @param ties boolean to indicate the use of a ties correction. Default is TRUE
#' @param cc boolean to indicate the use of a continuity correction. Default is FALSE
#' 
#' @returns
#' A dataframe with:
#' \item{n}{the sample size}
#' \item{statistic}{test statistic}
#' \item{p-value}{significance (p-value)}
#' \item{test}{description of the test used}
#' 
#' @details
#' The formula used is:
#' \deqn{z = \frac{N_Y - N_X}{2\times\sqrt{SS_X + SS_Y - M_X\times M_Y}}}
#' With:
#' \deqn{SS_X = \sum_{x\in X} \left(N_X - M_X\right)^2, SS_Y = \sum_{y\in Y} \left(N_Y - M_Y\right)^2}
#' \deqn{M_X = \frac{N_X}{n_x}, M_Y = \frac{N_Y}{n_y}}
#' \deqn{N_X = \sum_{x \in X} N\left(x\right), N_Y = \sum_{y \in Y} N\left(y\right)}
#' \deqn{N\left(y\right) = \sum_{x\in X} f\left(y, x\right)}
#' \deqn{N\left(x\right) = \sum_{y\in Y} f\left(x, y\right)}
#' \deqn{f\left(a, b\right) = \begin{cases} 1 & \text{ if } a> b \\ 0 & \text{ if } a\leq b \end{cases}}
#' 
#' In case of a tie correction (Hollander et al., 2014, p. 146):
#' \deqn{z = \frac{\left|N_Y - N_X\right| - 0.5}{2\times\sqrt{SS_X + SS_Y - M_X\times M_Y}}}
#' \deqn{f\left(a, b\right) = \begin{cases} 1 & \text{ if } a> b \\ 0.5 & \text{ if } a = b \\ 0 & \text{ if } a\leq b \end{cases}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{X} the scores in the first category
#' \item \eqn{Y} the scores in the second category
#' \item \eqn{n_i} the number of scores in the i category
#' }
#' 
#' The test is described by Fligner and Policello (1981), and can also be found in Kloke and McKean (2015, p. 68)
#' 
#' @section Before, After and Alternatives:
#' Before running the test you might first want to get an impression using a cross table using \code{\link{tab_cross}}, or a stacked bar chart using \code{\link{vi_bar_stacked_multiple}} for a visualisation.
#' 
#' After the test you might want an effect size measure, either \code{\link{es_common_language_is}} for the CLES, \code{\link{es_hodges_lehmann_is}} for Hodges-Lehmann, or \code{\link{r_rank_biserial_is}} for the (Glass) rank biserial (Cliff delta).
#' 
#' Independent samples tests for a binary vs ordinal variable include Brunner-Munzel (\code{\link{ts_brunner_munzel}}), C-square (\code{\link{ts_c_square}}), Cliff-Delta (\code{\link{ts_cliff_delta_is}}), Fligner-Policello (\code{\link{ts_fligner_policello}}), Mann-Whitney U (\code{\link{ts_mann_whitney}}), Mood-Median (\code{\link{ts_mood_median}})
#' 
#' @references 
#' Fligner, M. A., & Policello, G. E. (1981). Robust rank procedures for the Behrens-Fisher problem. *Journal of the American Statistical Association, 76*(373), 162-168. https://doi.org/10.1080/01621459.1981.10477623
#' 
#' Hollander, M., Wolfe, D. A., & Chicken, E. (2014). *Nonparametric statistical methods* (3rd ed.). John Wiley & Sons, Inc.
#' 
#' Kloke, J., & McKean, J. W. (2015). *Nonparametric statistical methods using R*. CRC Press, Taylor & Francis.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' myLevels = c('Not scientific at all', 'Not too scientific', 'Pretty scientific', 'Very scientific')
#' ts_fligner_policello(df1[['sex']], df1[['accntsci']], levels = myLevels)
#' ts_fligner_policello(df1[['sex']], df1[['accntsci']], levels = myLevels, ties= FALSE, cc=TRUE)
#' ts_fligner_policello(df1[['sex']], df1[['accntsci']], levels = myLevels, ties= TRUE, cc=FALSE)
#' 
#' #Example 2: vectors
#' binary = c("apple", "apple", "apple", "peer", "peer", "peer", "peer")
#' ordinal = c(4, 3, 1, 6, 5, 7, 2)
#' ts_fligner_policello(binary, ordinal, categories=c("peer", "apple"))
#' 
#' @export
ts_fligner_policello <- function(catField, 
                                 ordField, 
                                 categories=NULL, 
                                 levels=NULL, 
                                 ties= TRUE, 
                                 cc=FALSE){
  testUsed = "Fligner-Policello test"
  
  #create a cross table
  ct = tab_cross(ordField, catField, order1 = levels, order2 = categories)
  nr = nrow(ct)
  
  fx = c()
  fy = c()
  
  if (ties) {
    fx[1] = 0.5*ct[1,2]
    fy[1] = 0.5*ct[1,1]
  }
  else{
    fx[1] = 0
    fy[1] = 0
  }
  
  nx = fx[1] * ct[1,1]
  ny = fy[1] * ct[1,2]
  
  n1 = ct[1,1]
  n2 = ct[1,2]
  
  for (i in 2:nr){
    fx[i] = fx[i-1] + ct[i-1,2]
    fy[i] = fy[i - 1] + ct[i - 1, 1]
    
    if (ties){
      fx[i] = fx[i] + 0.5 * ct[i, 2] - 0.5 * ct[i-1,2]
      fy[i] = fy[i] + 0.5 * ct[i, 1] - 0.5 * ct[i - 1, 1]
    }
    
    nx = nx + fx[i] * ct[i, 1]
    n1 = n1 + ct[i, 1]
    
    ny = ny + fy[i] * ct[i, 2] 
    n2 = n2 + ct[i, 2] 
  }
  
  MX = nx / n1
  MY = ny / n2
  
  ssx = ssy = 0
  for (i in 1:nr){
    ssx = ssx + ct[i, 1] * (fx[i] - MX) ^ 2
    ssy = ssy + ct[i, 2]  * (fy[i] - MY) ^ 2
  }
  
  if (cc){num = abs(nx - ny) - 0.5}
  else{num = nx - ny}
  
  z = num / (2 * sqrt(ssx + ssy + MX * MY))
  
  pValue = 2*(1 - pnorm(abs(z)))
  
  if (cc && ties){testUsed = "Fligner-Policello test, with continuity and ties correction"}
  else if (cc){testUsed = "Fligner-Policello test, with continuity correction"}
  else if (ties){testUsed = "Fligner-Policello test, with ties correction"}
  else{testUsed = "Fligner-Policello test"}
  
  results <- data.frame(n1+n2, z, pValue, testUsed)
  colnames(results) = c("n", "statistic", "p-value", "test")
  
  return (results)
  
}



