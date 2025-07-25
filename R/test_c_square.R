#' C-square Test
#' @description
#' An improved version of the Brunner-Munzel test. It tests so-called stochastic dominance. It tests if you pick a random score from the first category, will there be a higher (or lower) chance than 50% it will be higher than a random score from the other category. Note this is not the same as a median test. A very short example.
#' 
#' If we have the scores 40, 50 and 60 in group A, and the scores 30, 50, and 51 in group B, the median of each group is 50. However if we pick 40 from group A, there is a 1/3 chance it is higher than a value from group B. If we pick 50 there is also a 1/3 chance, and if we pick 60 there is a 3/3 chance. Together, there is a (1/3+1/3+3/3)/3 = 5/9 chance that if you pick a random value from group A and B, that the one in group A is higher. This is quite off from the 50%.
#' 
#' The test could therefor be used if we have a binary and an ordinal variable. The Brunner-Munzel test is known not to be suitable in small sample sizes. This C-square test has no such limitation.
#'  
#' @param catField A vector or dataframe with the group data
#' @param ordField A vector or dataframe with the scores data
#' @param categories optional list with the two categories to use from catField. If not set the first two found will be used
#' @param levels optional list with the scores in order
#' 
#' @returns
#' A dataframe with:
#' \item{var. est.}{the estimated overall variance}
#' \item{test-statistic}{the Cliff-Delta statistic}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the significance (p-value), two-tailed}
#' \item{categories}{the two categories that were used}
#' 
#' @details
#' The test statistic (Schuurhuis et al., 2025, p. 9, eq. 17):
#' \deqn{C^2 = \frac{4}{\hat{\sigma}_N^2} \times \hat{\theta}_N\times\left(1 - \hat{\theta}_N\right) \times \left(\hat{\theta}_N - \frac{1}{2}\right)^2}
#' 
#' the estimate of the overall variance of (Schuurhuis et al., 2025, p. 5, eq. 5):
#' \deqn{\hat{\sigma}_N^2 = \frac{1}{d_n}\times\left(\left(\sum_{k=1}^2 SS_k^*\right) - n_1 \times n_2 \times \left(\hat{\theta}_N\times\left(1 - \hat{\theta}_N\right) - \frac{\hat{\tau}_N}{4}\right)\right)}
#' 
#' the sum of squares for each of the two categories of the placement values:
#' \deqn{SS_k^\ast = \sum_{i=1}^{n_k} \left(R_{ik}^\ast - \bar{R}_{k}^*\right)^2}
#' 
#' the mean of placement values (Schuurhuis et al., 2025, p. 5):
#' \deqn{\bar{R}_{k}^\ast = \frac{\sum_{i=1}^{n_k} R_{ik}^\ast}{n_k}}
#' 
#' the placement values:
#' \deqn{R_{ik}^* = R_{ik} - R_{ik}^{(k)}}
#' 
#' Denominator (Schuurhuis et al., 2025, p. 5):
#' \deqn{d_n = n_1\times\left(n_1 - 1\right) \times n_2\times\left(n_2 - 1\right)}
#' 
#' the probability of ties in the overlap (Schuurhuis et al., 2025, p. 5, eq. 4):
#' \deqn{\hat{\tau}_N = \frac{1}{n_1}\times\left(\bar{R}_2^+ - \bar{R}_2^- - \left(\bar{R}_2^{(2)+} - \bar{R}_2^{(2)-}\right)\right)}
#' 
#' estimated Mann-Whitney effect (Schuurhuis et al., 2025, p. 4, eq. 1):
#' \deqn{\hat{\theta} = \frac{1}{n_1} \times \left(\bar{R}_2 - \frac{n_2 + 1}{2}\right)}
#' note that this is the same as \eqn{\hat{p}} in the Brunner-Munzel test.
#' 
#' some means from the second category:
#' \deqn{\bar{R}_2 = \frac{\sum_{i=1}^{n_2} R_{i2}}{n_2}, \bar{R}_2^+ = \frac{\sum_{i=1}^{n_2} R_{i2}^+}{n_2}, \bar{R}_2^{(2)+} = \frac{\sum_{i=1}^{n_2} R_{i2}^{(2)+}}{n_2}, \bar{R}_2^{(2)-} = \frac{\sum_{i=1}^{n_2} R_{i2}^{(2)-}}{n_2}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{R_{ik}}, the mid rank of the i-th score in category k, when using all combined scores
#' \item \eqn{R_{ik}^{(k)}}, the mid-rank of the i-th score in category k, when using only scores from category k
#' \item \eqn{R_{ik}^-}, the minimum rank of the i-th score in category k, when using all combined scores
#' \item \eqn{R_{ik}^{(k)-}}, the minimum-rank of the i-th score in category k, when using only scores from category k
#' \item \eqn{R_{ik}^+}, the maximum rank of the i-th score in category k, when using all combined scores
#' \item \eqn{R_{ik}^{(k)+}}, the maximum-rank of the i-th score in category k, when using only scores from category k
#' \item \eqn{N}, the total sample size
#' \item \eqn{n_{k}}, the number of scores in category k
#' }
#' 
#' @section Before, After and Alternatives:
#' Before running the test you might first want to get an impression using a cross table using \code{\link{tab_cross}}, or a stacked bar chart using \code{\link{vi_bar_stacked_multiple}} for a visualisation.
#' 
#' After the test you might want an effect size measure, either \code{\link{es_common_language_is}} for the CLES, \code{\link{es_hodges_lehmann_is}} for Hodges-Lehmann, or \code{\link{r_rank_biserial_is}} for the (Glass) rank biserial (Cliff delta).
#' 
#' Independent samples tests for a binary vs ordinal variable include Brunner-Munzel (\code{\link{ts_brunner_munzel}}), C-square (\code{\link{ts_c_square}}), Cliff-Delta (\code{\link{ts_cliff_delta_is}}), Fligner-Policello (\code{\link{ts_fligner_policello}}), Mann-Whitney U (\code{\link{ts_mann_whitney}}), Mood-Median (\code{\link{ts_mood_median}})
#' 
#' @references 
#' Schuurhuis, S., Konietschke, F., & Brunner, E. (2025). A new approach to the Nonparametric Behrens-Fisher problem with compatible confidence intervals (arXiv:2504.01796). arXiv. https://doi.org/10.48550/arXiv.2504.01796
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv"
#' dfr <- read.csv(dataFile, sep=";", na.strings=c("", "NA"))
#' coding = c("Far too little", "too little", "Enough", "Too much", "Far too much")
#' ts_c_square(dfr[['Gen_Gender']], dfr[['Mix_NrAct']], levels=coding)
#' 
#' @export
ts_c_square <- function(catField, ordField, categories=NULL, levels=NULL){
  
  #remove rows with missing values
  dfr = data.frame(ordField, catField)
  dfr = na.omit(dfr)
  colnames(dfr) = c("score", "group")
  
  #replace the ordinal values if levels is provided
  if (!is.null(levels)){
    dfr$score = factor(dfr$score, ordered = TRUE, levels = levels)        
  }
  dfr$score = as.numeric(dfr$score)
  
  #the two categories
  if (!is.null(categories)){
    cat1 = categories[1]
    cat2 = categories[2]
  }
  else {
    cat1 = names(table(dfr[ ,2]))[1]
    cat2 = names(table(dfr[ ,2]))[2]
  }
  dfr = na.omit(dfr)
  #seperate the scores for each category
  X = unname(unlist((subset(dfr, dfr[ ,2] == cat1)[1])))
  Y = unname(unlist((subset(dfr, dfr[ ,2] == cat2)[1])))
  
  n1 = length(X)
  n2 = length(Y)    
  N = n1 + n2
  minN = min(n1, n2)
  
  #combine this into one long list
  A = c(X, Y)
  
  
  # all ranks for within category 1
  Ri1 = rank(X)
  Ri1Min = rank(X, ties.method='min')
  Ri1Plus = rank(X, ties.method='max')
  # all ranks for within category 2
  Ri2 = rank(Y)
  Ri2Min = rank(Y, ties.method='min')
  Ri2Plus = rank(Y, ties.method='max')
  
  #pooled ranks
  Rik = rank(A)
  RikMin = rank(A, ties.method='min')
  RikPlus = rank(A, ties.method='max')
  
  #placement values
  Rast1 =Rik[1:n1] - Ri1
  Rast2 = Rik[(n1+1):N] - Ri2
  
  #means
  barR2 = sum(Rik[(n1+1):N])/n2
  barR2plus = sum(RikPlus[(n1+1):N])/n2
  barR2min = sum(RikMin[(n1+1):N])/n2
  barR2plus2 = mean(Ri2Plus)
  barR2min2 = mean(Ri2Min)
  barRast1 = mean(Rast1)
  barRast2 = mean(Rast2)
  
  SSrAst1 = sum((Rast1 - barRast1)**2)
  SSrAst2 = sum((Rast2 - barRast2)**2)
  
  # test
  thetaN = 1/n1*(barR2 - (n2 + 1)/2)
  tauN = 1/n1*(barR2plus - barR2min - (barR2plus2 - barR2min2))
  dn = n1*(n1 - 1)*n2*(n2 - 1)
  varN = 1/dn*(SSrAst1 + SSrAst2 - n1*n2*(thetaN*(1 - thetaN) - tauN/4))
  Csq = 1/varN*4*thetaN*(1 - thetaN)*(thetaN - 1/2)**2
  pVal = 1 - pchisq(Csq, 1)
  
  #results
  cats = paste0(cat1, ', ', cat2)
  
  results <- data.frame(varN, Csq, 1, pVal, cats)
  colnames(results) <- c("var. est.", "test statistic", "df", "p-value", "categories")
  
  return(results)  
}



