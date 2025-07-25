#' Cliff Delta Test
#' @description
#' The Cliff test, is a test for stochastic-equivelance. This means that even if the medians are equal between two independent samples, this test could be significant.
#' 
#' Lets say we have one group A that scored 1, 2, 2, 5, 6, 6, 7, and another group B that scored 4, 4, 4, 5, 10, 10, 12. Each group has the same median (i.e. 5), and are symmetric around the median, but if a high score is positive, most people would rather be in group B than in group A. This is where 'stochastic equality' comes in. It looks at the chance if you pick a random person from group A and B each, the one from group A scores lower than the one from group B, and add half the chance that their equal. In this example that's about 0.68.
#' 
#' Cliff (1993) used a standard normal distribution, while Vargha and Delaney later used a t-distribution (Vargha, 2000; Vargha & Delaney, 2000, Delaney & Vargha, 2002).
#' 
#' @param catField A vector or dataframe with the group data
#' @param ordField A vector or dataframe with the scores data
#' @param categories optional list with the two categories to use from catField. If not set the first two found will be used
#' @param levels optional list with the scores in order
#' @param var_ver c("unbiased", "consistent"), optional version of estimated variance to use.
#' @param test_ver c('cliff', 'fligner-policello', 'brunner-munzel', 'fp', 'bm'), optional version of the test to use.
#' @param round_df bool, optional round degrees of freedom
#' 
#' @returns
#' A dataframe with:
#' \item{var}{the estimated overall variance}
#' \item{test-statistic}{the Cliff-Delta statistic}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the significance (p-value), two-tailed}
#' \item{categories}{the two categories that were used}
#' \item{var used}{the variance version used}
#' \item{test used}{description of the test used and distribution}
#' 
#' @details
#' The test-statistic for this test is (Cliff, 1993, p. 500):
#' \deqn{C = \frac{d}{\sqrt{\hat{\sigma}_d^2}}}
#' 
#' Which follows a normal distribution (Cliff, 1993, p. 500), i.e. 
#' \deqn{p = 2 \times \left(1 - \Phi\left(\left|C\right|\right)\right)}
#' 
#' With:
#' \deqn{d = \frac{\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} d_{x_i, y_j}}{n_1 \times n_2}}
#' \deqn{d_{x_i, y_j} = S_{i,j}}
#' \deqn{S_{i,j} = \text{sign}\left(x_i - y_j\right) = \begin{cases} 1 & x_i >t y_j \\ 0 & x_i = y_j \\ -1 & x_i < y_j \end{cases}}
#' \deqn{\hat{\sigma}_d^2 = \max{\left(s_d^2, s_{d, min}^2\right)}}
#' \deqn{s_{d, min}^2 = \frac{1 - d^2}{n_1\times n_2 - 1}}
#' \deqn{s_d^2 = \frac{n_2^2 \times SS_{\hat{d}_{1}} + n_1^2 \times SS_{\hat{d}_{2}} - SS_{d_{xy}}}{n_1\times n_2 \times \left(n_1 - 1\right) \times \left(n_2 - 1\right)}}
#' \deqn{SS_{\hat{d}_{1}} = \sum_{i=1}^{n_1} \left(\hat{d}_{i1} - d\right)^2, SS_{\hat{d}_{2}} = \sum_{j=1}^{n_2} \left(\hat{d}_{j2} - d\right)^2, SS_{d_{xy}} = \sum_{i=1}^{n_1}\sum_{j=1}^{n_2} \left(d_{x_i, y_j} - d\right)^2}
#' \deqn{\hat{d}_{i1} = \frac{1}{n_2} \times \sum_{j=1}^{n_2} \begin{cases} 1 & x_i > y_j \\ 0 & x_i = y_j \\ -1 & x_i < y_j \end{cases} = \frac{1}{n_2} \times \sum_{j=1}^{n_2} S_{i,j}}
#' \deqn{\hat{d}_{j2} = \frac{1}{n_1} \times \sum_{i=1}^{n_1} \begin{cases} 1 & x_i > y_j \\ 0 & x_i = y_j \\ -1 & x_i < y_j \end{cases} = \frac{1}{n_1} \times \sum_{i=1}^{n_1} S_{i,j}}
#' 
#' Alternatively, but with the same result, the sample variance of d, can be calculated with:
#' \deqn{s_d^2 = \frac{n_2}{n_1} \times \frac{s_{\hat{d}_1}^2}{n_2 - 1} + \frac{n_1}{n_2} \times \frac{s_{\hat{d}_2}^2}{n_1 - 1} - \frac{s_{d_{xy}}^2}{n_1 \times n_2}}
#' \deqn{s_{\hat{d}_1}^2 = \frac{SS_{\hat{d}_{1}}}{n_1 -1}, s_{\hat{d}_2}^2 = \frac{SS_{\hat{d}_{2}}}{n_2 -1}, s_{d_{xy}}^2 = \frac{SS_{d_{xy}}}{\left(n_1 - 1\right) \times \left(n_2 - 1\right)}}
#' 
#' A different estimate (a 'consistent') is given by (Cliff, 1993, p. 499, eq. 7):
#' \deqn{\hat{\sigma}_d^2 = \frac{\left(n_2 - 1\right) \times s_{\hat{d}_1}^2 + \left(n_1 - 1\right) \times s_{\hat{d}_2}^2 + s_{d_{xy}}^2}{n_1 \times n_2}}
#' 
#' Vargha (2000, p. 280) and also Vargha and Delaney (2000, p. 7, eq. 9) use a t-distribution, instead of the standard normal distribution. They use the same test-statistic, and as degrees of freedom they use:
#' \deqn{p = 2 \times \left(1 - t\left(\left|C\right|, df\right)\right)}
#' \deqn{df = \frac{\left(a + b\right)^2}{\frac{a^2}{n_1 - 1} + \frac{b^2}{n_2 - 1}}}
#' 
#' With:
#' \deqn{a_{BM} = \frac{1}{n_1} \times \frac{s_{R_1^\ast}^2}{n_2^2}, b_{BM} = \frac{1}{n_2} \times \frac{s_{R_2^*}^2}{n_1^2}}
#' \deqn{s_{R_1^\ast}^2 = \frac{SS_{R_1^\ast}}{n_1 -1}, s_{R_2^\ast}^2 = \frac{SS_{R_2^*}}{n_2 -1}}
#' \deqn{SS_{R_1^\ast} = \sum_{i=1}^{n_1} \left(r_{i1}^\ast - \bar{R}_1^\ast\right)^2, SS_{R_2^\ast} = \sum_{j=1}^{n_2} \left(r_{j2}^\ast - \bar{R}_2^\ast\right)^2}
#' \deqn{\bar{R}_{1}^\ast = \frac{\sum_{i=1}^{n_1} r_{i1}^\ast}{n_1}, \bar{R}_2^\ast = \frac{\sum_{j=1}^{n_2} r_{j2}^\ast}{n_2}}
#' \deqn{r_{i,1}^\ast = \sum_{j=1}^{n_2} \begin{cases} 1 & s_{i,j} = 1 \\  0.5 & s_{i,j} = 0 \\  0 & s_{i,j} = -1 \end{cases}}
#' \deqn{r_{i,2}^\ast = \sum_{i=1}^{n_1} \begin{cases} 1 & s_{i,j} = -1 \\  0.5 & s_{i,j} = 0 \\  0 & s_{i,j} = 1 \end{cases}}
#' 
#' This is in line with the Fligner-Policello test, and used when *test_used='fligner-policello*. 
#' 
#' Delaney and Vargha (2002) proposed also an alternative degrees of freedom, in line with the Brunner-Munzel test, as:
#' \deqn{a_{FP} = \frac{s_{R_1^\ast}^2}{n_1}, b_{FP} = \frac{s_{R_2^*}^2}{n_2}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_{k}}, the number of scores in category k
#' \item \eqn{\Phi\left(\dots\right)}, the cumulative distribution function of the standard normal distribution
#' \item \eqn{t\left(\dots\right)}, the cumulative distribution function of the t distribution
#' }
#' 
#' A few additional notes.
#' In Vargha and Delaney (2000, p. 7, eq. 9) they mention "df is rounded to the nearest integer", which is why this function allows for a rounding of the degrees of freedom (df).
#' 
#' The \eqn{r^*} values, are so-called placement values. These can also be calculated by subtracting the within-mid-rank from the pooled mid-rank.
#' 
#' The \eqn{d} value is known as Cliff Delta, and is also the same as the (Glass) rank biserial correlation coefficient. This in itself is sometimes used as an effect size measure, and various methods are available to calculate it.
#' 
#' For the \eqn{\hat{d}_{j2}} Cliff (1993) mentions: " \eqn{ d_{.j} } represents the proportion of scores from the first population that lies **above** a given score from the second, minus the reverse " (p. 499), which is the one used here. However, Delaney and Vargha (2002) wrote: " \eqn{d_{.j}} denotes the proportion of the \eqn{X} scores that lie **below** \eqn{Y_j} minus the proportion that lie above" (p. 9), but this is most likely a mistake.
#' 
#' @section Before, After and Alternatives:
#' Before running the test you might first want to get an impression using a cross table using \code{\link{tab_cross}}, or a stacked bar chart using \code{\link{vi_bar_stacked_multiple}} for a visualisation.
#' 
#' After the test you might want an effect size measure, either \code{\link{es_common_language_is}} for the CLES, \code{\link{es_hodges_lehmann_is}} for Hodges-Lehmann, or \code{\link{r_rank_biserial_is}} for the (Glass) rank biserial (Cliff delta).
#' 
#' Independent samples tests for a binary vs ordinal variable include Brunner-Munzel (\code{\link{ts_brunner_munzel}}), C-square (\code{\link{ts_c_square}}), Cliff-Delta (\code{\link{ts_cliff_delta_is}}), Fligner-Policello (\code{\link{ts_fligner_policello}}), Mann-Whitney U (\code{\link{ts_mann_whitney}}), Mood-Median (\code{\link{ts_mood_median}})
#' 
#' 
#' @references 
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions. *Psychological Bulletin, 114*(3), 494-509. doi:10.1037/0033-2909.114.3.494
#' 
#' Delaney, H. D., & Vargha, A. (2002). Comparing several robust tests of stochastic equality with ordinally scaled variables and small to moderate sized samples. *Psychological Methods, 7*(4), 485-503. doi:10.1037/1082-989X.7.4.485
#' 
#' Vargha A. (2000). Ket pszichologiai populacio sztochasztikus egyenlosegenek ellenorzesere alkalmas statisztikai probak osszehasonlito vizsgalata. *Magyar Pszichologiai Szemle, 55*(2-3), Article 2-3. doi:10/1/mpszle.55.2000.2-3.5.pdf
#' 
#' Vargha, A., & Delaney, H. D. (2000). Comparing several robust tests of stochastic equality. https://eric.ed.gov/?id=ED441836
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv"
#' dfr <- read.csv(dataFile, sep=";", na.strings=c("", "NA"))
#' coding = c("Far too little", "too little", "Enough", "Too much", "Far too much")
#' ts_cliff_delta_is(dfr[['Gen_Gender']], dfr[['Mix_NrAct']], levels=coding)
#' 
#' @export
ts_cliff_delta_is <- function(catField, ordField, categories=NULL, levels=NULL, var_ver='unbiased', test_ver='cliff', round_df=FALSE){
  
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
  
  diff <- outer(X, Y, "-")
  rownames(diff) <- names(X)
  colnames(diff) <- names(Y)
  S <- sign(diff)
  di1 <- rowSums(S) / ncol(S)
  dj2 <- colSums(S) / nrow(S)
  d = mean(di1)
  
  SSd1 = sum((di1 - d)**2)
  SSd2 = sum((dj2 - d)**2)
  SSd = sum((S - d)^2)
  
  if (var_ver=='unbiased'){
    var_used = 'unbiased'
    var_d_v1 = (n2^2 * SSd1 + n1^2 * SSd2 - SSd) / (n1*n2*(n1 - 1)*(n2 - 1))
    var_d_min = (1 - d^2)/(n1*n2 - 1)
    var = max(var_d_v1, var_d_min)}
  else if (var_ver=='consistent'){
    var_used = 'consistent'
    var_d1 = SSd1 / (n1 - 1)
    var_d2 = SSd2 / (n2 - 1)
    var_d = SSd / ((n1 - 1)*(n2 - 1))
    var = ((n2 - 1)*var_d1 + (n1 - 1)*var_d2 + var_d)/(n1*n2)}
  
  C = d/(var**0.5)    
  
  if (test_ver=='cliff'){
    df = 'n.a.'
    test_used = 'Cliff Delta test with standard normal distribution'
    
    p = 2 * (1 - pnorm(abs(C))) 
  }
  else{                
    r_ast_1 <- rowSums(S == 1) + 0.5 * rowSums(S == 0)
    r_ast_2 <- colSums(S == -1) + 0.5 * colSums(S == 0)
    
    var_R_ast_1 = var(r_ast_1)
    var_R_ast_2 = var(r_ast_2)
    
    if (test_ver=='fligner-policello' || test_ver=='fp'){
      test_used = 'Cliff Delta test with t distribution, and Fligner-Policello df'
      a = var_R_ast_1/n1
      b = var_R_ast_2/n2}
    else if (test_ver=='brunner-munzel' || test_ver=='bm'){
      test_used = 'Cliff Delta test with t distribution, and Brunner-Munzel df'
      a = 1/n1*var_R_ast_1/n2^2
      b = 1/n2*var_R_ast_2/n1^2}
    
    # degrees of freedom
    df = (a + b)^2 /(a^2/(n1 - 1) + b^2/(n2 - 1))
    if (round_df){
      df = round(df,0)}
    
    p = 2*(1 - pt(abs(C), df)) 
  }
  cats = paste0(cat1, ', ', cat2)
  
  results <- data.frame(var, C, df, p, cats, var_used, test_used)
  colnames(results) <- c("var", "test statistic", "df", "p-value", "categories", 'var used', "test used")
  
  return(results)  
}



