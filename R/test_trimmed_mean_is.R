#' Independent Samples Trimmed/Yuen Mean Test
#' @description
#' A test to compare two means. The null hypothesis would be that the means of each category are equal in the population.
#' 
#' There are four similar tests, with different assumptions. 
#' 
#' |test|equal variance|normality|
#' |-------|-----------|---------|
#' |Student| yes | yes|
#' |Welch | no | yes|
#' |Trimmed | yes | no | 
#' |Yuen-Welch |no | no |
#' 
#' The Student and Welch are available as separate functions. The Trimmed Means and Yuen-Welch test are available in this one.
#' 
#' @param catField A vector with the categorical data
#' @param scaleField A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' @param dmu Optional difference according to null hypothesis (default is 0)
#' @param trimProp Optional proportion to trim in total for each category. If for example set to 0.1 then 0.05 from each side for each category will be trimmed. Default is 0.1.
#' @param se Optional to indicate which standard error to use. Either "yuen" (default) or "yuen-dixon".
#' 
#' @returns 
#' A dataframe with:
#' \item{n cat. 1}{the sample size of the first category}
#' \item{n cat. 2}{the sample size of the second category}
#' \item{trim mean cat. 1}{the sample trimmed mean of the first category}
#' \item{trim mean cat. 2}{the sample trimmed mean of the second category}
#' \item{diff.}{difference between the two sample means}
#' \item{hyp. diff.}{hypothesized difference between the two population means}
#' \item{statistic}{the test statistic (t-value)}
#' \item{df}{degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' \item{test}{name of test used}
#'  
#' @details
#' **YUEN**
#' 
#' The default `se="yuen"` will perform a Yuen-Welch test.
#' 
#' The formula used is (Yuen, 1974, p. 167):
#' \deqn{t = \frac{\bar{x}_{t,1} - \bar{x}_{t,2}}{SE}}
#' \deqn{sig = 2\times\left(1 - T\left(\left|t\right|, df\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{s_{w,1}^2}{m_1} + \frac{s_{w,2}^2}{m_2}}}
#' \deqn{s_{w,i}^2 = \frac{SSD_{w,i}}{m_i - 1}}
#' \deqn{df = \frac{1}{\frac{c^2}{m_1 - 1} + \frac{\left(1 - c\right)^2}{m_2 -1}}}
#' \deqn{c = \frac{\frac{s_{w,1}^2}{m_1}}{\frac{s_{w,1}^2}{m_1} + \frac{s_{w,2}^2}{m_2}}}
#' \deqn{\bar{x}_{t,i} = \frac{\sum_{j=g_i+1}^{n_i - g_i}y_{i,j}}{}}
#' \deqn{g_i = \lfloor n_i\times p_t\rfloor}
#' \deqn{m_i = n_ - 2\times g_i}
#' \deqn{SSD_{w,i} = g_i\times\left(y_{i,g_i+1} - \bar{x}_{wi}\right)^2 + g_i\times\left(y_{i,n_i-g_i} - \bar{x}_{w,i}\right)^2 + \sum_{j=g+1}^{n_i - g_i} \left(y_{i,j} - \bar{x}_{w,i}\right)^2}
#' \deqn{\bar{x}_{w,i} = \frac{\bar{x}_{t,i}\times m_i + g_i\times\left(y_{i, g_i+1} + y_{i, n_i-g_i}\right)}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{t,i}} the trimmed mean of the scores in category i
#' \item \eqn{x_{w,i}} The Winsorized mean of the scores in category i
#' \item \eqn{SSD_{w,i}} the sum of squared deviations from the Winsorized mean of category i
#' \item \eqn{m_i} the number of scores in the trimmed data set from category i
#' \item \eqn{y_{i,j}} the j-th score after the scores in category i, after they are sorted from low to high
#' \item \eqn{p_t} the proportion of trimming on each side, we can define
#' }
#' 
#' **YUEN-DIXON**
#' 
#' If `se="yuen-dixon` a trimmed means test will be performed.
#' 
#' The formula used is (Yuen & Dixon, 1973, p. 394):
#' \deqn{t = \frac{\bar{x}_{t,1} - \bar{x}_{t,2}}{SE}}
#' \deqn{sig = 2\times\left(1 - T\left(\left|t\right|, df\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{SSD_{w,1} + SSD_{w,2}}{m_1 + m_2 - 2}\times\left(\frac{1}{m_1} + \frac{1}{m_2}\right)}}
#' \deqn{df = m_1 + m_2 - 2}
#' \deqn{\bar{x}_{t,i} = \frac{\sum_{j=g_i+1}^{n_i - g_i}y_{i,j}}{}}
#' \deqn{g_i = \lfloor n_i\times p_t\rfloor}
#' \deqn{m_i = n_ - 2\times g_i}
#' \deqn{SSD_{w,i} = g_i\times\left(y_{i,g_i+1} - \bar{x}_{w,i}\right)^2 + g_i\times\left(y_{i,n_i-g_i} - \bar{x}_{w,i}\right)^2 + \sum_{j=g+1}^{n_i - g_i} \left(y_{i,j} - \bar{x}_{w,i}\right)^2}
#' \deqn{\bar{x}_{w,i} = \frac{\bar{x}_{t,i}\times m_i + g_i\times\left(y_{i, g_i+1} + y_{i, n_i-g_i}\right)}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{ti}} the trimmed mean of the scores in category i
#' \item \eqn{x_{wi}} The Winsorized mean of the scores in category i
#' \item \eqn{SSD_{wi}} the sum of squared deviations from the Winsorized mean of category i
#' \item \eqn{m_i} the number of scores in the trimmed data set from category i
#' \item \eqn{y_{i,j}} the j-th score after the scores in category i, after they are sorted from low to high
#' \item \eqn{p_t} the proportion of trimming on each side, we can define
#' }
#' 
#' @references
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population variances. *Biometrika, 61*(1), 165-170. https://doi.org/10.1093/biomet/61.1.165
#'  
#' Yuen, K. K., & Dixon, W. J. (1973). The approximate behaviour and performance of the two-sample trimmed t. *Biometrika, 60*(2), 369-374. https://doi.org/10.2307/2334550
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_trimmed_mean_is <- function(catField, 
                                 scaleField, 
                                 categories=NULL, 
                                 dmu=0, 
                                 trimProp = 0.1, 
                                 se = c("yuen", "wilcox")){
  
  if (length(se)>1){se="yuen"}
  
  #remove rows with missing values
  df = data.frame(scaleField, catField)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  df$score = as.numeric(df$score)
  
  #the two categories
  if (!is.null(categories)){
    cat1 = categories[1]
    cat2 = categories[2]
  }
  else {
    cat1 = names(table(df$group))[1]
    cat2 = names(table(df$group))[2]
  }
  
  scores1 = df$score[df$group == cat1]
  scores2 = df$score[df$group == cat2]
  
  n1 = length(scores1)
  g1 = floor(n1*trimProp/2)
  n1t = n1 - 2*g1
  m1t = mean(scores1, trim=trimProp/2)
  #note that R base mean(..., trim) will round down
  
  n2 = length(scores2)
  g2 = floor(n2*trimProp/2)
  n2t = n2 - 2*g2
  m2t = mean(scores2, trim=trimProp/2)
  
  #Winsorize the data
  sc1Sort = sort(scores1)
  minReplace1 = sc1Sort[g1+1]
  maxReplace1 = sc1Sort[n1 - g1]
  scores1 = replace(scores1, scores1<minReplace1, minReplace1)
  scores1 = replace(scores1, scores1>maxReplace1, maxReplace1)
  
  sc2Sort = sort(scores2)
  minReplace2 = sc2Sort[g2+1]
  maxReplace2 = sc2Sort[n2 - g2]
  scores2 = replace(scores2, scores2<minReplace2, minReplace2)
  scores2 = replace(scores2, scores2>maxReplace2, maxReplace2)
  var1 = var(scores1)
  var2 = var(scores2)
  
  if (se=="yuen"){
    ssw1 = var1*(n1-1)
    ssw2 = var2*(n2-1)
    var1w = ssw1 / (n1t - 1)
    var2w = ssw2 / (n2t - 1)
    se = (var1w / n1t + var2w / n2t)**0.5
    
    c = (var1w / n1t) / (var1w / n1t + var2w / n2t)
    df = 1 / (c ^ 2 / (n1t - 1) + (1 - c) ^ 2 / (n2t - 1))
    testUsed = "Yuen-Welch independent samples t-test"
  }
  else{
    s2 = ((n1 - 1)*var1 + (n2 - 1)*var2)/((n1t - 1) + (n2t - 1))
    se = sqrt(s2 * (1/n1t + 1/n2t))
    df = n1t + n2t - 2
    testUsed = "Trimmed Mean independent samples t-test"
  }
  
  tValue = (m1t - m2t - dmu) / se
  pValue = 2*(1 - pt(abs(tValue), df))  
  statistic=tValue
  
  
  results <- data.frame(n1, n2, m1t, m2t, m1t - m2t, dmu, statistic, df, pValue, testUsed)
  colnames(results) = c(paste("n", cat1), paste("n", cat2), paste("trim mean", cat1), paste("trim mean", cat2), "diff.", "hyp. diff.", "statistic", "df", "p-value", "test")
  
  return(results)
  
}



