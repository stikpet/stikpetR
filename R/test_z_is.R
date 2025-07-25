#' Independent Samples Z Test
#' @description
#' A test to compare two means. It requires the population variances, but if these are unknown for large enough sample sizes, the sample variances can be used instead.
#' 
#' For smaller sample sizes a t-test (Student, Welch or Trimmed Means) could be used instead.
#' 
#' @param catField A vector with the categorical data
#' @param scaleField A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' @param dmu Optional difference according to null hypothesis (default is 0)
#' @param sigma1 Optional population standard deviation of the first group, if NULL sample results will be used
#' @param sigma2 Optional population standard deviation of the second group, if NULL sample results will be used
#' 
#' @returns 
#' A dataframe with:
#' \item{n cat. 1}{the sample size of the first category}
#' \item{n cat. 2}{the sample size of the second category}
#' \item{mean cat. 1}{the sample mean of the first category}
#' \item{mean cat. 2}{the sample mean of the second category}
#' \item{diff.}{difference between the two sample means}
#' \item{hyp. diff.}{hypothesized difference between the two population means}
#' \item{statistic}{the test statistic (z-value)}
#' \item{pValue}{the significance (p-value)}
#' \item{test}{name of test used}
#' 
#' @details
#' 
#' The formula used is:
#' \deqn{z = \frac{\bar{x}_1 - \bar{x}_2}{SE}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}}}
#' \deqn{\sigma_i^2 \approx s_i^2 = \frac{\sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' }
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['age']
#' ex1 = replace(ex1, ex1=="89 OR OLDER", "90")
#' ts_z_is(df1['sex'], ex1)
#' 
#' #Example 2: vectors
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("nat.","int.","int.","nat.","int.", "int.","nat.","nat.","int.",
#' "int.","int.","int.","int.","int.","nat.", "int." ,NA,"nat.","int.","int.")
#' ts_z_is(groups, scores)
#' 
#' @export
ts_z_is <- function(catField, scaleField, categories=NULL, dmu=0, sigma1=NULL, sigma2=NULL){
  
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
  
  X1 = df$score[df$group == cat1]
  X2 = df$score[df$group == cat2]
  
  n1 = length(X1)
  n2 = length(X2)
  
  if (is.null(sigma1)) {
    var1 = var(X1)
  }
  else{
    var1 = sigma1**2
  }
  
  if (is.null(sigma2)) {
    var2 = var(X2)
  }
  else{
    var2 = sigma2**2
  }
  
  sse = var1/n1 + var2/n2
  se = sqrt(sse)
  
  m1 = mean(X1)
  m2 = mean(X2)
  
  z = (m1 - m2 - dmu)/se
  pValue = 2*(1-pnorm(abs(z)))
  statistic = z
  testUsed = "independent samples z-test"
  
  results <- data.frame(n1, n2, m1, m2, m1 - m2, dmu, statistic, pValue, testUsed)
  colnames(results) = c(paste("n", cat1), paste("n", cat2), paste("mean", cat1), paste("mean", cat2), "diff.", "hyp. diff.", "statistic", "p-value", "test")
  
  return(results)
  
}



