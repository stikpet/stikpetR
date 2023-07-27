#' Welch t Test (Independent Samples)
#' @description
#' A test to compare two means. The null hypothesis would be that the means of each category are equal in the population.
#' 
#' Unlike the Student t-test, the Welch test does not assume the variances of the two categories to be equal in the population. Ruxten (2006) even argues that the Welch t-test should always be prefered over the Student t-test.
#' 
#' There are four similar tests, with different assumptions.
#' |test|equal variance|normality|
#' |-------|-----------|---------|
#' |Student| yes | yes|
#' |Welch | no | yes|
#' |Trimmed | yes | no | 
#' |Yuen-Welch |no | no |
#' 
#' The Trimmed and Yuen-Welch can be found in the **ts_trimmed_mean_is()**, and the Student t-test with the **ts_student_t_is().**
#' 
#' @param catField A vector with the categorical data
#' @param scaleField A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' @param dmu Optional difference according to null hypothesis (default is 0)
#' 
#' @returns 
#' A dataframe with:
#' \item{n cat. 1}{the sample size of the first category}
#' \item{n cat. 2}{the sample size of the second category}
#' \item{mean cat. 1}{the sample mean of the first category}
#' \item{mean cat. 2}{the sample mean of the second category}
#' \item{diff.}{difference between the two sample means}
#' \item{hyp. diff.}{hypothesized difference between the two population means}
#' \item{statistic}{the test statistic (t-value)}
#' \item{df}{degrees of freedom}
#' \item{pValue}{the significance (p-value)}
#' \item{test}{name of test used}
#' 
#' 
#' @details
#' 
#' The formula used is:
#' \deqn{t = \frac{\bar{x}_1 - \bar{x}_2}{SE}}
#' \deqn{df = \frac{SE^4}{\frac{\left(s_1^2\right)^2}{n_1^2\times\left(n_1 - 1\right)} + \frac{\left(s_2^2\right)^2}{n_2^2\times\left(n_2 - 1\right)}}}
#' \deqn{sig. = 2\times\left(1 - T\left(\left|t\right|, df\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' }
#' 
#' @references
#' Ruxton, G. D. (2006). The unequal variance t-test is an underused alternative to Student’s t-test and the Mann–Whitney U test. *Behavioral Ecology, 17*(4), 688–690. https://doi.org/10.1093/beheco/ark016
#' 
#' Welch, B. L. (1947). The generalization of `Student’s’ problem when several different population variances are involved. *Biometrika, 34*(1/2), 28–35. https://doi.org/10.2307/2332510
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
#' ts_welch_t_is(df1['sex'], ex1)
#' 
#' #Example 2: vectors
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("nat.","int.","int.","nat.","int.", "int.","nat.","nat.","int.","int.","int.","int.","int.","int.","nat.", "int." ,NA,"nat.","int.","int.")
#' ts_welch_t_is(groups, scores)
#' 
#' @export
ts_welch_t_is <- function(catField, scaleField, categories=NULL, dmu=0){
  
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
  
  var1 = var(X1)
  var2 = var(X2)
  
  sse = var1/n1 + var2/n2
  se = sqrt(sse)
  
  m1 = mean(X1)
  m2 = mean(X2)
  
  t = (m1 - m2 - dmu)/se
  
  df = sse**2/(var1**2/(n1**2*(n1 - 1)) + var2**2/(n2**2*(n2 - 1)))
  
  pValue = 2*(1-pt(abs(t), df))
  
  statistic = t
  testUsed = "Welch independent samples t-test"
  
  results <- data.frame(n1, n2, m1, m2, m1 - m2, dmu, statistic, df, pValue, testUsed)
  colnames(results) = c(paste("n", cat1), paste("n", cat2), paste("mean", cat1), paste("mean", cat2), "diff.", "hyp. diff.", "statistic", "df", "p-value", "test")
  
  return(results)
  
}