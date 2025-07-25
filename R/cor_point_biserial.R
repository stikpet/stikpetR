#' Point Biserial Correlation Coefficient
#' @description
#' This can be seen as coding a binary variable with the groups into 0 and 1, and then calculates a (Pearson) correlation coefficient between the those values and the scores (Tate, 1954, p. 603).
#' 
#' As the name implies a correlation coefficient indicates how two variables co-relate, i.e. if one goes up is it likely for the other to go up or down. A zero would indicate there is not (linear) relation, while a -1 would mean a perfect negative correlation (if one goes up, the other goes down, and vice versa), and a +1 a perfect positive correlation (if one goes up, the other also goes up, and vice versa).
#' 
#' With two categories we could read this more as if the score go up and there is a positive correlation, it is more likely that it came from a category 1 case, rather than a category 0.
#' 
#' Note that if the two categories come from a so-called latent normally distributed variable, the *biserial correlation* might be better. This is the case if scores were categorized and then compared to some other numeric scores (e.g. grades being categorized into pass/fail, and then use this pass/fail to correlate with age). A separate function is available on the biserial correlation.
#' 
#' @param catField A vector with the categorical data
#' @param scaleField A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' 
#' @returns 
#' A dataframe with:
#' \item{cat. 1}{the category that was used as category 1}
#' \item{cat. 2}{the category that was used as category 2}
#' \item{mean 1}{the arithmetic mean of the scores from category 1}
#' \item{mean 2}{the arithmetic mean of the scores from category 2}
#' \item{r_pb}{the point-biserial correlation coefficient}
#' 
#' @details
#' The formula used is (Tate, 1955, p. 1081):
#' \deqn{r_{pb} = \frac{\bar{x}_2 - \bar{x}_1}{\sigma_x} \times \sqrt{p \times q}}:
#' 
#' With:
#' \deqn{p = \frac{n_1}{n}, q = \frac{n_2}{n}}
#' \deqn{\bar{x}_1 = \frac{\sum_{i=1}^{n_2} x_{i,1}}{n_1}}
#' \deqn{\bar{x}_2 = \frac{\sum_{i=1}^{n_2} x_{i,2}}{n_2}}
#' \deqn{\sigma = \sqrt{\frac{SS}{n}}}
#' \deqn{SS = \sum_{j=1}^{2} \sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}\right)^2}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_1}, the sample size of the first category
#' \item \eqn{n_2}, the sample size of the second category
#' \item \eqn{n}, the total sample size, i.e. \eqn{n = n_1 + n_2}
#' \item \eqn{x_{i,j}} is the \eqn{i}-th score in category \eqn{j}
#' }
#' 
#' The oldest formula I could find is from Soper (1914, p. 384), which somewhat re-written is:
#' \deqn{r_{pb} = \frac{\bar{x}_2 - \bar{x}}{\sigma_x} \times \frac{\sqrt{p \times q}}{q}}
#' 
#' Tate also gave another formula (Tate, 1954, p. 606):
#' \deqn{r_{pb} = \frac{\bar{x}_2 - \bar{x}_1}{\sqrt{SS}} \times \frac{n_1 \times n_2}{n}}
#' 
#' Friedman (1968, p. 245) uses the degrees of freedom and test-statistic from the Student t-test for independent samples:
#' \deqn{r_{pb} = \sqrt{\frac{t^2}{t^2 + df}}}
#' 
#' As mentioned in the introduction, it can also be calculated by converting the categories to binary values, and then determine the Pearson product-moment correlation coefficient between these binary values and the scores.
#' 
#' Note that all of these should give the same result.
#'  
#' @references 
#' Friedman, H. (1968). Magnitude of experimental effect and a table for its rapid estimation. *Psychological Bulletin, 70*(4), 245-251. https://doi.org/10.1037/h0026258
#' 
#' Soper, H. E. (1914). On the probable error of the bi-serial expression for the correlation coefficient. *Biometrika, 10*(2/3), 384-390. https://doi.org/10.2307/2331789
#' 
#' Tate, R. F. (1954). Correlation between a discrete and a continuous variable. Point-biserial correlation. *The Annals of Mathematical Statistics, 25*(3), 603-607. https://doi.org/10.1214/aoms/1177728730
#' 
#' Tate, R. F. (1955). Applications of correlation models for biserial data. *Journal of the American Statistical Association, 50*(272), 1078-1095. https://doi.org/10.2307/2281207
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' file = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' dfr = read.csv(file, sep=';', na.strings=c("", "NA"))
#' r_point_biserial(dfr[['Gen_Gender']], dfr[['Over_Grade']])
#' 
#' @export
r_point_biserial <- function(catField, scaleField, categories=NULL){
  
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
  combined = c(X1, X2)
  
  # sample sizes
  n1 = length(X1)
  n2 = length(X2)
  n = n1 + n2
  
  # sample proportions
  p = n1/n
  q = n2/n
  
  # means and overall population standard deviation
  m1 = mean(X1)
  m2 = mean(X2)
  s = sqrt((n-1)/n) * sd(combined)
  
  # point-biserial correlation
  r_pb = (m2 - m1)/s * (p*q)**0.5
  
  results <- data.frame(cat1, cat2, m1, m2, r_pb)
  colnames(results) = c("cat. 1", "cat. 2", "mean 1", "mean 2", "r_pb")
  
  return(results)
  
}



