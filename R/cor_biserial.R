#' Biserial Correlation Coefficient
#' @description
#' This is an extension of the point-biserial correlation coefficient, if the categories come from a so-called latent normally distributed scale. This is the case if scores were categorized and then compared to some other numeric scores (e.g. grades being categorized into pass/fail, and then use this pass/fail to correlate with age).
#' 
#' As the name implies a correlation coefficient indicates how two variables co-relate, i.e. if one goes up is it likely for the other to go up or down. A zero would indicate there is not (linear) relation, while a -1 would mean a perfect negative correlation (if one goes up, the other goes down, and vice versa), and a +1 a perfect positive correlation (if one goes up, the other also goes up, and vice versa).
#' 
#' With two categories we could read this more as if the score go up and there is a positive correlation, it is more likely that it came from a category 1 case, rather than a category 0.
#' 
#' There is a warning though that if one of the two categories has very small sample size compared to the other, this coefficient will not be very accurate (Soper, 1914, p.390; Jacobs & Viechtbauer, 2017, p. 165). Soper (1914, p. 390) warns to use this if one category is 4% or less from the combined sample size. On a website ChangingMinds someone posted as limit 10% (ChangingMinds, n.d.), unfortunately without a source.
#' 
#' The coefficient is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Correlations/Biserial.html)
#' 
#' @param catField A vector with the categorical data
#' @param scaleField A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' 
#' @returns 
#' A dataframe with:
#' \item{cat. 0}{the category that was used as category 0}
#' \item{cat. 1}{the category that was used as category 1}
#' \item{n1/n}{the proportion of scores in the category 1}
#' \item{mean 0}{the arithmetic mean of the scores from category 0}
#' \item{mean 1}{the arithmetic mean of the scores from category 1}
#' \item{r_b}{the biserial correlation coefficient}
#' 
#' @details
#' The formula used is (Tate, 1955a, p. 1087):
#' \deqn{r_b = \frac{p \times q \times \left(\bar{x}_2 - \bar{x}_1\right)}{\sigma_x \times p_{z_p}}}
#' 
#' With:
#' \deqn{p = \frac{n_1}{n}, q = \frac{n_0}{n}}
#' \deqn{\bar{x}_0 = \frac{\sum_{i=1}^{n_0} x_{i,0}}{n_0}}
#' \deqn{\bar{x}_1 = \frac{\sum_{i=1}^{n_1} x_{i,1}}{n_1}}
#' \deqn{\sigma = \sqrt{\frac{SS}{n}}}
#' \deqn{SS = \sum_{j=1}^{2} \sum_{i=1}^{n_j} \left(x_{i,j} - \bar{x}\right)^2}
#' \deqn{z_p = \Phi^{-1}\left(p\right)}
#' \deqn{p_{z_p} = \phi\left(z_p\right)}
#' 
#' Symbols used:
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_0}, the sample size of the first category
#' \item \eqn{n_1}, the sample size of the second category
#' \item \eqn{n}, the total sample size, i.e. \eqn{n = n_1 + n_2}
#' \item \eqn{x_{i,j}} is the \eqn{i}-th score in category \eqn{j}
#' }
#' 
#' The oldest formula I could find is from Pearson (1909, p. 97), which somewhat re-written is:
#' \deqn{r_b = \frac{\frac{\bar{x}_1 - \bar{x}}{\sigma_x}}{\frac{p_{z_p}}{p}}}
#' 
#' Since divide by a fraction is multiplying by its inverse, Soper (1914, p. 384) has:
#' \deqn{r_b = \frac{\bar{x}_1 - \bar{x}}{\sigma_x} \times \frac{p}{p_{z_p}}}
#' 
#' If we were to create binary values of the categories, then Tate (1955a, p. 1079; 1955b, p. 207) used the covariance between these and the scores:
#' \deqn{r_b = \frac{\sigma_{bx}}{\sigma_x \times p_{z_p}}}
#' 
#' Not too surprising, since it can be shown that \(\sigma_{bx} = p \times q \times \left(\bar{x}_1 - \bar{x}_0\right)\)
#' 
#' Tata (1955a, p. 1087; 1955b, p. 207) also show a conversion using the point-rank biserial:
#' \deqn{r_b = r_{pb} \times \frac{\sigma_b}{p_{z_p}}}
#' 
#' Note that all of these should give the same result.
#' 
#' @references 
#' ChangingMinds. (n.d.). Biserial Correlation Coefficient. Retrieved July 18, 2025, from https://changingminds.org/explanations/research/analysis/biserial.htm
#' 
#' Jacobs, P., & Viechtbauer, W. (2017). Estimation of the biserial correlation and its sampling variance for use in meta‐analysis. *Research Synthesis Methods, 8*(2), 161–180. https://doi.org/10.1002/jrsm.1218
#' 
#' Pearson, K. (1909). On a new method of determining correlation between a measured character A, and a character B. *Biometrika, 7*(1/2), 96–105. https://doi.org/10.2307/2345365
#' 
#' Soper, H. E. (1914). On the probable error of the bi-serial expression for the correlation coefficient. *Biometrika, 10*(2/3), 384–390. https://doi.org/10.2307/2331789
#' 
#' Tate, R. F. (1955a). Applications of correlation models for biserial data. *Journal of the American Statistical Association, 50*(272), 1078–1095. https://doi.org/10.2307/2281207
#' 
#' Tate, R. F. (1955b). The theory of correlation between two continuous variables when one is dichotomized. *Biometrika, 42*(1/2), 205–216. https://doi.org/10.2307/2333437
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
r_biserial <- function(catField, scaleField, categories=NULL){
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
  
  # the normal distribution part
  z_p = qnorm(p)
  p_zp = dnorm(z_p) 
  
  # biserial correlation
  r_b = p*q*(m2 - m1)/(s * p_zp)
  
  #the results
  results <- data.frame(cat1, cat2, q, m1, m2, r_b)
  colnames(results) = c("cat. 0", "cat. 1", "n1/n", "mean 0", "mean 1", "r_b")
  
  return(results)
  }