#' Trinomial Test (Paired Samples)
#' @description 
#' A similar test as the sign test, but also includes the pairs that are tied.
#' 
#' @param field1 the numeric scores of the first variable
#' @param field2 the numeric scores of the second variable
#' @param levels vector, optional. the levels from field1 and field2
#' @param dmu float, optional. The difference according to the null hypothesis (default is 0)
#' 
#' @returns 
#' A dataframe with:
#' \item{n pos}{the number of scores with a positive difference}
#' \item{n neg}{the number of scores with a negative difference}
#' \item{n 0}{the number of scores with a no difference}
#' \item{p-Value}{the significance (p-value)}
#' 
#' @details
#' The formula used  (Bian et al., 2009, p. 6):
#' \deqn{sig. = 2\times \text{TRI}\left(\left(n_{pos}, n_{neg}, n_0\right), \left(p_{pos}, p_{neg}, p_0\right)\right)}
#' 
#' With:
#' \deqn{n_{pos} = \sum_{i=1}^n \begin{cases} 1 & \text{ if } d_i>d_{H0} \\ 0 & \text{ if } d_i \le d_{H0} \end{cases}}
#' \deqn{n_{neg} = \sum_{i=1}^n \begin{cases} 0 & \text{ if } d_i \ge d_{H0} \\ 1 & \text{ if } d_i < d_{H0} \end{cases}}
#' \deqn{n_{0} = \sum_{i=1}^n \begin{cases} 1 & \text{ if } d_i = d_{H0} \\ 0 & \text{ if } d_i \ne d_{H0} \end{cases}}
#' \deqn{d_i = x_i - y_i}
#' \deqn{p_0 = \frac{n_0}{n}}
#' \deqn{p_{pos}=p_{neg}=\frac{1 - p_0}{2}}
#' 
#' The cumulative mass function of the trinomial distribution is then calculated using:
#' \deqn{\text{TRI}\left(\left(n_{pos}, n_{neg}, n_0\right), \left(p_{pos}, p_{neg}, p_0\right)\right) = \sum_{i=n_d}^n\sum_{j=0}^{\left\lfloor\frac{n-i}{2}\right\rfloor} \text{tri}\left(\left(j, j+i, n-j-\left(j+i\right)\right), \left(p_{pos}, p_{neg}, p_0\right)\right)}
#' \deqn{n_d = \left|n_{pos} - n_{neg}\right|}
#' 
#' The probability mass function of the trinomial distribution is (Bian et al., 2009, p. 5):
#' \deqn{\text{tri} = \left(\left(n_a, n_b, n_c\right), \left(p_{a}, p_b, p_c\right)\right) = \frac{n!}{a!\times b! \times c!}\times p_a^{n_a}\times p_b^{n_b}\times p_c^{n_c}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n} is the number of pairs with a difference unequal to zero
#' \item \eqn{n_{pos}} the number of pairs with a difference greater than the null hypothesis
#' \item \eqn{n_{neg}} the number of pairs with a difference greater than the null hypothesis
#' \item \eqn{n_{0}} the number of pairs with no difference with the null hypothesis
#' \item \eqn{d_{H0}} the difference according to the null hypothesis, usually 0
#' \item \eqn{x_i} the i-th score from the first variable
#' \item \eqn{y_i} the i-th score from the second variable
#' \item \eqn{\text{tri}\left(\dots,\dots\right)} the probability mass function of the trinomial distribution
#' }
#' 
#' **Alternatives**
#' 
#' *library(EMT)*
#' 
#' datFrame = na.omit(data.frame(ord1,ord2))
#' 
#' d = datFrame$ord1 - datFrame$ord2
#' 
#' pos = sum(d>0)
#' 
#' neg = sum(d<0)
#' 
#' ties = sum(d==0)
#' 
#' n = pos + neg + ties
#' 
#' p0 = ties/n
#' 
#' p1 = (1 - p0)/2
#' 
#' multinomial.test(c(pos, neg, ties), c(p1, p1, p0))
#' 
#' @references 
#' Bian, G., McAleer, M., & Wong, W.-K. (2009). A trinomial test for paired data when there are many ties. *SSRN Electronic Journal*. https://doi.org/10.2139/ssrn.1410589
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_trinomial_ps <- function(field1, field2, levels=NULL, dmu=0){
  
  if (!is.null(levels)){
    myFieldOrd = factor(field1, ordered = TRUE, levels = levels)
    field1 = as.numeric(myFieldOrd)
  }
  
  if (!is.null(levels)){
    myFieldOrd = factor(field2, ordered = TRUE, levels = levels)
    field2 = as.numeric(myFieldOrd)
  }
  
  datFrame = na.omit(data.frame(field1,field2))
  
  d = datFrame$field1 - datFrame$field2
  
  pos = sum(d>dmu)
  neg = sum(d<dmu)
  ties = sum(d==dmu)
  
  n = pos + neg + ties
  
  p0 = ties/n
  p1 = (1 - p0)/2
  k = abs(pos - neg)
  
  sig=0
  for (z in k:n) {
    for (i in 0:floor((n - z)/2)) {
      sig = sig + dmultinom(c(i, i+z, n - i - (i+z)), prob = c(p1, p1, p0))
    }
  }
  
  pValue = sig*2  
  
  results = data.frame(pos, neg, ties, pValue)
  colnames(results) = c("n pos", "n neg", "n 0", "p-value")
  return(results)
  
}