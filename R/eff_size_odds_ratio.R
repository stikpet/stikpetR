#' Odds Ratio
#' 
#' @description
#' Determines the odds ratio from a 2x2 table.
#' 
#' Odds can sometimes be reported as 'a one in five odds', but sometimes as 1 : 4. This later notation is less often seen, but means for every one event on the left side, there will be four on the right side.
#' 
#' The Odds is the ratio of that something will happen, over the probability that it will not. For the Odds Ratio, we compare the odds of the first category with the second group.
#' 
#' If the result is 1, it indicates that one variable has no influence on the other. A result higher than 1, indicates the odds are higher for the first category. A result lower than 1, indicates the odds are lower for the first.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' 
#' @returns 
#' Dataframe with:
#' \item{OR}{the odds ratio}
#' \item{n}{the sample size}
#' \item{statistic}{the test statistic (z-value)}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' 
#' The formula used is (Fisher, 1935, p. 50):
#' \deqn{OR = \frac{a/c}{b/d} = \frac{a\times d}{b\times c}}
#'  
#' *Symbols used:*
#' \itemize{
#' \item \eqn{a} the count in the top-left cell
#' \item \eqn{b} the count in the top-right cell
#' \item \eqn{c} the count in the bottom-left cell
#' \item \eqn{d} the count in the bottom-right cell
#' \item \eqn{\Phi\left(\dots\right)} the cumulative density function of the standard normal distribution
#' }
#' 
#' As for the test (McHugh, 2009, p. 123):
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{1}{a} + \frac{1}{b} + \frac{1}{c} + \frac{1}{d}}}
#' \deqn{z = \frac{\ln{\left(OR\right)}}{SE}}
#' 
#' The p-value is for the null-hypothesis that the population OR is 1.
#' 
#' The term Odds Ratio can for example be found in Cox (1958, p. 222).
#' 
#' @section Alternatives:
#' 
#' R's *stats* library has a function that also shows an odds ratio: *fisher.test()*
#' 
#' @seealso 
#' \code{\link{th_odds_ratio}}, rules of thumb for odds ratio
#' 
#' \code{\link{es_convert}}, to convert an odds ratio to Yule Q, Yule Y, or Cohen d.
#' 
#' @references
#' Cox, D. R. (1958). The regression analysis of binary sequences. *Journal of the Royal Statistical Society: Series B (Methodological), 20*(2), 215–232. https://doi.org/10.1111/j.2517-6161.1958.tb00292.x
#' 
#' Fisher, R. A. (1935). The logic of inductive inference. *Journal of the Royal Statistical Society, 98*(1), 39–82. https://doi.org/10.2307/2342435
#' 
#' McHugh, M. (2009). The odds ratio: Calculation, usage, and interpretation. *Biochemia Medica, 19*(2), 120–126. https://doi.org/10.11613/BM.2009.011
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_odds_ratio(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
es_odds_ratio <- function(field1, field2, categories1=NULL, categories2=NULL){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  #The odds ratio:
  or = (a/c) / (b/d)
  
  #Significance
  L = log(or)
  SE = sqrt(sum(1/ct))
  z = L/SE
  n = a + b + c +d
  pValue = 2*(1 - pnorm(abs(z)))
  
  # the results
  results = data.frame(or, n, z, pValue)
  colnames(results)<-c("OR", "n", "statistic", "p-value")
  
  return (results)
}