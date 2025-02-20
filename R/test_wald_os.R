#' One-Sample Wald Test
#' 
#' @description 
#' A one-sample score test could be used with binary data, to test if the two categories
#' have a significantly different proportion. It is an approximation of a binomial test, 
#' by using a standard normal distribution. Since the binomial distribution is discrete 
#' while the normal is continuous, a so-called continuity correction can (should?) be 
#' applied.
#' 
#' The null hypothesis is usually that the proportions of the two categories in the 
#' population are equal (i.e. 0.5 for each). If the p-value of the test is below the 
#' pre-defined alpha level (usually 5% = 0.05) the null hypothesis is rejected and 
#' the two categories differ in proportion significantly.
#' 
#' The input for the function doesn't have to be a binary variable. 
#' A nominal variable can also be used and the two categories to compare indicated.
#' 
#' A significance in general is the probability of a result as in the sample, 
#' or more extreme, if the null hypothesis is true. 
#' 
#' Some info on the different tests can be found in \href{https://youtu.be/jQ-nSPTGOgE}{video}.
#' This function is shown in this [YouTube video](https://youtu.be/1_KPlWx2vLk) and the test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tests/SignOneSample.html)
#' 
#' @param data A vector with the data
#' @param p0 Optional hypothesized proportion for the first category (default is 0.5)
#' @param p0Cat Optional the category for which p0 was used
#' @param codes Optional vector with the two codes to use
#' @param cc use of continuity correction (default is "none")
#' 
#' @returns 
#' Dataframe with:
#' \item{n}{the sample size}
#' \item{statistic}{the test value}
#' \item{pValue}{two-sided p-value}
#' \item{test}{a description of the test used}
#' 
#' @details 
#' To decide on which category is associated with p0 the following is used:
#' \itemize{
#' \item If codes are provided, the first code is assumed to be the category for the p0.
#' \item If p0Cat is specified that will be used for p0 and all other categories will be considered as category 2, this means if there are more than two categories the remaining two or more (besides p0Cat) will be merged as one large category.
#' \item If neither codes or p0Cat is specified and more than two categories are in the data a warning is printed and no results.
#' \item If neither codes or p0Cat is specified and there are two categories, p0 is assumed to be for the category closest matching the p0 value (i.e. if p0 is above 0.5 the category with the highest count is assumed to be used for p0)
#' }
#' 
#' This test differs from the one-sample score test in the calculation of the standard error. 
#' For the ‘regular’ version this is based on the expected proportion, 
#' while for the Wald version it is done with the observed proportion.
#' 
#' The formula used (Wald, 1943):
#' \deqn{z=\frac{x - \mu}{SE}}
#' With:
#' \deqn{\mu = n\times p_0}
#' \deqn{SE = \sqrt{x\times\left(1 - \frac{x}{n}\right)}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x} is the number of successes in the sample
#' \item \eqn{p_0} the expected proportion (i.e. the proportion according to the null hypothesis)
#' }
#' 
#' If the Yates continuity correction is used the formula changes to (Yates, 1934, p. 222):
#' \deqn{z_{Yates} = \frac{\left|x - \mu\right| - 0.5}{SE}}
#' 
#' The formula used in the calculation is the one from IBM (2021, p. 997).
#' IBM refers to Agresti, most likely Agresti (2013, p. 10), who in turn
#' refer to Wald (1943)
#' 
#' @section Before, After and Alternatives:
#' Before running the test you might first want to get an impression using a frequency table: tab_frequency
#' 
#' After the test you might want an effect size measure:
#' \code{\link{es_cohen_g}}, for Cohen g
#' \code{\link{es_cohen_h_os}}, for Cohen h'
#' \code{\link{es_alt_ratio}}, for Alternative Ratio
#' 
#' Alternatives for this test could be:
#' \code{\link{ts_binomial_os}}, for One-Sample Binomial Test
#' \code{\link{ts_score_os}}, for One-Sample Score Test
#' 
#' @references 
#' Agresti, A. (2013). *Categorical data analysis* (3rd ed.). Wiley.
#' 
#' IBM SPSS Statistics Algorithms. (2021). IBM.
#' 
#' Wald, A. (1943). Tests of statistical hypotheses concerning several parameters when the number of observations is large. *Transactions of the American Mathematical Society, 54*(3), 426–482. doi:10.2307/1990256
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. doi:10.2307/2983604
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: Numeric list
#' ex1 = c(1, 1, 2, 1, 2, 1, 2, 1)
#' ts_wald_os(ex1)
#' ts_wald_os(ex1, p0=0.3)
#' ts_wald_os(ex1, p0=0.3, cc="yates")
#' 
#' #Example 2: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ts_wald_os(df1['sex'])
#' ts_wald_os(df1['mar1'], codes=c("DIVORCED", "NEVER MARRIED"))
#' 
#' @export
ts_wald_os <- function(data,
                       p0 = 0.5,
                       p0Cat = NULL,
                       codes=NULL,
                       cc=c("none", "yates")){
  
  if (length(cc)>1) {cc = "none"}
  
  data = na.omit(data)
  
  #if no codes provided use first found
  if (is.null(codes)) {
    freq = table(data)
    
    if (is.null(p0Cat)){
      #check if there were exactly two categories or not
      if (length(freq) != 2){
        # unable to determine which category p0 would belong to, so print warning and end
        print("WARNING: data does not have two unique categories, please specify two categories using codes parameter")
        return
      }
      else{
        n1 = unname(freq[1])
        n2 = unname(freq[2])
        n = n1 + n2
        
        #determine p0 was for which category
        p0_cat = names(freq)[1]
        if (p0 > 0.5 & n1 < n2){
          n3 = n2
          n2 = n1
          n1 = n3
          p0_cat = names(table(data))[2]
        }
        
        cat_used = paste0(" (assuming p0 for ", p0_cat, ")")
      }
    }
    else {
      n = sum(table(data))
      n1 = sum(data==p0Cat)
      n2 = n - n1
      p0_cat = p0Cat
      cat_used = paste0(" (with p0 for ", p0Cat, ")")
    }
  }
  
  else{
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
    n = n1 + n2
    cat_used = paste0(" (with p0 for ", codes[1], ")")
  }
  
  minCount = n1
  ExpProp = p0
  if (n2 < n1){
    minCount = n2
    ExpProp = 1 - ExpProp}
  
  #Wald approximation
  if (cc == "none"){
    p = minCount / n
    q = 1 - p
    se = (p * q / n) ^ 0.5
    Z = (p - ExpProp) / se
    pValue = 2 * (1 - pnorm(abs(Z)))
    statistic = Z
    testUsed = "one-sample Wald"}
  
  else if (cc == "yates"){
    #Wald approximation with continuity correction
    p = (minCount + 0.5) / n
    q = 1 - p
    se = (p * q / n) ^ 0.5
    Z = (p - ExpProp) / se
    pValue = 2 * (1 - pnorm(abs(Z)))
    statistic = Z
    testUsed = "one-sample Wald with Yates continuity correction"}
  
  testUsed = paste0(testUsed, cat_used)
  testResults <- data.frame(n, statistic, pValue, testUsed)
  colnames(testResults)<-c("n", "statistic", "p-value (2-sided)", "test")
  
  return (testResults)
}