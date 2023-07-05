#' Cohen's h' 
#' 
#' @description 
#' An adaptation of Cohen h (\code{\link{es_cohen_h}}) for a one-sample case. It is an effect size measure that could 
#' be accompanying a one-sample binomial, score or Wald test.
#' 
#' A [YouTube](https://youtu.be/ddWe94VKX_8) video on Cohen h'.
#' 
#' @param data a vector with the data
#' @param codes optional vector with the two codes to use
#' @param p0 optional the hypothesized proportion for the first category (default is 0.5)
#' 
#' @return Cohen's h'
#' 
#' @details 
#' Formula used (Cohen, 1988, p. 202):
#' \deqn{h'=\phi_{1}-\phi_{h_0}}
#' With:
#' \deqn{\phi_{i}=2\times\textup{arcsin}\sqrt{p_{i}}}
#' \deqn{p_i = \frac{F_i}{n}}
#' \deqn{n = \sum_{i=1}^k F_i}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{F_i} is the (absolute) frequency (count) of category i
#' \item \eqn{n} is the sample size, i.e. the sum of all frequencies
#' \item \eqn{p_i} the proportion of cases in category i
#' \item \eqn{p_{h_0}} the expected proportion (i.e. the proportion according to the null hypothesis)
#' }
#' 
#' @section Alternatives:
#' 
#' The *'pwr'* library has a similar function: *ES.h()*
#' 
#' @seealso 
#' \code{\link{es_convert}}, to convert Cohen h' to Cohen h, use fr="cohenhos" and to=cohenh
#' 
#' \code{\link{th_cohen_h}}, rules-of-thumb for Cohen g
#'   
#' @references 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: Numeric list
#' ex1 = c(1, 1, 2, 1, 2, 1, 2, 1)
#' es_cohen_h_os(ex1)
#' es_cohen_h_os(ex1, p0=0.3)
#' 
#' #Example 2: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_cohen_h_os(df1['sex'])
#' 
#' @export
es_cohen_h_os <- function(data, codes=NULL, p0=0.5){
  
  data = na.omit(data)
  
  if (is.null(codes)){
    freq <- table(data)
    n1 <- sum(data==rownames(freq)[1])
    n <- sum(freq)
    n2 = n - n1    
  } else {
    #Determine number of successes
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
    n<-n1 + n2
  }
  
  p1 <- n1/n
  
  phi1 = 2 * asin(sqrt(p1))
  phic = 2 * asin(sqrt(p0))
  
  h2 = phi1 - phic
  
  return (h2)

}
