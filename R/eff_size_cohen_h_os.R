#' Cohen's h' 
#' 
#' @param data A vector with the data
#' @param codes Optional vector with the two codes to use
#' @param p0 Optional the hypothesized proportion for the first category (default is 0.5)
#' @return Cohen's h'
#' 
#' @details 
#' An adaptation of Cohen h for a one-sample case.
#' 
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
#' For classification rule-of-thumb use a conversion to Cohen h: *es_convert(h2, from="cohenhos", to="cohenh")*
#' 
#' Then use: *th_cohen_h()*
#' 
#' **Alternative**
#' 
#' The *'pwr'* library has a similar function: *ES.h()*
#' 
#' @examples 
#' data <- c("Female", "Male", "Male", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Female", "Male", "Male")
#' es_cohen_h_os(data)
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' @export
es_cohen_h_os <- function(data, codes=NULL, p0=0.5){
  
  if (is.null(codes)){
    freq <- table(data)
    n1 <- sum(data==data[1])
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
