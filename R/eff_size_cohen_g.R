#' Cohen's g
#' 
#' @param data A vector with the data
#' @param codes Optional vector with the two codes to use
#' @return Cohen's g as a single value
#' 
#' @description 
#' Cohen’s g (Cohen, 1988) is an effect size measure that could be accompanying a one-sample binomial (see Rosnow & Rosenthal, 2003), 
#' score or Wald test. It is simply the difference of the sample proportion with 0.5. 
#' 
#' A video explanation of Cohen g can be found \href{https://youtu.be/tPZMvB8QrM0}(here on YouTube)
#' 
#' @details 
#' The formula used is (Cohen, 1988, p. 147):
#' \deqn{g=p-0.5}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{p} is the sample proportion
#' }
#'  
#' 
#' ## Alternative
#' 
#' I'm not aware of any alternative library that has this function.
#' 
#' @examples
#' data <- c("Female", "Male", "Male", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Female", "Male", "Male")
#' es_cohen_g(data)
#' 
#' @seealso 
#' This effect size could be used with a one-sample binomial test (\code{\link{ts_binomial_os}}), 
#' score test (\code{\link{ts_score_os}}), or Wald test (\code{\link{ts_wald_os}}).
#' 
#' Alternatives for Cohen g are the Alternative Ratio (\code{\link{es_alt_ratio}}) or Cohen's h' (\code{\link{es_cohen_h_os}})
#' 
#' For a rule of thumb interpretation use \code{\link{th_cohen_g}}
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#'  
#' @author 
#' P. Stikker
#' 
#' Please visit: \href{https://PeterStatistics.com}(PeterStatistics.com)
#' 
#' YouTube channel: \href{https://www.youtube.com/stikpet}(stikpet)
#' 
#' @export
es_cohen_g <- function(data, codes=NULL){
  
  if (is.null(codes)){
    freq <- table(data)
    prop <- freq / sum(freq)
    p1 <- unname(prop[1])
    
  } else {
    #Determine number of successes
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
    
    #Determine total sample size
    n<-n1 + n2
    
    p1 <- n1/n
    
  }
  
  g <- p1 - 0.5
  
  
  return (g)
  
}