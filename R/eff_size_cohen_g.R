#' Cohen's g
#' 
#' @param data A vector with the data
#' @param codes Optional vector with the two codes to use
#' @return Cohen's g as a single value
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
#' For a classification use: *th_cohen_g(g)*
#' 
#' **Alternative**
#' 
#' I'm not aware of any alternative library that has this function.
#' 
#' @examples
#' data <- c("Female", "Male", "Male", "Female", "Male", "Male", "Female", "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Female", "Male", "Female", "Male", "Male")
#' es_cohen_g(data)
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
