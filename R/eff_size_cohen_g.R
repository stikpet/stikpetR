#' Cohen's g
#' 
#' @description 
#' Cohen’s g (Cohen, 1988) is an effect size measure that could be accompanying a one-sample binomial (see Rosnow & Rosenthal, 2003), 
#' score or Wald test. It is simply the difference of the sample proportion with 0.5. 
#' 
#' A video explanation of Cohen g can be found \href{https://youtu.be/tPZMvB8QrM0}{here on YouTube}
#' 
#' @param data vector with the data
#' @param codes optional vector with the two codes to use
#' 
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
#' @seealso 
#' \code{\link{th_cohen_g}}, rules-of-thumb for Cohen g
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @examples
#' data <- c("Female", "Male", "Male", "Female", "Male", "Male")
#' es_cohen_g(data)
#' 
#' @export
es_cohen_g <- function(data, codes=NULL){
  
  data = na.omit(data)
  
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