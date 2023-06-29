#' Alternative Ratio
#' 
#' @description 
#' The Alternative Ratio is an effect size measure that could be accompanying a one-sample binomial, score or Wald test.It is simply the sample proportion (percentage), divided by the expected population proportion (often set at 0.5)
#' 
#' The Alternative Ratio is only mentioned in the documentation of a program called PASS from NCSS (n.d.), and referred to as Relative Risk by JonB (2015).
#' 
#' @param data vector with the data
#' @param codes optional vector with the two codes to use
#' @param p0 optional the hypothesized proportion for the first category (default is 0.5)
#' @param category optional category to label as 'success', otherwise the first category found is used.
#' 
#' @returns 
#' Dataframe with:
#' \item{AR1}{the alternative category for one category}
#' \item{AR2}{the alternative category for the other category}
#' 
#' @details 
#' 
#' If codes and category are not provided the first category will be the first data point.
#' 
#' If codes only are provided the first category in the codes is used.
#' 
#' The formula used is:
#' \deqn{AR=\frac{p}{\pi}}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{p} is the sample proportion of one of the categories
#' \item \eqn{\pi} the expected proportion
#' } 
#' 
#' @section Alternatives:
#' 
#' I'm not aware of any alternative library that has this function.
#' 
#' @examples
#' data <- c("Female", "Male", "Male", "Female", "Male", "Male")
#' es_alt_ratio(data)
#' es_alt_ratio(data, category="Male")
#' es_alt_ratio(data, c("Male", "Female"))
#' 
#' @references
#' JonB. (2015, October 14). Effect size of a binomial test and its relation to other measures of effect size. StackExchange - Cross Validated. https://stats.stackexchange.com/q/176856
#' 
#' NCSS. (n.d.). Tests for one proportion. In PASS Sample Size Software (pp. 100-1-100–132). Retrieved November 10, 2018, from https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Tests_for_One_Proportion.pdf
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
es_alt_ratio <- function(data, codes=NULL, p0=0.5, category=NULL){
  
  data = na.omit(data)
  
  if (is.null(codes)){
    freq <- table(data)
    n <- sum(freq)
    
    if (is.null(category)){
      n1 <- sum(data==rownames(freq)[1])}
    else{
      n1<-sum(data==category)}
    
    n2 = n - n1
  }
  else{
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
    n<-n1 + n2
    
    if(!is.null(category)){
      if(codes[2]==category){
        n3 = n1
        n1 = n2
        n2 = n3}
    }
  }
  
  p1 = n1 / n
  p2 = n2 / n
  
  AR1 = p1 / p0
  AR2 = p2 / (1 - p0)
  
  results <- data.frame(AR1, AR2)
  colnames(results)<-c("Alt. Ratio Cat. 1", "Alt. Ratio Cat. 2")
  
  return (results)
  
}