#' Alternative Ratio
#' 
#' @description 
#' The Alternative Ratio is an effect size measure that could be accompanying a one-sample binomial, score or Wald test.It is simply the sample proportion (percentage), divided by the expected population proportion (often set at 0.5)
#' 
#' The Alternative Ratio is only mentioned in the documentation of a program called PASS from NCSS (n.d.), and referred to as Relative Risk by JonB (2015).
#' 
#' This function is shown in this [YouTube video](https://youtu.be/IwZph3C9xFc) and the effect size is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/AlternativeRatio.html)
#' 
#' @param data vector with the data
#' @param p0 Optional hypothesized proportion for the first category (default is 0.5)
#' @param p0Cat Optional the category for which p0 was used
#' @param codes Optional vector with the two codes to use
#' 
#' @returns 
#' Dataframe with:
#' \item{AR1}{the alternative category for one category}
#' \item{AR2}{the alternative category for the other category}
#' \item{comment}{the category for which p0 was}
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
#' The formula used is:
#' \deqn{AR=\frac{p}{\pi}}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{p} is the sample proportion of one of the categories
#' \item \eqn{\pi} the expected proportion
#' } 
#' 
#' @section Before, After and Alternatives:
#' Before this effect size you might first want to perform a test:
#' \code{\link{ts_binomial_os}}, for One-Sample Binomial Test
#' \code{\link{ts_score_os}}, for One-Sample Score Test
#' \code{\link{ts_wald_os}}, for One-Sample Wald Test
#' 
#' Unfortunately I'm not aware of any rule-of-thumb for this measure.
#' 
#' Alternatives for this effect size could be:
#' \code{\link{es_cohen_g}}, for Cohen g
#' \code{\link{es_cohen_h_os}}, for Cohen h'
#' \code{\link{r_rosenthal}}, for Rosenthal Correlation if a z-value is available
#' 
#' @references
#' JonB. (2015, October 14). Effect size of a binomial test and its relation to other measures of effect size. StackExchange - Cross Validated. https://stats.stackexchange.com/q/176856
#' 
#' NCSS. (n.d.). Tests for one proportion. In PASS Sample Size Software (pp. 100-1-100-132). Retrieved November 10, 2018, from https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Tests_for_One_Proportion.pdf
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @examples
#' # Example 1: Numeric list
#' ex1 = c(1, 1, 2, 1, 2, 1, 2, 1)
#' es_alt_ratio(ex1)
#' es_alt_ratio(ex1, p0=0.3)
#' 
#' # Example 2: Text list
#' ex2 = c("Female", "Male", "Male", "Female", "Male", "Male")
#' es_alt_ratio(ex2)
#' es_alt_ratio(ex2, p0Cat='Female')
#' es_alt_ratio(ex2, codes=c('Male', 'Female'))
#' 
#' # Example 3: dataframe
#' file1 <- "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(file1, sep=",", na.strings=c("", "NA"))
#' es_alt_ratio(df1['sex'])
#' es_alt_ratio(df1['mar1'], codes=c("DIVORCED", "NEVER MARRIED"))
#' 
#' @export
es_alt_ratio <- function(data, p0=0.5, p0Cat=NULL, codes=NULL){
  
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
        
        cat_used = paste0("assuming p0 for ", p0_cat)
      }
    }
    else {
      n = sum(table(data))
      n1 = sum(data==p0Cat)
      n2 = n - n1
      p0_cat = p0Cat
      cat_used = paste0("with p0 for ", p0Cat)
    }
  }
  
  else{
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
    n = n1 + n2
    cat_used = paste0("with p0 for ", codes[1])
  }
  
  p1 = n1 / n
  p2 = n2 / n
  
  AR1 = p1 / p0
  AR2 = p2 / (1 - p0)
  
  results <- data.frame(AR1, AR2, cat_used)
  colnames(results)<-c("Alt. Ratio Cat. 1", "Alt. Ratio Cat. 2", "comment")
  
  return(results)
}



