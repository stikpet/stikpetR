#' Cohen's g
#' 
#' @description 
#' Cohen’s g (Cohen, 1988) is an effect size measure that could be accompanying a one-sample binomial (see Rosnow & Rosenthal, 2003), 
#' score or Wald test. It is simply the difference of the sample proportion with 0.5. 
#' 
#' A video explanation of Cohen g can be found \href{https://youtu.be/tPZMvB8QrM0}{here on YouTube}
#' 
#' @param data vector with the data
#' @param p0Cat Optional the category for which p0=0.5 was used
#' @param codes Optional vector with the two codes to use
#' 
#' @returns 
#' Dataframe with:
#' \item{g for cat 1}{Cohen g for category 1}
#' \item{g for cat 2}{Cohen g for category 2}
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
#' The formula used is (Cohen, 1988, p. 147):
#' \deqn{g=p-0.5}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{p} is the sample proportion
#' }
#'  
#' @section Before, After and Alternatives:
#' Before this effect size you might first want to perform a test:
#' \code{\link{ts_binomial_os}}, for One-Sample Binomial Test
#' \code{\link{ts_score_os}}, for One-Sample Score Test
#' \code{\link{ts_wald_os}}, for One-Sample Wald Test
#' 
#' After this, you might want a rule-of-thumb:
#' \code{\link{th_cohen_g}}, for rules-of-thumb for Cohen g
#' 
#' Alternatives for this effect size could be:
#' \code{\link{es_cohen_h_os}}, for Cohen h'
#' \code{\link{es_alt_ratio}}, for Alternative Ratio
#' \code{\link{r_rosenthal}}, for Rosenthal Correlation if a z-value is available
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example 1: Numeric list
#' ex1 = c(1, 1, 2, 1, 2, 1, 2, 1)
#' es_cohen_g(ex1)
#' 
#' #Example 2: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' es_cohen_g(df1['sex'])
#' #Using two specific categories:
#' es_cohen_g(df1['mar1'], codes=c("DIVORCED", "NEVER MARRIED"))
#' 
#' @export
es_cohen_g <- function(data, p0Cat=NULL, codes=NULL){
  
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
        
        #assume p0 was for first category
        cat1_lbl = names(freq)[1]
        cat2_lbl = names(freq)[2]
        
      }
    }
    else {
      n = sum(table(data))
      n1 = sum(data==p0Cat)
      n2 = n - n1
      cat1_lbl = p0Cat
      if (length(freq)==2){
        if (cat1_lbl == names(freq)[1]){
          cat2_lbl = names(freq)[2]
        }
        else{
          cat2_lbl = names(freq)[1]
        }
      }
      else{
        cat2_lbl = "all other"
      }
    }
  }
  
  else{
    n1<-sum(data==codes[1])
    n2<-sum(data==codes[2])
    n = n1 + n2
    cat1_lbl = codes[1]
    cat2_lbl = codes[2]
  }
  
  p1 = n1/n
  p2 = 1 - p1
  g1 <- p1 - 0.5
  g2 <- p2 - 0.5
  
  results <- data.frame(g1, g2)
  colnames(results)<-c(paste0("g for ", cat1_lbl), paste0("g for ", cat2_lbl))
  
  return (results)
}