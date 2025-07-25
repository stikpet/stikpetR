#' Cohen's h' 
#' 
#' @description 
#' An adaptation of Cohen h (\code{\link{es_cohen_h}}) for a one-sample case. It is an effect size measure that could 
#' be accompanying a one-sample binomial, score or Wald test.
#' 
#' A [YouTube](https://youtu.be/ddWe94VKX_8) video on Cohen h'.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/sGfFB7Zzeas) and the effect size is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/CohenH.html)
#' 
#' @param data a vector with the data
#' @param p0 Optional hypothesized proportion for the first category (default is 0.5)
#' @param p0Cat Optional the category for which p0 was used
#' @param codes Optional vector with the two codes to use
#' 
#' @returns 
#' Dataframe with:
#' \item{Cohen h'}{the Cohen h' value}
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
#' @section Before, After and Alternatives:
#' Before this effect size you might first want to perform a test:
#' \code{\link{ts_binomial_os}}, for One-Sample Binomial Test
#' \code{\link{ts_score_os}}, for One-Sample Score Test
#' \code{\link{ts_wald_os}}, for One-Sample Wald Test
#' 
#' After this, you might want a rule-of-thumb or first convert this to a 'regular' Cohen h:
#' \code{\link{es_convert}}, to convert Cohen h' to Cohen h, use fr="cohenhos" and to=cohenh
#' \code{\link{th_cohen_h}}, for rules-of-thumb for Cohen h
#' 
#' Alternatives for this effect size could be:
#' \code{\link{es_cohen_g}}, for Cohen g
#' \code{\link{es_alt_ratio}}, for Alternative Ratio
#' \code{\link{r_rosenthal}}, for Rosenthal Correlation if a z-value is available
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
es_cohen_h_os <- function(data, p0=0.5, p0Cat=NULL, codes=NULL){
  
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
  
  p1 <- n1/n
  
  phi1 = 2 * asin(sqrt(p1))
  phic = 2 * asin(sqrt(p0))
  
  h2 = phi1 - phic
  
  results <- data.frame(h2, cat_used)
  colnames(results)<-c("Cohen h'", "comment")
  
  return(results)

}



