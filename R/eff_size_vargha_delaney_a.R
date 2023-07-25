#' Vargha and Delaney A
#' @description
#' A generalisation of the Common Language Effect size for ordinal data. Taking a random case from each of the two categories, Vargha-Delaney A is then the probability that the score of the first category is higher than the second.
#' 
#' @param catField A vector with the scores data
#' @param ordField A vector with the group data
#' @param categories optional vector with categories to use and order for the categorical field. Otherwise the first two found will be used.
#' @param levels optional vector with the labels of the ordinal field in order.
#' 
#' @returns
#' A dataframe with:
#' \item{A1}{Vargha-Delaney A for the first category}
#' \item{A2}{Vargha-Delaney A for the second category}
#' 
#' @details
#' 
#' The formula used is (Vargha & Delaney, 2000, p. 107):
#' \deqn{A = \frac{1}{n_j}\times\left(\frac{R_i}{n_i} - \frac{n_i + 1}{2}\right)}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_i} the number of scores in category i
#' \item \eqn{R_i} the sum of the ranks in category i
#' }
#' 
#' This effect size is an adaptation of the Common Language effect size, 
#' adapted for ordinal data.
#' 
#' It could also be calculated from the Mann-Whitney U value:
#' \deqn{A = \frac{U}{n_1\times n_2}}
#' 
#' Note that the difference between the two options (using category 1 or category 2) will be the deviation from 0.5. If all scores in the first category are lower than the scores in the second, A will be 0 using the first category, and 1 for the second.
#' 
#' If the number of scores in the first category higher than the second, is the same as the other way around, A (no matter which category used) will be 0.5.
#' 
#' The VDA can be converted to a Rank Biserial (= Cliff delta) using the **es_convert()** function. This can then be converted to a Cohen d, and then the rules-of-thumb for Cohen d could be used (**th_cohen_d()**)
#' 
#' @seealso 
#' \code{\link{es_convert}}, to convert an VDA to a rank biserial use `fr="vda", to="rb"`. To convert the result to Cohen d, use `fr="rb", to="cohend"`.
#' 
#' \code{\link{th_cohen_d}}, rules of thumb for Cohen d
#' 
#' 
#' @references 
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. *Journal of Educational and Behavioral Statistics, 25*(2), 101–132. https://doi.org/10.3102/10769986025002101
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' myLevels = c('Not scientific at all', 'Not too scientific', 'Pretty scientific', 'Very scientific')
#' es_vargha_delaney_a(df1['sex'], df1['accntsci'], levels=myLevels)
#' 
#' #Example 2: vectors
#' binary = c("apple", "apple", "apple", "peer", "peer", "peer", "peer")
#' ordinal = c(4, 3, 1, 6, 5, 7, 2)
#' es_vargha_delaney_a(binary, ordinal, categories=c("peer", "apple"))
#' 
#' 
#' @export
es_vargha_delaney_a <- function(catField, ordField, categories=NULL, levels=NULL){
  
  #remove rows with missing values
  df = data.frame(ordField, catField)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  #replace the ordinal values if levels is provided
  if (!is.null(levels)){
    df$score = factor(df$score, ordered = TRUE, levels = levels)        
  }
  df$score = as.numeric(df$score)
  
  #the two categories
  if (!is.null(categories)){
    cat1 = categories[1]
    cat2 = categories[2]
  }
  else {
    cat1 = names(table(df$group))[1]
    cat2 = names(table(df$group))[2]
  }
  
  #seperate the scores for each category
  scoresCat1 = unname(unlist((subset(df, group == cat1)[1])))
  scoresCat2 = unname(unlist((subset(df, group == cat2)[1])))
  
  n1 = length(scoresCat1)
  n2 = length(scoresCat2)    
  n = n1 + n2
  
  #combine this into one long list
  allScores = c(scoresCat1, scoresCat2)
  
  #get the ranks
  allRanks = rank(allScores)
  
  #get the ranks per category
  cat1Ranks = allRanks[1:n1]
  cat2Ranks = allRanks[(n1+1):n]
  
  r1 = sum(cat1Ranks)
  r2 = sum(cat2Ranks)
  
  a1 = 1/n2 * (r1/n1 - (n1 + 1)/2)
  a2 = 1/n1 * (r2/n1 - (n2 + 1)/2)  
  
  results <- data.frame(a1, a2)
  colnames(results)<-c(paste("A-", cat1),paste("A-", cat2))
  
  return(results)  
}
