#' Rank Biserial Correlation
#' @description 
#' This function will calculate Rank biserial correlation coefficient (independent-samples)
#' 
#' Cureton (1956) was perhaps the first to mention this term and provided a formula. His formula actually yields the same result as Goodman-Kruskal gamma (Goodman & Kruskal, 1954). Glass (1965; 1966) also developed a formula, but only for cases when there are no ties between the two categories. His formula will yield the same result as Somers'd (1962) and Cliff delta (1993). Cureton (1968) responded to Glass and gave his formula in an alternative form. Willson (1976) showed the link with Cureton formula and the Mann-Whitney U statistic. For more details on this see the article from Rubia (2022).
#' 
#' @param catField A vector with the scores data
#' @param ordField A vector with the group data
#' @param categories optional vector with categories to use and order for the categorical field. Otherwise the first two found will be used.
#' @param levels optional vector with the labels of the ordinal field in order.
#' @param version optional the method to use to calculate rank-biserial correlation. Either "cureton" (default), "glass"
#' 
#' @returns
#' Rank Biserial Correlation
#' 
#' @details
#' If version='cureton', the formula from Cureton (1968, p. 68) is used:
#' \deqn{r_{rb} = \frac{\bar{R}_1 - \left(n + 1\right)/2}{n_2/2 - B/n_1}}
#' 
#' If version='glass', the formula from Glass (1965, p. 91; 1966, p. 626) is used:
#' \deqn{r_b = \frac{2\times\left(\bar{R}_1 - \bar{R}_2\right)}{n}}
#' With:
#' \deqn{B = \frac{\sum_{i=1}^c t_{i,1} \times t_{i,2}}{2}}
#' \deqn{\bar{R}_i=\frac{R_i}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{\bar{R}_i} the average of ranks in category i
#' \item \eqn{R_i} the sum of ranks in category i
#' \item \eqn{n} the total sample size
#' \item \eqn{n_i} the number of scores in category i
#' \item \eqn{t_{i,j}}, the i-th number of tied scores in j
#' }
#' 
#' If one category has two scores of 3 and the other has three scores of 3, then \eqn{t_{1,1} = 2, t_{1,2} = 3}, if the first category has also one score of 4 and the second has two scores of 4, then \eqn{t_{2,1} = 1, t_{2,2} = 2}, etc.
#' 
#' Cureton's version is the same as Goodman-Kruskal gamma, while Glass's version is the same as Somers' d (1962, p. 804) and Cliff Delta (1993, p. 495).
#' 
#' The rank biserial can be converted to a Cohen d (using the **es_convert()** function), and then the rules-of-thumb for Cohen d could be used (**th_cohen_d()**)
#' 
#' @seealso 
#' \code{\link{es_convert}}, to convert to Cohen d, use `fr="rb", to="cohend"`.
#' 
#' \code{\link{th_cohen_d}}, rules of thumb for Cohen d
#' 
#' @references 
#' Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer ordinal questions. *Psychological Bulletin, 114*(3), 494-509. https://doi.org/10.1037/0033-2909.114.3.494
#' 
#' Cureton, E. E. (1956). Rank-biserial correlation. *Psychometrika, 21*(3), 287-290. https://doi.org/10.1007/BF02289138
#' 
#' Glass, G. V. (1966). Note on rank biserial correlation. *Educational and Psychological Measurement, 26*(3), 623-631. https://doi.org/10.1177/001316446602600307
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' myLevels = c('Not scientific at all', 'Not too scientific', 'Pretty scientific', 'Very scientific')
#' r_rank_biserial_is(df1['sex'], df1['accntsci'], levels=myLevels)
#' 
#' #Example 2: vectors
#' binary = c("apple", "apple", "apple", "peer", "peer", "peer", "peer")
#' ordinal = c(4, 3, 1, 6, 5, 7, 2)
#' r_rank_biserial_is(binary, ordinal, categories=c("peer", "apple"))
#' 
#' @export
r_rank_biserial_is <- function(catField, ordField, categories=NULL, levels=NULL, version="cureton"){
  
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
    cat1 = names(table(df[ ,2]))[1]
    cat2 = names(table(df[ ,2]))[2]
  }
  
  #seperate the scores for each category
  scoresCat1 = unname(unlist((subset(df, df[ ,2] == cat1)[1])))
  scoresCat2 = unname(unlist((subset(df, df[ ,2] == cat2)[1])))
  
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
  
  r1Avg = r1/n1
  r2Avg = r2/n2
  
  if (version=='glass'){
    rb = 2*(r1Avg - r2Avg)/n}
  else if (version=='cureton'){
    #determine bracket ties
    b = 0
    for (i in unique(cat1Ranks)) {
      b <- b + sum(cat2Ranks == i) * sum(cat1Ranks == i)
    }
    # rb using Cureton
    rb = (r1Avg - (n + 1) / 2) / (n2 / 2 - (b / 2) / n1)
  }
  return(rb)  
}



