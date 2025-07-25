#' Mann-Whitney U Test
#' @description
#' The Mann-Whitney U and Wilcoxon Rank Sum test are the same. Mann and Whitney simply expanded on the ideas from Wilcoxon.
#' 
#' The test will compare the distribution of ranks between two categories. The assumption is that the two categories have the same mean rank (which often is stated simplified as having the same median in the population).
#' 
#' @param catField A vector or dataframe with the group data
#' @param ordField A vector or dataframe with the scores data
#' @param categories : optional list with the two categories to use from catField. If not set the first two found will be used
#' @param levels optional list with the scores in order
#' @param method c("exact", "appr") exact method or normal approximation
#' @param cc boolean to indicate the use of a continuity correction
#' 
#' @returns
#' A dataframe with:
#' \item{n}{the sample size}
#' \item{U1}{the Mann-Whitney U score of the first category}
#' \item{U2}{the Mann-Whitney U score of the second category}
#' \item{statistic}{test statistic}
#' \item{pValue}{significance (p-value)}
#' \item{test}{description of the test used}
#' 
#' @details
#' The formula used is (Mann & Whitney, 1947, p. 51):
#' \deqn{U_i = R_i - \frac{n_i\times\left(n_i + 1\right)}{2}}
#' With:
#' \deqn{R_i = \sum_{j=1}^{n_i} r_{i,j}}
#' 
#' For an approximation the following is used:
#' \deqn{sig. = 2\times\left(1 - Z\left(z\right)\right)}
#' With:
#' \deqn{z = \frac{U_i - \frac{n_1\times n_2}{2}}{SE}}
#' \deqn{SE = \sqrt{\frac{n_1\times n_2}{n\times\left(n - 1\right)}\times\left(\frac{n^3 - n}{12} - \sum_i T_i\right)}}
#' \deqn{T_i = \frac{t_i^3 - t_i}{12}}
#' \deqn{n = n_1 + n_2}
#' 
#' If a continuity correction is used the z-value is calculated using:
#' \deqn{z_{cc} = z - \frac{0.5}{SE}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_i} the sample size of category i
#' \item \eqn{n} the total sample size
#' \item \eqn{r_{i,j}} the j-th rank of category i
#' }
#' 
#' The ties correction (\eqn{T}) can be found in Lehmann and D'Abrera (1975, p. 20)
#' 
#' For the exact distribution the Mann-Whitney-Wilcoxon distribution is used, from the **pwilcox()** function from R.
#' 
#' Wilcoxon (1945) had developed this test earlier for the case when both categories have the same sample size, and Mann and Whitney expanded on this.
#' 
#' @section Before, After and Alternatives:
#' Before running the test you might first want to get an impression using a cross table using \code{\link{tab_cross}}, or a stacked bar chart using \code{\link{vi_bar_stacked_multiple}} for a visualisation.
#' 
#' After the test you might want an effect size measure, either \code{\link{es_common_language_is}} for the CLES, \code{\link{es_hodges_lehmann_is}} for Hodges-Lehmann, or \code{\link{r_rank_biserial_is}} for the (Glass) rank biserial (Cliff delta).
#' 
#' Independent samples tests for a binary vs ordinal variable include Brunner-Munzel (\code{\link{ts_brunner_munzel}}), C-square (\code{\link{ts_c_square}}), Cliff-Delta (\code{\link{ts_cliff_delta_is}}), Fligner-Policello (\code{\link{ts_fligner_policello}}), Mann-Whitney U (\code{\link{ts_mann_whitney}}), Mood-Median (\code{\link{ts_mood_median}})
#' 
#'
#' @references 
#' Lehmann, E. L., & D'Abrera, H. J. M. (1975). *Nonparametrics: Statistical methods based on ranks*. Holden-Day.
#' 
#' Mann, H. B., & Whitney, D. R. (1947). On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other. *The Annals of Mathematical Statistics, 18*(1), 50-60. https://doi.org/10.1214/aoms/1177730491
#' 
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods. *Biometrics Bulletin, 1*(6), 80. https://doi.org/10.2307/3001968
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' myLevels = c('Not scientific at all', 'Not too scientific', 'Pretty scientific', 'Very scientific')
#' ts_mann_whitney(df1['sex'], df1['accntsci'], levels=myLevels)
#' 
#' #Example 2: vectors
#' binary = c("apple", "apple", "apple", "peer", "peer", "peer", "peer")
#' ordinal = c(4, 3, 1, 6, 5, 7, 2)
#' ts_mann_whitney(binary, ordinal, categories=c("peer", "apple"))
#' 
#' @export
ts_mann_whitney <- function(catField, ordField, categories=NULL, levels=NULL, method="exact", cc=TRUE){
  
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
  df = na.omit(df)
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
  
  R1 = sum(cat1Ranks)
  R2 = sum(cat2Ranks)
  
  #The U statistics
  U1 = R1 - n1 * (n1 + 1) / 2
  U2 = R2 - n2 * (n2 + 1) / 2
  U = min(U1, U2)
  
  #The count of each rank
  counts = table(allRanks)
  
  #check if ties exist
  if (max(counts)>1 && method=="exact"){
    print("ties exist, swith to approximate")
    method="approx"
  }
  
  if (method=="exact") {
    testUsed = "Mann-Whitney U exact"
    if (U2==U){
      nf = n1
      n1 = n2
      n2 = nf
    }
    
    pValue = pwilcox(U, n1, n2)*2
    statistic = "n.a."
  }
  else{
    testUsed = "Mann-Whitney U normal approximation"
    t = 0
    for (i in counts) {
      t = t + i^3 - i
    }
    t = t/12
    
    se = sqrt(n1*n2/(n*(n-1))*((n^3 - n)/12 - t))
    z = (U - n1*n2/2)/se
    zabs = abs(z)
    
    if (cc) {
      zabs = zabs - 0.5/se
      testUsed = "Mann-Whitney U normal approximation, with continuity correction"
    }
    
    #still need abs since cc could make it negative
    pValue = 2*(1 - pnorm(abs(zabs)))
    statistic = zabs
  }
  
  results <- data.frame(n, U1, U2, statistic, pValue, testUsed)
  colnames(results)<-c("n", "U1", "U2", "statistic", "p-value", "test")
  
  return(results)  
}



