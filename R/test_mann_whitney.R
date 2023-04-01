#' Mann-Whitney U Test
#' 
#' @param dataVar A vector with the scores data
#' @param groupVar A vector with the group data
#' @param method c("exact", "appr") exact method or normal approximation
#' @param corr boolean to indicate the use of a continuity correction
#' @return dataframe with U, the test statistic, p-value, and the test used
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
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references 
#' Lehmann, E. L., & D’Abrera, H. J. M. (1975). *Nonparametrics: Statistical methods based on ranks*. Holden-Day.
#' 
#' Mann, H. B., & Whitney, D. R. (1947). On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other. *The Annals of Mathematical Statistics, 18*(1), 50–60. https://doi.org/10.1214/aoms/1177730491
#' 
#' @examples 
#' scores = c(5, 12, 3, 4, 6, 1, 11, 13, NA)
#' groups = c("A","A","A","B","B","B","B", NA, "C")
#' ts_mann_whitney(scores, groups)
#' 
#' @export
ts_mann_whitney <- function(dataVar, groupVar, method="exact", corr=TRUE){
  
  #make sure data is numeric
  scores = as.numeric(dataVar)
  
  #remove rows with missing values
  df = data.frame(scores, groupVar)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  n = length(df$group)
  n1 = sum(df$group==df$group[1])
  n2 = n - n1
  
  rankScores = rank(df$score)
  
  R1 = sum(rankScores[df$group==df$group[1]])
  R2 = n*(n+1)/2 - R1
  
  U1 = n1*n2 + n1*(n1 + 1)/2 - R1
  U2 = n1*n2 + n2*(n2 + 1)/2 - R2
  
  U = min(U1, U2)
  
  ties=FALSE
  if (length(unique(df$score))!=n) {
    ties=TRUE
  }
  
  if (method=="exact" && ties) {
    print("ties present, cannot use exact method")
    method="approx"
  }
  
  if (method=="exact") {
    pVal = pwilcox(U, n1, n2)*2
    testUsed = "Mann-Whitney U exact"
    testResults <- data.frame(U, pVal, testUsed)
    
  }
  else{
    testUsed = "Mann-Whitney U normal approximation"
    t = 0
    freq = table(rankScores)
    for (i in freq) {
      t = t + i^3 - i
    }
    t = t/12
    
    se = sqrt(n1*n2/(n*(n-1))*((n^3 - n)/12 - t))
    z = (U - n1*n2/2)/se
    zabs = abs(z)
    
    if (corr) {
      zabs = zabs - 0.5/se
      testUsed = "Mann-Whitney U normal approximation, with continuity correction"
    }
    
    #still need abs since cc could make it negative
    pValue = 2*(1 - pnorm(abs(zabs)))
    
    statistic = z
    testResults <- data.frame(U, statistic, pValue, testUsed)
    
  }
  
  return(testResults)
  
}