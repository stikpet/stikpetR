#' Paired Samples Wilcoxon Signed Rank Test
#' 
#' @param ord1 the scores on the first variable
#' @param ord2 the scores on the second variable
#' @param appr c("wilcoxon", "exact", "imanz", "imant") optional which method to use for approximation (default is "wilcoxon")
#' @param noDiff c("wilcoxon", "pratt", "zsplit") optional method to deal with scores equal on both variables (default is "wilcoxon")
#' @param ties optional boolean to use a tie correction (default is True)
#' @param cc optional boolean to use a continuity correction (default is False)
#' @returns 
#' A dataframe with:
#' \item{W}{the Wilcoxon W}
#' \item{statistic}{the test statistic}
#' \item{df}{degrees of freedom (only applicable for Iman t approximation)}
#' \item{pValue}{the significance (p-value)}
#' \item{testUsed}{the test that was used}
#' 
#' @details 
#' The unadjusted test statistic is given by:
#' \deqn{W=\min\left(W_{neg}, W_{pos}\right)}
#' With:
#' \deqn{W_{pos} = \sum_{i=1}^n \begin{cases} r_i & \text{ if } d_i > 0 \\ 0 & \text{ if } d_i \le 0 \end{cases}}
#' \deqn{W_{neg} = \sum_{i=1}^n \begin{cases} r_i & \text{ if } d_i < 0 \\ 0 & \text{ if } d_i \ge 0 \end{cases}}
#' \deqn{d_i = x_i - y_i}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n} the number of scores (equal for each variable)
#' \item \eqn{x_i} the i-th score on the first variable
#' \item \eqn{y_i} the i-th score on the second variable
#' \item \eqn{r_i} the i-th rank of the absolute differences (\eqn{d_i})
#' \item \eqn{W_{pos}} is the number of ranks with a positive difference
#' \item \eqn{W_{neg}} is the number of ranks with a negative difference
#' }
#' 
#' The distribution and test for \eqn{W} can now be performed the same way as
#' for the one-sample case. See *ts_wilcoxon_os()* for details on the calculations.
#' The \eqn{d_i} scores are now the one-sample, and the hypothesized median would be 0.
#' 
#' **Alternatives**
#' 
#' *R's stats library *
#' 
#' wilcox.test(ord1, ord2, paired=TRUE, exact=FALSE, correct=TRUE)
#' 
#' wilcox.test(ord1, ord2, paired=TRUE, exact=FALSE, correct=FALSE)
#' 
#' *library(coin)*
#' 
#' wilcoxsign_test(ord1 ~ ord2, zero.method = "Wilcoxon")
#' 
#' wilcoxsign_test(ord1 ~ ord2, zero.method = "Pratt")
#' 
#' @examples 
#' ord1 = c(50,45,33,22,99,79,4,36,62,51,27,15,26,83,86)
#' ord2 = c(47,45,31,24,78,76,13,46,45,44,23,14,34,79,81)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="wilcoxon", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="wilcoxon", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="wilcoxon", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="wilcoxon", ties = FALSE, cc = TRUE)
#' 
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="pratt", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="pratt", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="pratt", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="pratt", ties = FALSE, cc = TRUE)
#' 
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="zsplit", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="zsplit", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="zsplit", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="wilcoxon", noDiff="zsplit", ties = FALSE, cc = TRUE)
#' 
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="wilcoxon", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="wilcoxon", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="wilcoxon", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="wilcoxon", ties = FALSE, cc = TRUE)
#' 
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="pratt", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="pratt", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="pratt", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="pratt", ties = FALSE, cc = TRUE)
#' 
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="zsplit", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="zsplit", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="zsplit", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imanz", noDiff="zsplit", ties = FALSE, cc = TRUE)
#' 
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="wilcoxon", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="wilcoxon", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="wilcoxon", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="wilcoxon", ties = FALSE, cc = TRUE)
#' 
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="pratt", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="pratt", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="pratt", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="pratt", ties = FALSE, cc = TRUE)
#' 
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="zsplit", ties = TRUE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="zsplit", ties = TRUE, cc = TRUE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="zsplit", ties = FALSE, cc = FALSE)
#' ts_wilcoxon_ps(ord1, ord2, appr="imant", noDiff="zsplit", ties = FALSE, cc = TRUE)
#' 
#' ord1 = c(5, 8, 6, 3, 9, 9, 18, 13, 16)
#' ord2 = c(4, 11, 8, 7, 3, 16, 5, 18, 7)
#' ts_wilcoxon_ps(ord1, ord2, appr="exact")
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods. *Biometrics Bulletin, 1*(6), 80. https://doi.org/10.2307/3001968
#'  
#' @export
ts_wilcoxon_ps <- function(ord1, ord2,
                           appr=c("wilcoxon", "exact", "imant", "imanz"), 
                           noDiff=c("wilcoxon", "zsplit", "pratt"), 
                           ties = TRUE, 
                           cc = FALSE){
  
  if (length(appr)>1) {
    appr="wilcoxon"
  }
  
  if (length(noDiff)>1) {
    noDiff="wilcoxon"
  }
  
  def="n.a."
  
  #remove missing values
  df = na.omit(data.frame(ord1, ord2))
  
  #calculate differences and absolute differences
  df$d = df$ord1 - df$ord2
  df$dAbs = abs(df$d)
  
  #remove zero differences if noDiff is set to 'wilcoxon'
  if (noDiff=="wilcoxon") {
    df = df[df$dAbs != 0,]
    
  }
  #rank the absolute differences
  df$rank = rank(df$dAbs)
  posSum = sum(df$rank[df$d > 0])
  negSum = sum(df$rank[df$d < 0])
  nr = nrow(df)
  
  W = min(posSum, negSum)

  if (appr=="exact") {
    #check if ties existed
    if (max(table(df$rank)) > 1) {
      return ("ties exist, cannot compute exact method") }
    else {
      nr = nrow(df)
      pVal = psignrank(W, nr)*2
      statistic = W
      testUsed = "paired-sample Wilcoxon signed rank exact test"
    }
  }
  
  else{
    #add half the equal to median ranks if zsplit is used
    nD0 = sum(df$dAbs==0)
    if (noDiff == "zsplit") {
      W = W + sum(df$rank[df$d == 0])/2}
    
    rAvg = nr*(nr + 1)/4
    s2 = nr * (nr+1) * (2*nr + 1)/24
    
    #adjust if Pratt method is used
    if (noDiff == "pratt") {
      #normal approximation adjustment based on Cureton (1967)
      s2 = s2 - nD0 * (nD0 + 1) * (2 * nD0 + 1) / 24
      rAvg = (nr * (nr + 1) - nD0 * (nD0 + 1)) / 4
    }
    
    #the ties correction
    t = 0
    if (ties) {
      #remove ranks of scores equal to hypMed for Pratt (wilcoxon already removed)
      if (noDiff == "pratt") {
        df = df[df$d != 0,]
      }
      
      for (i in table(df$rank)){
        t = t + (i^3 - i)
      }
      t = t/48
      s2 = s2 - t
    }
    
    se = sqrt(s2)
    
    num = abs(W - rAvg)
    #apply continuity correction if needed
    if (cc) {
      num = num - 0.5}
    
    
    if (appr=="imant") {
      statistic = num / sqrt((s2 * nr - (W - rAvg) ^ 2) / (nr - 1))
      def = nr - 1
      pVal = 2 * (1 - pt(abs(statistic), def))}
    
    else{
      statistic = num / se
      if (appr == "imanz"){
        statistic = statistic / 2 * (1 + sqrt((nr - 1) / (nr - statistic ^ 2)))}
      pVal = 2 * (1 - pnorm(abs(statistic)))  
      
    }
    
    testUsed = "paired Wilcoxon signed rank test"
    
    if (ties && cc){
      testUsed = paste0(testUsed, ", with ties and continuity correction")}
    else if (ties){
      testUsed = paste0(testUsed, ", with ties correction")}
    else if (cc){
      testUsed = paste0(testUsed, ", with continuity correction")}
    
    if (appr == "imant") {
      testUsed = paste0(testUsed, ", using Iman (1974) t approximation")}
    else if (appr == "imanz"){
      testUsed = paste0(testUsed, ", using Iman (1974) z approximation")}
    
    if (noDiff == "pratt"){
      testUsed = paste0(testUsed, ", Pratt method for equal to hyp. med. (inc. Cureton adjustment for normal approximation)")}
    else if (noDiff == "zsplit"){
      testUsed = paste0(testUsed, ", z-split method for equal to hyp. med.")}
    
  }
  
  df = def
  testResults <- data.frame(W, statistic, df, pVal, testUsed)
  
  return (testResults)
  
}
