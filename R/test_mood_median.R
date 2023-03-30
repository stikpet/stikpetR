#' Mood Median Test
#' 
#' 
#' @param scores the variable with the scores
#' @param groups the variable with the groups
#' @param cc optional continuity correction to use
#' @returns 
#' Returns a dataframe with:statistic, df, pValue, minExp, propBelow5
#' \item{statistic}{the test statistic}
#' \item{df}{the degrees of freedom}
#' \item{pValue}{the two-tailed significance, a.k.a. p-value}
#' \item{minExp}{the minimum expected count}
#' \item{propBelow5}{the proportion of cells with expected count below 5}
#' 
#' @details 
#' This test simply counts how many scores fall above the median in each category and 
#' how many equal or below, and then perform a Pearson chi-square test on the table
#' 
#' The formula used is:
#' \deqn{\chi_M^2 = \sum_{i=1}^2 \sum_{j=1}^k \frac{\left(F_{i,j} - E_{i,j}\right)^2}{E_{i,j}}}
#' \deqn{df = k - 1}
#' \deqn{sig. = 1 - \chi^2\left(\chi_M^2, df\right)}
#' With:
#' \deqn{E_{i,j} = \frac{R_i\times C_j}{n}}
#' \deqn{R_i = \sum_{j=1}^k F_{i,j}}
#' \deqn{C_i = \sum_{i=1}^2 F_{i,j}}
#' \deqn{n = \sum_{i=1}^2 \sum_{j=1}^k F_{i,j}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{F_{1,j}} the number of scores is category j that are above the overall median
#' \item \eqn{F_{2,j}} the number of scores is category j that are equal or below the overall median
#' \item \eqn{E_{i,j}} the expected count for row i, column j
#' \item \eqn{R_i} the row total of row i
#' \item \eqn{C_i} the column total of column i
#' \item \eqn{n} the total sample size
#' \item \eqn{k} the number of categories
#' \item \eqn{\chi^2\left(\dots, \dots\right)} the cumulative distribution function of the chi-square distribution
#' }
#' 
#' The Yates correction (yates) is calculated using (Yates, 1934, p. 222):
#' 
#' Use instead of \eqn{F_{i,j}} the adjusted version defined by:
#' \deqn{F_{i,j}^\ast = \begin{cases} F_{i,j} - 0.5 & \text{ if } F_{i,j}>E_{i,j}  \\ F_{i,j} & \text{ if } F_{i,j}= E_{i,j}\\ F_{i,j} + 0.5 & \text{ if } F_{i,j}<E_{i,j} \end{cases}}
#' 
#' The Pearson correction (pearson) is calculated using (E.S. Pearson, 1947, p. 157):
#' \deqn{\chi_{PP}^2 = \chi_{P}^{2}\times\frac{n - 1}{n}}
#' 
#' The Williams correction (williams) is calculated using (Williams, 1976, p. 36):
#' \deqn{\chi_{PW}^2 = \frac{\chi_{P}^2}{q}}
#' With:
#' \deqn{q = 1 + \frac{\left(n\times\left(\sum_{i=1}^r \frac{1}{R_i}\right)-1\right) \times \left(n\times\left(\sum_{j=1}^c \frac{1}{C_j}\right)-1\right)}{6\times n\times df}}
#' 	
#' Often Mood (1950) is used as a source, but I found the explanation clearer in Brown and Mood (1951)
#' 
#' @references 
#' Brown, G. W., & Mood, A. M. (1951). On median tests for linear hypotheses. *Proceedings of the Second Berkeley Symposium on Mathematical Statistics and Probability*, 2, 159–167. 
#' 
#' Mood, A. M. (1950). *Introduction to the theory of statistics*. McGraw-Hill.
#' 
#' Pearson, E. S. (1947). The choice of statistical tests illustrated on the Interpretation of data classed in a 2 × 2 table. *Biometrika, 34*(1/2), 139–167. https://doi.org/10.2307/2332518
#' 
#' Williams, D. A. (1976). Improved likelihood ratio tests for complete contingency tables. *Biometrika, 63*(1), 33–37. https://doi.org/10.2307/2335081
#' 
#' Yates, F. (1934). Contingency tables involving small numbers and the chi square test. *Supplement to the Journal of the Royal Statistical Society, 1*(2), 217–235. https://doi.org/10.2307/2983604
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' scores = c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, NA, 3, 4, 2, 4, 2, 1, 3, 2, 2, 4, 1, 1, 3, 1, 2, 4, 1, 5, 4, 2, 3, 4, 1, 2, 5, 1, 1, 3, 3, 3, 1, 4, 3, 1, 1, 2, 3, 1)
#' groups = c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
#' ts_mood_median(scores, groups)
#' ts_mood_median(scores, groups, cc="yates")
#' ts_mood_median(scores, groups, cc="pearson")
#' ts_mood_median(scores, groups, cc="williams")
#' 
#' @export
ts_mood_median <- function(scores, groups, cc=c(NULL,"yates", "pearson", "williams")){
  
  if (length(cc)>1) {
    cc = cc[1]
  }
  
  datFrame = na.omit(data.frame(groups, scores))
  med = median(datFrame$scores)
  
  datTable = table(datFrame)
  nr = nrow(datTable) #the number of groups
  
  obs = matrix(1, nrow=2, ncol=nr)
  for (i in 1:nr) {
    obs[1,i] = sum(datFrame$scores[datFrame$groups == row.names(datTable)[i]] > med)
    obs[2,i] = sum(datFrame$scores[datFrame$groups == row.names(datTable)[i]] <= med)
  }
  
  R = rowSums(obs)
  C = colSums(obs)
  n = sum(R)
  
  expCount = matrix(1, nrow=2, ncol=nr)
  for (j in 1:nr) {
    for (i in 1:2){
      expCount[i,j] = R[i]*C[j]/n
    }
  }
  
  if (!is.null(cc) && cc=="yates") {
    #Yates correction
    for (j in 1:nr) {
      for (i in 1:2){
        if (obs[i,j] > expCount[i,j]) {
          obs[i,j] = obs[i,j] - 0.5
        }
        else if (obs[i,j] < expCount[i,j]) {
          obs[i,j] = obs[i,j] + 0.5
        }
      }
    }
  }
  
  chiVal = sum((obs - expCount)^2/expCount)
  
  if (!is.null(cc) && cc=="pearson") {
    #Pearson Correction
    chiVal = (n - 1)/n * chiVal
  }
  else if (cc=="williams"){
    #Williams Correction
    q = 1 + (n*sum(1/R)-1) * (n*sum(1/C)-1) / (6*n*(nr - 1)*(2 - 1))
    chiVal = chiVal/q
  }
  
  df = nr - 1
  
  pValue = 1 - pchisq(chiVal, df)
  
  minExp = min(expCount)
  propBelow5 = sum(expCount < 5)
  
  statistic = chiVal
  results <- data.frame(statistic, df, pValue, minExp, propBelow5)
  
  return(results)
  
}