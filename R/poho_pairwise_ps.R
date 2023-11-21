#' Post-Hoc Pairwise Paired Samples Tests
#' @description 
#' This function simply performs pairwise paired tests: Sign, Wilcoxon and Trinomial. It then adds a Bonferroni correction.
#' 
#' These could be used with a Friedman test, but other post-hoc tests are also available in the ph_friedman() function (a Dunn, Nemenyi and Conover test).
#' 
#' @param data dataframe. A column for each variable
#' @param levels vector, optional. Indication of what the levels are in order
#' @param method string, optional. Post-Hoc method to use. Either "dunn" (default), "conover", "nemenyi"
#' @param test string, optional. Test to use in pairwise comparisons. Either "sign" (details), "wilcoxon", "trinomial".
#' @param appr string, optional. Option for sign and wilcoxon test. Default for wilcoxon is wilcoxon, for sign is appr. Either "exact", "appr", "wilcoxon", "imanz", "imant"
#' @param noDiff string, optional. Method to deal with scores equal to mu. Either "wilcoxon" (default), "pratt", "zsplit". Only applies if test="wilcoxon"
#' @param ties boolean, optional. Apply a ties correction. Default is True
#' @param cc boolean, optional. use a continuity correction. Default is False. Only applies if test="wilcoxon"
#' 
#' @returns 
#' res, a dataframe with the test results and:
#' \item{var 1}{the name of the first variable in the pair}
#' \item{var 2}{the name of the second variable in the pair}
#' \item{adj. p-value}{the Bonferroni adjusted p-value}
#' 
#' @details
#' This function creates each possible pair of the variables (columns) and then uses the requested paired samples test.
#' 
#' See for the calculations:
#' \itemize{
#' \item Sign test -> ts_sign_ps()
#' \item Wilcoxon signed rank test -> ts_wilcoxon_ps()
#' \item Trinomial test -> ts_trinomial_ps()
#' }
#' 
#' The Bonferroni adjustment is done using:
#' \deqn{sig._{adj} = \min \left(sig. \times n_c, 1\right)}
#' \deqn{n_c = \frac{k\times\left(k-1\right)}{2}}
#' 
#' Where \eqn{n_c} is the number of comparisons (pairs)
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_pairwise_ps <- function(data, levels=NULL, test = "sign", 
                           appr = "wilcoxon", 
                           noDiff = "wilcoxon", 
                           ties = TRUE, 
                           cc = FALSE){
  
  #remove missing values
  df = na.omit(data)
  n = nrow(df)
  k = length(colnames(df))
  
  #the pairwise comparisons
  ncomp = k*(k - 1)/2
  colNames = colnames(df)
  useRow=1
  
  if (test=="sign"){
    resCol = 7
    resLbls = c("var 1", "var 2", "n pos", "n neg", "statistic", "p-value", "adj. p-value")
  }
  if (test=="wilcoxon"){
    resCol = 9
    resLbls = c("var 1", "var 2", "nr", "mu", "W", "statistic", "df", "p-value", "adj. p-value")
  }
  if (test=="trinomial"){
    resCol = 7
    resLbls = c("var 1", "var 2", "n pos", "n neg", "n 0", "p-value", "adj. p-value")
  }
  
  res = data.frame(matrix(nrow = ncomp, ncol = resCol))
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      cat1 = colNames[i]
      cat2 = colNames[j]
      res[useRow,1] = cat1
      res[useRow,2] = cat2 
      
      if (test=="sign"){
        testRes = ts_sign_ps(df[, i], df[, j], levels=levels, method=appr)
        res[useRow,3] = testRes[1,1]
        res[useRow,4] = testRes[1,2]
        res[useRow,5] = testRes[1,3]
        res[useRow,6] = testRes[1,4]
        res[useRow,7] = testRes[1,4]*ncomp
        if (res[useRow,7] > 1){
          res[useRow,7] = 1
        }
      }
      
      else if (test=="wilcoxon"){
        testRes = ts_wilcoxon_ps(df[, i], df[, j], levels=levels, appr=appr, noDiff=noDiff, ties=ties, cc=cc)
        res[useRow,3] = testRes[1,1]
        res[useRow,4] = testRes[1,2]
        res[useRow,5] = testRes[1,3]
        res[useRow,6] = testRes[1,4]
        res[useRow,7] = testRes[1,5]
        res[useRow,8] = testRes[1,6]
        res[useRow,9] = testRes[1,6]*ncomp
        if (res[useRow,9] > 1){
          res[useRow,9] = 1
        }
      }
      
      else if (test=="trinomial"){
        testRes = ts_trinomial_ps(df[, i], df[, j], levels=levels)
        res[useRow,3] = testRes[1,1]
        res[useRow,4] = testRes[1,2]
        res[useRow,5] = testRes[1,3]
        res[useRow,6] = testRes[1,4]
        res[useRow,7] = testRes[1,4]*ncomp
        if (res[useRow,7] > 1){
          res[useRow,7] = 1
        }
      }
      
      useRow=useRow+1
      
    }
  }
  
  colnames(res) = resLbls
  return (res)
  
}