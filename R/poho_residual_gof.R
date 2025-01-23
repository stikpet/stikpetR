#' Post-Hoc Residuals Goodness-of-Fit Test
#' 
#' @description 
#' This function will perform a standardized residuals post-hoc test for each of the categories in a nominal field.
#' 
#' The unadjusted p-values and Bonferroni adjusted p-values are both determined.
#' 
#' @param data dataframe with scores
#' @param expCount optional dataframe with categories and expected counts
#' 
#' @returns
#' a dataframe with:
#' \item{category}{the label of the first category}
#' \item{z-statistic}{the standardized residuals}
#' \item{p-value}{the unadjusted significance}
#' \item{adj. p-value}{the adjusted significance}
#' 
#' @details
#' The formula used is:
#' \deqn{z = \frac{F_i - E_i}{\sqrt{E_i}}}
#' \deqn{sig = 2\times\left(1 - \Phi\left(\left|z\right|\right)\right)}
#' 
#' With:
#' \itemize{
#' \item \eqn{F_i}, the observed count for category $i$
#' \item \eqn{E_i}, the expected count for category $i$
#' \item \eqn{\Phi\left(\dots\right)}, the cumulative distribution function of the standard normal distribution
#' }
#' 
#' If no expected counts are provide it is assumed they are all equal for each category, i.e. \eqn{E_i = \frac{n}{k}}
#' 
#' The Bonferroni adjustment is calculated using:
#' \deqn{p_{adj} = \min \left(p \times k, 1\right)}
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_residual_gof <- function(data, expCount=NULL){
  data = na.omit(data)
  
  #the sample size n
  n = length(data)
  
  #determine the observed counts
  if (is.null(expCount)){
    #generate frequency table
    freqTable<-table(data)
    
    #number of categories to use (k)
    k = dim(freqTable)
    
    #rewrite frequency table in long format
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in 1:k){
      freq[nrow(freq) + 1,] = c(names(freqTable)[i], freqTable[i])
    }
  }
  
  else{
    #if expected counts are given
    freq = data.frame(matrix(nrow=0, ncol=2))
    for (i in expCount[,1]){
      freq[nrow(freq) + 1,] = c(i, sum(data==i))
    }
    
    #number of categories to use (k)
    k = nrow(expCount)    
  }
  
  n = sum(as.numeric(freq[,2]))  
  
  
  #the true expected counts
  if (is.null(expCount)){
    #assume all to be equal
    expC = rep(n/k,k)
    #number of expected counts is simply sample size
    nE = n  
  }
  else{
    #sum of the given expected counts
    nE = sum(expCount[,2])
    #check if categories match
    expC = c()
    for (i in 1:k){
      for (j in 1:k){
        if (expCount[i,1]==freq[j,1]){
          expC = append(expC, expCount[i,2]/nE * n)
        }
      }
    }
    
  }
  
  columns=c("category", "z-statistic","p-value", "adj. p-value")
  
  res = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(res) = columns
  adjFactor = k*(k - 1)/2
  for (i in 1:k){
    cat = freq[i,1]
    z = (as.numeric(freq[i,2]) - expC[i])/(expC[i]**0.5)
    sig = 2*(1 - pnorm(abs(z)))
    adjSig = min(sig*k, 1)
    res[i,] = c(cat, z, sig, adjSig)
  }
  
  return (res)
}