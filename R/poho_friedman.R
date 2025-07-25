#' Post-Hoc Tests for a Friedman Test
#' @description 
#' A post-hoc test after a Friedman test can be used to determine which variables differ significantly.
#' 
#' This function provides three options: Dunn, Conover, and Nemenyi.
#' 
#' @param data dataframe. A column for each variable
#' @param levels vector, optional. Indication of what the levels are in order
#' @param method string, optional. Post-Hoc method to use. Either "dunn" (default), "conover", "nemenyi"
#' @param ties boolean, optional. Apply a ties correction. Default is True
#' 
#' @returns
#' res : dataframe with the following columns
#' \item{field 1}{label of first column in pair}
#' \item{field 2}{label of second column in pair}
#' \item{n}{sample size}
#' \item{statistic}{test statistic used}
#' \item{df}{degrees of freedom (if applicable)}
#' \item{p-value}{the p-value (significance)}
#' \item{adj. p-value}{Bonferroni adjusted p-value}
#' 
#' @details
#' **Conover**
#' 
#' Bartz-Beielstein et al. (2010, p. 319) attributes this to Conover (1999) (but also seen sites refering to Conover (1980), just different editions of the book) and uses as a formula:
#' \deqn{t_{1,2} = \frac{\left|R_1 - R_2\right|}{SE}}
#' \deqn{SE = \sqrt{\frac{2\times n\times\left(1-\frac{\chi_F^2}{n\times\left(k-1\right)}\right)\times\left(\sum_{i=1}^n\sum_{j=1}^k r_{i,j}^2 - \frac{n\times k\times\left(k+1\right)^2}{4}\right)}{\left(n-1\right)\times\left(k-1\right)}}}
#' 
#' With:
#' \deqn{R_j = \sum_{i=1}^n r_{i,j}}
#' 
#' In the original source it mentions \eqn{2\times k} in \eqn{SE} instead of \eqn{2\times n}, this was indeed an error (pers. comm. with Conover).
#' 
#' Gnambs (n.d.) and BrightStat (n.d.) show a different formula, that gives the same result, and is the one the function uses:
#' \deqn{SE = \sqrt{\frac{2\times\left(\sum_{i=1}^n\sum_{j=1}^k r_{i,j}^2 - \sum_{j=1}^k R_j^2\right)}{\left(n-1\right)\times\left(k-1\right)}}}
#' 
#' The significance is then determined using:
#' \deqn{sig. = 2\times\left(1 - T\left(\left|t_{i,2}\right|, df\right)\right)}
#' 
#' Note that in the calculation \eqn{SE} is determined using all ranks, including those not in the comparison.
#' 
#' **Nemenyi**
#' 
#' Pohlert (2016, p. 15) shows the formula from Nemenyi (1963) as well as in Demsar (2006, pp. 11-12):
#' \deqn{q_{1,2} = \frac{\left|R_1 - R_2\right|}{\sqrt{\frac{k\times\left(k+1\right)}{6\times n}}}\times\sqrt{2}}
#' \deqn{df = n - k}
#' 
#' This follows then a studentized range distribution with:
#' \deqn{sig. = 1 - Q\left(q_{1,2}, k, df\right)}
#' 
#' **Dunn**
#' 
#' Benavoli et. al (2016, pp. 2-3) and IBM SPSS (2021, p. 814):
#' \deqn{z_{1,2} = \frac{\left|R_1 - R_2\right|}{SE}}
#' \deqn{SE = \sqrt{\frac{k\times\left(k+1\right)}{6\times n}}}
#' 
#' This follows a standard normal distribution:
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_{i,2}\right|\right)\right)}
#' 
#' **Bonferroni adjustment**
#' 
#' The Bonferroni adjustment is done using:
#' \deqn{sig._{adj} = \min \left(sig. \times n_c, 1\right)}
#' \deqn{n_c = \frac{k\times\left(k-1\right)}{2}}
#' 
#' *Symbols Used*
#' \itemize{
#' \item \eqn{n}, the number of cases
#' \item \eqn{k}, the number of variables
#' \item \eqn{r_{i,j}}, the rank of case i, in variable j. The ranks are determined for each case.
#' \item \eqn{\Phi\left(\dots\right)}, the standard normal cumulative distribution function.
#' \item \eqn{Q\left(\dots\right)}, the studentized range distribution cumulative distribution function.
#' \item \eqn{T\left(\dots\right)}, the Student t cumulative distribution function.
#' \item \eqn{n_c}, the number of comparisons (pairs)
#' }
#' 
#' @references
#' Benavoli, A., Corani, G., & Mangili, F. (2016). Should we really use post-hoc tests based on mean-ranks? *Journal of Machine Learning Research, 17*, 1-10. doi:10.48550/ARXIV.1505.02288
#' 
#' BrightStat. (n.d.). Friedman test. BrightStat. Retrieved November 5, 2023, from https://secure.brightstat.com/index.php?p=c&d=1&c=2&i=9
#' 
#' Conover, W. J. (1980). *Practical nonparametric statistics* (2nd ed.). Wiley.
#' 
#' Demsar, J. (2006). Statistical comparisons of classifiers over multiple data sets. *The Journal of Machine Learning Research, 7*, 1-30. doi:10.5555/1248547.1248548
#' 
#' Gnambs, T. (n.d.). SPSS Friedman. http://timo.gnambs.at/sites/default/files/spss_friedmanph.sps
#' 
#' IBM. (2021). IBM SPSS Statistics Algorithms. IBM.
#' 
#' Nemenyi, P. (1963). *Distribution-free Multiple Comparisons*. Princeton University.
#' 
#' Pohlert, T. (2016). The pairwise multiple comparison of mean ranks package (PMCMR). https://cran.r-hub.io/web/packages/PMCMR/vignettes/PMCMR.pdf
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ph_friedman <- function(data, levels=NULL, method="dunn", ties=TRUE){
  
  #remove missing values
  df = na.omit(data)
  n = nrow(df)
  k = length(colnames(df))
  
  if (!is.null(levels)){
    for (i in 1:k){
      myFieldOrd = factor(df[ ,i], ordered = TRUE, levels = levels)
      df[, i] = as.numeric(myFieldOrd)            
    }
  }
  
  #ranks per row
  for (i in 1:n){
    df[i, ] = rank(df[i, ])
  }
  
  
  rs = sum(df)
  rm = rs/(n*k)
  
  #Determine for each variable the average rank, and
  #the squared difference of this average with the overall average rank.
  rmj = colSums(df)/n
  rs2 = sum(df**2)
  sst = n*sum((rmj - rm)**2)
  sse = sum((df-rm)**2)/(n*k)*k/(k-1)
  
  if (ties){
    # Hollander and Wolfe (1999, p. 274) rewritten:
    qadj = sst / sse}
  else{
    # Friedman (1937, p. 679):
    qadj = 12 / (n * k * (k + 1)) * rs2 - 3 * n * (k + 1)}
  
  ncomp = k*(k - 1)/2
  #the pairwise comparisons
  varNames = colnames(df)
  resRow=1
  res = data.frame(matrix(nrow = ncomp, ncol = 7))
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      cat1 = varNames[i]
      cat2 = varNames[j]
      res[resRow, 1] = cat1
      res[resRow, 2] = cat2
      res[resRow, 3] = n  
      resRow = resRow + 1  
      
      
      
    }
  }
  
  useRow=1
  if (method=="dunn"){
    se = (k * (k + 1) / (6 * n))**0.5
    df = n - k  
    for (i in 1:(k-1)){
      for (j in (i+1):k){
        z = (rmj[i] - rmj[j])/se
        res[useRow,4] = z
        res[useRow,5] = NA
        res[useRow,6] = 2 * (1 - pnorm(abs(z))) 
        
        useRow = useRow + 1
      }
    }
  }
  
  else if (method == "conover"){
    r2 = sum((rmj * n)**2)
    se = (2 * (n * rs2 - r2) / ((n - 1) * (k - 1)))**0.5
    df = (n - 1) * (k - 1)
    for (i in 1:(k-1)){
      for (j in (i+1):k){
        tVal = (rmj[i] - rmj[j]) * n / se
        res[useRow,4] = tVal
        res[useRow,5] = df
        res[useRow,6] = 2 * (1 - pt(abs(tVal), df)) 
        
        useRow = useRow + 1
      }
    }
  }
  
  else if (method == "nemenyi"){
    se = (k * (k + 1) / (6 * n))**0.5
    df = n - k
    for (i in 1:(k-1)){
      for (j in (i+1):k){
        q = (rmj[i] - rmj[j]) / se * (2**0.5)
        res[useRow,4] = q
        res[useRow,5] = df
        res[useRow,6] = 1 - ptukey(abs(q), k, df, 1)
        
        useRow = useRow + 1
      }
    }
  }
  
  useRow=1
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      res[useRow,7] = res[useRow, 6]*ncomp
      if (res[useRow,7] > 1){
        res[useRow,7] = 1}
      useRow = useRow + 1
    }
  }
  
  return (res)
}



