#' Friedman Test
#' @description 
#' A test to determine if any of the variables has a significant different average ranking than any of the others.
#' 
#' It is a paired-samples version of a Kruskal-Wallis test. If the p-value is below a pre-defined threshold (usually 0.05) it indicates at least one variable (column) is different than another.
#' 
#' @param data dataframe. A column for each variable
#' @param levels vector, optional. Indication of what the levels are in order
#' @param ties boolean, optional. Apply a ties correction. Default is True
#' @param dist string, optional. Distribution to use. Either "chi" (default), "f", "normal"
#' 
#' @returns
#' res : dataframe with the following columns
#' \item{n}{sample size}
#' \item{statistic}{test statistic used}
#' \item{df, df1, df2}{degrees of freedom (if applicable)}
#' \item{p-value}{the p-value (significance)}
#' 
#' @details 
#' The formula used in case of no ties (Friedman, 1937, p. 679):
#' \deqn{\chi_F^2 = \left(\frac{12}{n\times k\times\left(k+1\right)}\times\sum_{j=1}^k R_j^2\right)-3\times n\times\left(k+1\right)}
#' \deqn{df = k - 1}
#' 
#' With:
#' \deqn{R_j = \sum_{i=1}^n r_{i,j}}
#' 
#' In case a ties correction is used (Hollander & Wolfe, 1999, p. 274):
#' \deqn{\chi_{Fadj}^2 = \frac{12\times \sum_{j=1}^k R_j^2 - 3\times n^2\times\left(k +1\right)^2}{n\times\left(k+1\right) - \frac{\left(\sum t_{i,j}^3\right)-k}{k-1}}}
#' 
#' The ties correction used by IBM SPSS (2021, p. 811) will give the same result:
#' \deqn{\chi_{Fadj}^2 = \frac{\chi_F^2}{1 - \frac{\sum t_{i,j}^3 - t_{i,j}}{n\times\left(k^3-k\right)}}}
#' 
#' The function uses more of a one-way ANOVA approach in case of ties, but then on the ranks. It leads to the same result:
#' \deqn{\chi_{Fadj}^2 = \frac{n\times\sum_{j=1}^k\left(\bar{r}_j -\bar{r}\right)^2}{\left(\frac{\sum_{j=1}^k\sum_{i=1}^n \left(r_{i,j}-\bar{r}\right)^2}{n\times\left(k-1\right)}\right)}}
#' 
#' With:
#' \deqn{\bar{r}_j = \frac{R_j}{n}}
#' \deqn{\bar{r} = \frac{\sum_{j=1}^k R_j}{n\times k} = \frac{n\times\left(k+1\right)}{2}}
#' 
#' The significance is then determined using:
#' \deqn{sig. = 1 - \chi^2\left(\chi_F^2, df\right)}
#' 
#' A normal distribution approximation was proposed by Friedman (1937, p. 695; 1939, p. 109):
#' \deqn{z_F = \frac{\chi_F^2-\left(k-1\right)}{\sqrt{2\times\frac{n-1}{n}\times\left(k-1\right)}}}
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_F\right|\right)\right)}
#' 
#' And an F distribution by Iman and Davenport (1980, p. 573):
#' \deqn{F_F = \frac{\left(n-1\right)\times\chi_F^2}{n\times\left(k-1\right)-\chi_F^2}}
#' \deqn{df_1 = k - 1}
#' \deqn{df_2 = \left(k-1\right)\times\left(n - 1\right)}
#' \deqn{sig. = 1 - F\left(F_{F}, df_1, df_2\right)}
#' 
#' Some might refer to Conover for this F-distribution, but in Conover (1980, p. 300) it seems Conover credits Iman and Davenport himself. 
#' 
#' *Symbols Used*
#' 
#' \itemize{
#' \item \eqn{n}, the number of cases
#' \item \eqn{k}, the number of variables
#' \item \eqn{r_{i,j}}, the rank of case i, in variable j. The ranks are determined for each case.
#' \item \eqn{t_{i,j}}, the frequency of unique rank j, in case i. For each row the frequencies of each rank is determined in the calculations.
#' }
#' 
#' @references
#' Conover, W. J. (1980). *Practical nonparametric statistics* (2nd ed.). Wiley.
#' 
#' Friedman, M. (1937). The use of ranks to avoid the assumption of normality implicit in the analysis of variance. *Journal of the American Statistical Association, 32*(200), 675–701. doi:10.2307/2279372
#' 
#' Friedman, M. (1939). A correction. *Journal of the American Statistical Association, 34*(205), 109–109. doi:10.1080/01621459.1939.10502372
#' 
#' Hollander, M., & Wolfe, D. A. (1999). *Nonparametric statistical methods* (2nd ed.). Wiley.
#' 
#' IBM. (2021). IBM SPSS Statistics Algorithms. IBM.
#' 
#' Iman, R., & Davenport, J. (1980). Approximations of the critical region of the Friedman statistic. *Communications in Statistics-Theory and Methods, 9*, 571–595.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
ts_friedman <- function(data, levels=NULL, ties=TRUE, dist="chi"){
  
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
  rs2 = sum((rmj*n)**2)
  sst = n*sum((rmj - rm)**2)
  sse = sum((df-rm)**2)/(n*k)*k/(k-1)
  
  if (ties){
    # Hollander and Wolfe (1999, p. 274) rewritten:
    qadj = sst / sse}
  else{
    # Friedman (1937, p. 679):
    qadj = 12 / (n * k * (k + 1)) * rs2 - 3 * n * (k + 1)}
  
  df = k - 1
  
  if (dist=="f"){
    # Iman-Davenport F Distribution (1980, p. 573)
    fVal = (n - 1) * qadj / (n * (k - 1) - qadj)
    df1 = df
    df2 = (k - 1) * (n - 1)
    p = 1 - pf(fVal, df1, df2)
    
    res = data.frame(n, fVal, df1, df2, p)
    colnames(res) = c("n", "statistic", "df1", "df2", "p-value")
  }
  
  else if (dist=="normal"){        
    # Friedman Normal Distribution (1937, p. 695; 1939, p. 109):
    z = (qadj - (k - 1)) / (2 * (n - 1) / n * (k - 1))**0.5
    p = 2 * (1 - pnorm(abs(z))) 
    res = data.frame(n, z, p)
    colnames(res) = c("n", "statistic", "p-value")
  }
  
  else{
    #Friedman Chi-Square Distribution
    p = 1 - pchisq(qadj, df)
    res = data.frame(n, qadj, df, p)
    colnames(res) = c("n", "statistic", "df", "p-value")
  }
  
  return (res)
  
}