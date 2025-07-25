#' Pearson Correlation Coefficient
#' @importFrom gsl hyperg_2F1
#' 
#' @description 
#' A measure of linear correlation. A -1 indicates a perfect negative linear correlation (i.e. a straight line going down, if the score in one field goes up, the other one goes down), a 0 would indicate no correlation, and a +1 a perfect positive linear correlation (i.e. a straight line going up, if the score in one field goes up, the other one goes up as well).
#' 
#' Various tests can be used to determine if the coefficient is significantly different from zero. See notes for details.
#' 
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param corr string, optional. Which adjustment to make if any (default is "none")
#' @param test string, optional. Which test to use (see details). Either "t" (default), "z"
#' 
#' @returns 
#' A dataframe with:
#' \item{r}{the Pearson Correlation Coefficient}
#' \item{statistic}{the test statistic}
#' \item{df}{degrees of freedom (only applicable for t test)}
#' \item{p-value}{the significance (p-value)}
#' 
#' @details 
#' This function makes use of the *hyperg_2F1* function from the *gsl* library.
#' 
#' The formula used (Pearson, 1896, p. 265):
#' \deqn{r = \frac{\sum_{i=1}^n \left(x_i - \bar{x}\right) \times \left(y_i - \bar{y}\right)}{SS_x\times SS_y}}
#' With:
#' \deqn{SS_x = \sum_{i=1}^n \left(x_i - \bar{x}\right)^2}
#' \deqn{SS_y = \sum_{i=1}^n \left(y_i - \bar{y}\right)^2}
#' \deqn{\bar{x} = \frac{\sum_{i=1}^n x_i}{n}}
#' \deqn{\bar{y} = \frac{\sum_{i=1}^n y_i}{n}}
#' 
#' **Symbols used:**
#' \itemize{
#' \item \eqn{n} the number of pairs (sample size)
#' \item \eqn{x_i} the i-th score in the first variable
#' \item \eqn{y_i} the i-th score in the second variable
#' }
#' 
#' The test if test="t" is used is from Pugh and Winslow (1966, pp. 196,199):
#' \deqn{sig. = 2\times\left(1 - \text{T}\left(\left|t_r\right|, df\right)\right)}
#' With:
#' \deqn{t_r = r\times\sqrt{\frac{n - 2}{1 - r^2}}}
#' \deqn{df = n - 2}
#' 
#' The test if test="z" is used is based on a Fisher transformation (Fisher, 1915, p. 521):
#' \deqn{sig. = 2\times\left(1 - \Phi\left(\left|z_r\right|\right)\right)}
#' With:
#' \deqn{z_r = \text{atanh}\left(r\right)\times\sqrt{n - 3}}
#' This is derived since the Fisher transformation has a standard error of:
#' \deqn{SE = \frac{1}{n - 3}}
#' As a source for this standard error Fisher (1921) is sometimes reported, but couldn't clearly find it 
#' in there. It can for example be found in Steiger (1980, p. 246) who refers to Olkin and Siotani (1964).
#' 
#' The correlation coefficient is biased and can be adjusted. There are many different adjustments suggested.
#' For a great overview see Raju et al. (1997).
#' 
#' Fisher (1915, p. 521) - adj="fisher":
#' \deqn{r_{adj} = r\times\left(1 + \frac{1 - r^2}{2 \times n}\right)}
#' 
#' Smith (Ezekiel, 1929, p. 100) - adj="smith":
#' \deqn{r_{adj} = \sqrt{1 - \frac{1 - r^2}{1 - \frac{2}{n}}} = \sqrt{1 - \frac{n}{n - 2}\times\left(1 - r^2\right)}}
#'  
#' Wherry (1931, p. 451) - adj="wherry":
#' \deqn{r_{adj} = \sqrt{\frac{\left(n - 1\right)\times r^2 - 1}{n - 2}} = \sqrt{1 - \left(1 - r^2\right)\times\frac{n - 1}{n - 2}}}
#' 
#' Ezekiel (1930 as cited in Raju et al., 1997, p. 295) - adj="ezekiel":
#' \deqn{r_{adj} = \sqrt{1 - \frac{n - 1}{n - 3}\times\left(1 - r^2\right)}}
#' 
#' Olkin-Pratt (1958, p. 211) - adj="olkin-pratt-1":
#' \deqn{r_{adj} = r\times \text{HG}\left(\frac{1}{2}, \frac{1}{2}, \frac{n-1}{2}, 1 - r^2\right)}
#' 
#' Olkin-Pratt (1958, p. 203) - adj="olkin-pratt-2"
#' \deqn{r_{adj} = r\times\left(1 + \frac{1 - r^2}{2 \times \left(n - 3\right)}\right)}
#' 
#' Cattin (1980a, p. 64; 1980b, p. 409) - adj="cattin":
#' \deqn{r_{adj} = \sqrt{1 - \left(1 - r^2\right)\times\left(1 + \frac{2\times\left(1 - r^2\right)}{n-1} + \frac{8\times\left(1 - r^2\right)^2}{\left(n - 3\right)\times\left(n +1\right)}\right)}}
#' 
#' Pratt (1964, as cited in Claudy, 1978, p. 597) - adj="pratt":
#' \deqn{r_{adj} = \sqrt{1 - \left(1 - r^2\right)\times\left(1 + \frac{2\times\left(1 - r^2\right)}{n-4.3}\right)}}
#' 
#' Herzberg (1969, p. 5) - adj="herzberg":
#' \deqn{r_{adj} = \sqrt{1 - \left(1 - r^2\right)\times\left(1 + \frac{2\times\left(1 - r^2\right)}{n-1}\right)}}
#' 
#' Claudy (1978, p. 603) - adj="claudy":
#' \deqn{r_{adj} = \sqrt{1 - \frac{n-4}{n-2}\times\left(1 - r^2\right)\times\left(1 + \frac{2\times\left(1 - r^2\right)}{n-1}\right)}}
#' 
#' **Alternatives**
#' 
#' R's *stats* library has a similar function.
#' 
#' cor.test(var1, var2)
#' 
#' @references 
#' Cattin, P. (1980a). Note on the estimation of the squared cross-validated multiple correlation of a regression model. *Psychological Bulletin, 87*(1), 63-65. https://doi.org/10.1037/0033-2909.87.1.63
#' 
#' Cattin, P. (1980b). Estimation of the predictive power of a regression model. *Journal of Applied Psychology, 65*(4), 407-414. https://doi.org/10.1037/0021-9010.65.4.407
#' 
#' Claudy, J. G. (1978). Multiple regression and validity estimation in one sample. *Applied Psychological Measurement, 2*(4), 595-607. https://doi.org/10.1177/014662167800200414
#' 
#' Ezekiel, M. (1929). The application of the theory of error to multiple and curvilinear correlation. *Journal of the American Statistical Association, 24*(165), 99-104. https://doi.org/10.2307/2277015
#' 
#' Ezekiel, M. (1941). Methods of correlation analysis (2nd ed.). John Wiley & Sons.Fisher, R. A. (1915). Frequency distribution of the values of the correlation coefficient in samples from an indefinitely large population. *Biometrika, 10*(4), 507-521. https://doi.org/10.2307/2331838
#' 
#' Fisher, R. A. (1921). On the "probable error" of a coefficient of correlation deduced from a small sample. *Metron, 1*, 3-32.
#' 
#' Herzberg, P. A. (1969). The parameters of cross-validation. *Psychometrika Monograph Supplement, 34*(2), 1-70.
#' 
#' Olkin, I., & Pratt, J. W. (1958). Unbiased estimation of certain correlation coefficients. *The Annals of Mathematical Statistics, 29*(1), 201-211. https://doi.org/10.1214/aoms/1177706717
#' 
#' Olkin, I., & Siotani, M. (1964). Asymptotic distribution of functions of a correlation matrix. Laboratory for Quantitative Research in Education, Department of Statistics, Stanford University.
#' 
#' Pearson, K. (1896). Mathematical Contributions to the Theory of Evolution. III. Regression, Heredity, and Panmixia. *Philosophical Transactions of the Royal Society of London. *(A.), 1896, 253-318.
#' 
#' Pugh, E. M., & Winslow, G. H. (1966). *The analysis of physical measurements*. Addison-Wesley.
#' 
#' Raju, N. S., Bilgic, R., Edwards, J. E., & Fleer, P. F. (1997). Methodology review: Estimation of population validity and cross-validity, and the use of equal weights in prediction. *Applied Psychological Measurement, 21*(4), 291-305. https://doi.org/10.1177/01466216970214001
#' 
#' Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. *Psychological Bulletin, 87*(2), 245-251. https://doi.org/10.1037/0033-2909.87.2.245
#' 
#' Wherry, R. J. (1931). A new formula for predicting the shrinkage of the coefficient of multiple correlation. *The Annals of Mathematical Statistics, 2*(4), 440-457. https://doi.org/10.1214/aoms/1177732951
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
r_pearson <- function(field1, 
                      field2, 
                      corr=c("none", "wherry", "fisher", "olkin-pratt-1", "olkin-pratt-2", 
                             "olkin-pratt-3", "smith", "cattin", "pratt", "herzberg"),
                      test=c("t", "z")
                      ){
  
  if (length(test)>1) {
    test="t"
  }
  
  if (length(corr)>1) {
    corr="none"
  }
  
  dFr = na.omit(data.frame(field1, field2))
  
  n = nrow(dFr)
  
  varX = var(dFr$field1)
  varY = var(dFr$field2)
  
  mX = mean(dFr$field1)
  mY = mean(dFr$field2)
  
  sxy = sum((dFr$field1 - mX)*(dFr$field2 - mY))
  
  r = sxy/((n-1)*sqrt(varX*varY))
  
  if (corr=="fisher") {
    r = r*(1+(1 - r**2)/(2*n))
  }
  else if (corr=="smith") {
    r = sqrt(1 - n/(n - 2)*(1 - r**2))
  }
  else if (corr=="wherry") {
    r = sqrt(1 - (n - 1)/(n - 2)*(1 - r**2))
  }
  else if (corr=="ezekiel") {
    r = sqrt(1 - (n - 1)/(n - 3)*(1 - r**2))
  }
  else if (corr=="olkin-pratt-1") {
    r = r*hyperg_2F1(1/2, 1/2, (n - 1)/2, 1 - r**2)
  }
  else if (corr=="olkin-pratt-2") {
    r = r*(1+(1 - r**2)/(2*(n-3)))
  }
  else if (corr=="olkin-pratt-3") {
    r = sqrt(1 - (1 - r**2)*gsl::hyperg_2F1(1, 1, (n - 1)/2, 1 - r**2))
  }
  else if (corr=="cattin") {
    r = sqrt(1 - (1 - r**2)*(1 + 2*(1 - r**2)/(n - 1) + 8*(1 - r**2)**2/((n - 3)*(n + 1))))
  }
  else if (corr=="pratt") {
    r = sqrt(1 - (1 - r**2)*(1 + 2*(1 - r**2)/(n - 4.3)))
  }
  else if (corr=="herzberg") {
    r = sqrt(1 - (1 - r**2)*(1 + 2*(1 - r**2)/(n - 1)))
  }
  else if (corr=="claudy") {
    r = sqrt(1 - (n - 4)*(1 - r**2)/(n - 3)*(1 + 2*(1 - r**2)/(n - 1)))
  }
  
  
  if (test=="t") {
    t = r*sqrt((n - 2)/(1 - r**2))
    df = n - 2
    pValue = 2*(1 - pt(abs(t), df))
    statistic = t    
  }
  else if (test=="z") {
    df = NA  
    z = abs(atanh(r))*sqrt(n - 3)
    pValue = 2*(1 - pnorm(abs(z)))
    statistic = z
  }
  
  results = data.frame(r, statistic, df, pValue)
  colnames(results) = c("r", "statistic", "df", "p-value")  
  return(results)
  
}



