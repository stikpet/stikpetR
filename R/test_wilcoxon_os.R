#' One-Sample Wilcoxon Signed Rank Test
#' 
#' @param data A vector with the data as numbers
#' @param mu optional hypothesized median, otherwise the midrange will be used
#' @param ties optional boolean to use a tie correction (default is True)
#' @param appr optional which method to use for approximation (default is "wilcoxon")
#' @param eqMed optional method to deal with scores equal to hypMed (default is "wilcoxon")
#' @param cc optional boolean to use a continuity correction (default is FALSE)
#' @returns
#' A dataframe with:
#' \item{mu}{the hypothesized median according to the null}
#' \item{W}{the Wilcoxon W value}
#' \item{statistic}{test statistic}
#' \item{df}{degrees of freedom (only applicable for Iman t approximation)}
#' \item{pValue}{significance (p-value)}
#' \item{testUsed}{description of the test used}
#' 
#' @description 
#' The one-sample Wilcoxon signed rank test is often considered the non-parametric version of a one-sample t-test.
#' It can be used to determine if the median is significantly different from an hypothesized value. It actually 
#' doesn't always tests this specifically, but more if the mean rank is significantly different.
#' 
#' If the p-value is the probability of a result as in the sample, or more extreme, if the assumption about the 
#' population would be true. If this is below a certain threshold (usually 0.05) the assumption about the 
#' population is rejected. For this test the assumed median for the population is then incorrect.
#' 
#' Results in software packages for this test can vary, since there are a few different approaches. Especially if 
#' there are so-called ties. See the details for more information.
#' 
#' @details 
#' The unadjusted test statistic is given by:
#' \deqn{W=\sum_{i=1}^{n_{r}^{+}}r_{i}^{+}}
#' With:
#' \deqn{r=\textup{rank}(|d|)}
#' \deqn{d_{i}=y_{i}-\theta}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_{r}^{+}} is the number of ranks with a positive deviation from the hypothesized median
#' \item \eqn{r_{i}^{+}} the i-th rank of the ranks with a positive deviation from the hypothesized median
#' \item \eqn{\theta} is the median tested (the hypothesized median).
#' \item \eqn{y_i} is the i-th score of the variable after removing scores that were equal to \eqn{\theta}
#' }
#' 
#' If there are no ties, an exact method can be used, using the Sign Rank Distribution. 
#' R has this available with *psignrank()*. 
#' The exact test can be found in Zaiontz (n.d.)
#' 
#' **Approximations**
#' 
#' If the sample size is large enough, we can use a normal approximation. 
#' What is large enough varies quite per author. A few examples: n > 8 (slideplayer, 2015), 
#' n > 15 (SigMaxl, n.d.), n > 20 (Wikipedia, n.d.), n > 25 (Harris & Hardin, 2013), 
#' n > 30 (Winthrop, n.d.) . 
#' 
#' The z-statistic is given by (appr="wilcoxon", ties=FALSE, cc=FALSE):
#' \deqn{Z = \frac{W - \mu_w}{\sigma_w}}
#' or with a ties correction (appr="wilcoxon", ties=TRUE, cc=FALSE):
#' \deqn{Z_{adj} = \frac{W - \mu_w}{\sigma_w^*}}
#' With:
#' \deqn{\mu_w = \frac{n_r\times\left(n_r + 1\right)}{4}}
#' \deqn{\sigma_w^2 = \frac{n_r\times\left(n_r + 1\right)\times\left(2\times n_r + 1\right)}{24}}
#' \deqn{\sigma_w^{*2} = \sigma_w^2 - A}
#' \deqn{A = \frac{\sum_{i=1}^k \left(t_i^3 - t_i\right)}{48}}
#' 
#' *Additional symbols used*
#' \itemize{
#' \item \eqn{n_{r}} is the number of ranks used
#' \item \eqn{k} the number of unique ranks
#' \item \eqn{t_i} the frequency of the i-th unique rank
#' }
#' 
#' A Yates continuity correction can simply be applied:
#' In case of no ties (appr="wilcoxon", ties=FALSE, cc=TRUE):
#' \deqn{Z = \frac{\left|W - \mu_w\right| - 0.5}{\sigma_w}}
#' In case of ties (appr="wilcoxon", ties=TRUE, cc=TRUE):
#' \deqn{Z_{adj} = \frac{\left|W - \mu_w\right| - 0.5}{\sigma_w^*}}
#' 
#' 
#' An alternative approximation using the Student t distribution is given by Iman (1974, p. 799).
#' The formula is (appr="imant", ties=FALSE, cc=FALSE):
#' \deqn{t = \frac{W - \mu_w}{\sqrt{\frac{\sigma_w^2\times n_r - \left(W - \mu_w\right)^2}{n_r - 1}}}}
#' or with the ties correction (appr="imant", ties=TRUE, cc=FALSE):
#' \deqn{t = \frac{W - \mu_w}{\sqrt{\frac{\sigma_w^{*2}\times n_r - \left(W - \mu_w\right)^2}{n_r - 1}}}}
#' 
#' The two versions for with a continuity correction are:
#' No ties correction, but continuity (appr="imant", ties=FALSE, cc=TRUE):
#' \deqn{t = \frac{\left|W - \mu_w\right| - 0.5}{\sqrt{\frac{\sigma_w^2\times n_r - \left(\left|W - \mu_w\right| - 0.5\right)^2}{n_r - 1}}}}
#' Both corrections (appr="imant", ties=TRUE, cc=TRUE):
#' \deqn{t = \frac{\left|W - \mu_w\right| - 0.5}{\sqrt{\frac{\sigma_w^{*2}\times n_r - \left(\left|W - \mu_w\right| - 0.5\right)^2}{n_r - 1}}}}
#' 
#' Iman (1974, p. 803) also provides a combination of the t-approximation and the regular z-approximation
#' The equation is given by (appr="imanz"):
#' \deqn{Z_{I} = \frac{Z}{2}\times\left(1 + \sqrt{\frac{n_r - 1}{n_r - Z^2}}\right)}
#' The \eqn{Z} is any of the previous methods.
#' 
#' **Ties with mu**
#' 
#' The default (eqMed="wilcoxon") removes first any scores that are equal to the hypothesized median.
#' There are two alternative methods for this.Both re-define \eqn{d_i} to:
#' \deqn{d_i = x_i - \theta}
#' Where \eqn{x_i} is simply the i-th score.
#' 
#' For the z-split method we only need to re-define:
#' \deqn{W = \frac{\sum_{i=1}^{n_{d_0}}r_{i,0}}{2} + \sum_{i=1}^{n_{r}^{+}}r_{i}^{+}}
#' Where \eqn{n_{d_0}} is the number of scores that equal the hypothesized median, and
#' \eqn{r_{i,0}} is the rank of the i-th score that equals the hypothesized median.
#' 
#' In essence we added half the sum of the ranks that were equal to the hypothesized median.
#' 
#' For the z-split method all other calculations than go the same.
#' 
#' For the Pratt (1959) method we also re-define:
#' \deqn{\mu_w = \frac{n_r\times\left(n_r + 1\right) - n_{d_0}\times\left(n_{d_0} + 1\right)}{4}}
#' \deqn{\sigma_w^2 = \frac{n_r\times\left(n_r + 1\right)\times\left(2\times n_r + 1\right) - n_{d_0}\times\left(n_{d_0} + 1\right)\times\left(2\times n_{d_0} + 1\right)}{24}}
#' 
#' For the Pratt method, the ties correction still excludes the ties 
#' for the scores that equal the hypothesized median, 
#' but for the z-split method it will include them.
#' 
#' For both methods now \eqn{n_r=n}, where n is the number of scores.
#' 
#' The Pratt (1959) method and z-split method were found in Python’s documentation 
#' for scipy’s Wilcoxon function (scipy, n.d.). They also refer to Cureton (1967) for the Pratt method.
#' 
#' @examples 
#' data <- c(1, 2, 5, 1, 1, 5, 3, 1, 5, 1, 1, 5, 1, 1, 3, 3, 3, 4, 2, 4)
#' ts_wilcoxon_os(data)
#' ts_wilcoxon_os(data, ties=FALSE, appr="imanz", eqMed = "zsplit", cc = FALSE)
#' 
#' @section Alternatives:
#' 
#' R's stats library has a similar function: *wilcox.test()*
#' 
#' The *exactRankTests* library has a function for the exact test: *wilcox.exact()*
#' 
#' @seealso 
#' Alternative tests:
#' \itemize{
#' \item one-sample sign test, see \code{\link{ts_sign_os}}
#' \item one-sample trinomial test, see \code{\link{ts_trinomial_os}}
#' } 
#' 
#' As for an effect size measure, see:
#' \itemize{
#' \item Rosenthal Correlation, see \code{\link{r_rosenthal}}
#' \item Rank Biserial Correlation, see \code{\link{r_rank_biserial_os}}
#' \item Dominance, see \code{\link{es_dominance}}
#' }
#' 
#' @references 
#' Cureton, E. E. (1967). The normal approximation to the signed-rank sampling distribution when zero differences are present. *Journal of the American Statistical Association, 62*(319), 1068–1069. https://doi.org/10.1080/01621459.1967.10500917
#' 
#' Harris, T., & Hardin, J. W. (2013). Exact Wilcoxon Signed-Rank and Wilcoxon Mann–Whitney Ranksum Tests. *The Stata Journal, 13*(2), 337–343. https://doi.org/10.1177/1536867X1301300208
#' 
#' Iman, R. L. (1974). Use of a t-statistic as an approximation to the exact distribution of the wildcoxon signed ranks test statistic. *Communications in Statistics, 3*(8), 795–806. https://doi.org/10.1080/03610927408827178
#' 
#' Pratt, J. W. (1959). Remarks on zeros and ties in the Wilcoxon signed rank procedures. *Journal of the American Statistical Association, 54*(287), 655–667. https://doi.org/10.1080/01621459.1959.10501526
#' 
#' scipy. (n.d.). Scipy.stats.wilcoxon. Scipy. https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
#' 
#' SigMaxl. (n.d.). One Sample Wilcoxon Sign Test Exact. Retrieved August 30, 2020, from https://www.sigmaxl.com/OneSampleSignWilcoxonExact.shtml
#' 
#' slideplayer. (2015, June 13). Using statistics to make inferences 6.
#' 
#' Wikipedia. (n.d.). Wilcoxon signed-rank test. In Wikipedia. Retrieved August 30, 2020, from https://en.wikipedia.org/w/index.php?title=Wilcoxon_signed-rank_test&oldid=974561084
#' 
#' Wilcoxon, F. (1945). Individual comparisons by ranking methods. *Biometrics Bulletin, 1*(6), 80. https://doi.org/10.2307/3001968
#' 
#' Winthrop. (n.d.). The Wilcoxon signed rank test for one sample. Winthrop Univerisy Hospital. https://nyuwinthrop.org/wp-content/uploads/2019/08/wilcoxon-sign-rank-test-one-sample.pdf
#' 
#' Zaiontz, C. (n.d.). Wilcoxon signed ranks exact test. Real Statistics Using Excel. Retrieved January 25, 2023, from https://real-statistics.com/non-parametric-tests/wilcoxon-signed-ranks-test/wilcoxon-signed-ranks-exact-test/
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
ts_wilcoxon_os <- function(data, 
                           mu = NULL, ties = TRUE, 
                           appr = c("wilcoxon", "none", "imanz", "imant"), 
                           eqMed = c("wilcoxon", "zsplit", "pratt"), 
                           cc = FALSE){
  
  #set defaults if not provided
  if (length(appr)>1){appr="wilcoxon"}
  if (length(eqMed)>1){appr="wilcoxon"}
  
  #set hypothesized median to mid range if not provided
  if (is.null(mu)) {
    mu = (min(data) + max(data)) / 2
  }
  
  #sample size (n)
  n = length(data)
  
  #adjust sample size if wilcoxon method is used for equal distance
  if (eqMed == "wilcoxon") {
    nr = n - sum(data == mu, na.rm=TRUE)
  } else{
    nr = n
  }
  
  absDiffs = 0
  #remove scores equal to hypMed if eqMed is wilcoxon
  if (eqMed == "wilcoxon" || appr == "none"){
    data = data[data != mu]
  }
  
  #determine the absolute deviations from the hypMed
  diffs = data - mu
  absDiffs = abs(diffs)
  ranks = rank(absDiffs)
  
  W = sum(ranks[diffs > 0])
  
  df = "n.a."
  if (appr=="none"){
    #check if ties exist
    ranks = ranks[absDiffs != 0]
    
    if (max(table(ranks)) > 1 ){
      return ("ties exist, cannot compute exact method") }
    else {
      Wmin = sum(ranks[diffs < 0])
      pVal = psignrank(min(W, Wmin),length(ranks))*2
      statistic = min(W, Wmin)
      testUsed = "one-sample Wilcoxon signed rank exact test"
    }
  }
  else{
    #add half the equal to median ranks if zsplit is used
    nD0 = sum(diffs==0)
    if (eqMed == "zsplit") {
      W = W + sum(ranks[diffs == 0])/2}
    
    
    rAvg = nr*(nr + 1)/4
    s2 = nr * (nr+1) * (2*nr + 1)/24
    #adjust if Pratt method is used
    if (eqMed == "pratt") {
      #normal approximation adjustment based on Cureton (1967)
      s2 = s2 - nD0 * (nD0 + 1) * (2 * nD0 + 1) / 24
      rAvg = (nr * (nr + 1) - nD0 * (nD0 + 1)) / 4
    }
    
    #the ties correction
    t = 0
    if (ties) {
      #remove ranks of scores equal to hypMed for Pratt (wilcoxon already removed)
      if (eqMed == "pratt") {
        ranks = ranks[absDiffs != 0]
      }
      
      for (i in table(ranks)){
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
      if (cc){
        var = (s2 * nr - (abs(W - rAvg) - 0.5) ^ 2) / (nr - 1)
      }
      else{
        var = (s2 * nr - (W - rAvg) ^ 2) / (nr - 1)
      }
      statistic = num / sqrt(var)
      df = nr - 1
      pVal = 2 * (1 - pt(abs(statistic), df))}
    
    else{
      statistic = num / se
      if (appr == "imanz"){
        statistic = statistic / 2 * (1 + sqrt((nr - 1) / (nr - statistic ^ 2)))}
      pVal = 2 * (1 - pnorm(abs(statistic)))  
      
    }
    
    testUsed = "one-sample Wilcoxon signed rank test"
    
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
    
    if (eqMed == "pratt"){
      testUsed = paste0(testUsed, ", Pratt method for equal to hyp. med. (inc. Cureton adjustment for normal approximation)")}
    else if (eqMed == "zsplit"){
      testUsed = paste0(testUsed, ", z-split method for equal to hyp. med.")}
  }  
  
  pValue = pVal
  testResults <- data.frame(mu, W, statistic, df, pValue, testUsed)
  
  return (testResults)
}