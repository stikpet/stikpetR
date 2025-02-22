#' One-Sample Trinomial Test
#' 
#' @description
#' A test that could be used with ordinal data that includes ties
#' 
#' Similar as a sign-test but instead of ignoring scores that are tied with the hypothesized median they get included, hence instead of the binomial distribution, this will use the trinomial distribution.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/OiVUfX5lRww) and the test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tests/WilcoxonSignedRankOneSample.html).
#' 
#' @param data A vector or dataframe
#' @param levels optional list to indicate what values represent
#' @param mu optional hypothesized median, otherwise the midrange will be used
#' 
#' @returns 
#' A dataframe with:
#' \item{mu}{he hypothesized median}
#' \item{n-pos}{the number scores above mu}
#' \item{n-neg}{the number scores below mu}
#' \item{n-tied}{the number of scores tied with mu}
#' \item{p-value}{significance (p-value)}
#' \item{test}{description of the test used}
#' 
#' @details 
#' The p-value is calculated using (Bian et al., 2009, p. 6):
#' \deqn{p = 2\times \sum_{i=n_d}^n \sum_{j=0}^{\lfloor \frac{n - i}{2} \rfloor} \text{tri}\left(\left(j, j+i, n - i\right), \left(p_{pos}, p_{neg}, p_0\right) \right)}
#' 
#' With:
#' \deqn{p_0 = \frac{n_0}{n}}
#' \deqn{p_{pos} = p_{neg} = \frac{1 - p_0}{n}}
#' \deqn{\left|n_{pos} - n_{neg}\right|}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{n_0} the number of scores equal to the hypothesized median
#' \item \eqn{n_{pos}} the number of scores above the hypothesized median
#' \item \eqn{n_{neg}} the number of scores below the hypothesized median
#' \item \eqn{p_0} the probability of the a score in the sample being equal to the hypothesized median
#' \item \eqn{p_{pos}} the population proportion of a score being above the hypothesized median
#' \item \eqn{p_{neg}} the population proportion of a score being below the hypothesized median
#' \item \eqn{\text{tri}\left(…,… \right)} the trinomial probability mass function
#' }
#' 
#' The paired version of the test is described in Bian et al. (1941), while Zaiontz (n.d.) mentions it can also be used for one-sample situations.
#' 
#' @section Before, After and Alternatives:
#' Before this measure you might want an impression using a frequency table or a visualisation:
#' \code{\link{tab_frequency}}, for a frequency table
#' \code{\link{vi_bar_stacked_single}}, or Single Stacked Bar-Chart.
#' \code{\link{vi_bar_dual_axis}}, for Dual-Axis Bar Chart.
#' 
#' After this you might want to determine an effect size measure:
#' \code{\link{es_common_language_os}}, for the Common Language Effect Size. 
#' \code{\link{es_dominance}}, for the Dominance score.
#' \code{\link{r_rank_biserial_os}}, for the Rank-Biserial Correlation
#' 
#' Alternative tests:
#' \code{\link{ts_sign_os}}, for One-Sample Sign Test.
#' \code{\link{ts_wilcoxon_os}}, for One-Sample Wilcoxon Signed Rank Test.
#' 
#' @references 
#' Bian, G., McAleer, M., & Wong, W.-K. (2009). A trinomial test for paired data when there are many ties. *SSRN Electronic Journal*. https://doi.org/10.2139/ssrn.1410589
#' 
#' Zaiontz, C. (n.d.). Trinomial test. Real Statistics Using Excel. Retrieved March 2, 2023, from https://real-statistics.com/non-parametric-tests/trinomial-test/
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' #Example 1: Dataframe
#' ex1 = df2[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' ts_trinomial_os(ex1, levels=order)
#' 
#' #Example 2: Numeric data
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' ts_trinomial_os(ex2)
#' 
#' @export
ts_trinomial_os <- function(data, levels=NULL, mu=NULL){
  testUsed = "one-sample trinomial test"
  data = na.omit(data)
  
  if (!is.null(levels)){
    dataN = factor(data, ordered = TRUE, levels = levels)
    dataN = as.numeric(dataN)
  }
  else{dataN = data}
  
  #set hypothesized median to mid range if not provided
  if (is.null(mu)) {
    mu = (min(dataN) + max(dataN)) / 2
  }

  pos = sum(dataN>mu)
  neg = sum(dataN<mu)
  ties = sum(dataN==mu)
  
  n = pos + neg + ties
  
  p0 = ties/n
  p1 = (1 - p0)/2
  k = abs(pos - neg)
  
  sig=0
  for (z in k:n) {
    for (i in 0:floor((n - z)/2)) {
      sig = sig + dmultinom(c(i, i+z, n - i - (i+z)), prob = c(p1, p1, p0))
    }
  }
  
  pValue = sig*2  
  
  results = data.frame(mu, pos, neg, ties, pValue, testUsed)
  colnames(results)<-c("mu", "n-pos.", "n-neg.", "n-tied.", "p-value", "test")
  
  return(results)
  
}
