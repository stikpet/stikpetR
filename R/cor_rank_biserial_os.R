#' Rank biserial correlation coefficient (one-sample)
#' 
#' This function will calculate Rank biserial correlation coefficient (one-sample)
#' 
#' @param data vector with the numeric scores
#' @param levels optional vector with levels in order
#' @param mu optional parameter to set the hypothesized median. If not used the midrange is used
#' @return dataframe with the hypothesized median (mu) and the effect size measure
#' 
#' @details 
#' The formula used (Kerby, 2014, p. 5):
#' \deqn{r_{rb} = \frac{\left|R_{pos} - R_{neg}\right|}{R}}
#' This is actually the same as (King & Minium, 2008, p. 403):
#' \deqn{r_{rb} = \frac{4\times\left|R_{min} - \frac{R_{pos} + R_{min}}{2}\right|}{n\times\left(n + 1\right)}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{R_{pos}} the sum of the ranks with a positive deviation from the hypothesized median 
#' \item \eqn{R_{neg}} the sum of the ranks with a positive deviation from the hypothesized median 
#' \item \eqn{R_{min}} the minimum of \eqn{R_{pos}},\eqn{R_{neg}}
#' \item \eqn{n} the number of ranks with a non-zero difference with the hypothesized median 
#' \item \eqn{R} the sum of all ranks, i.e. \eqn{R_{pos} + R_{neg}}
#' }
#' 
#' If no hypothesized median is provided, the midrange is used, defined as:
#' \deqn{\frac{x_{max} - x_{min}}{2}}
#' Where \eqn{x_{max}} is the maximum value of the scores, and \eqn{x_{min}} the minimum
#' 
#' **Alternative**
#' 
#' The *effectsize* library has a similar function: *rank_biserial()*
#' 
#' @references 
#' Kerby, D. S. (2014). The simple difference formula: An approach to teaching nonparametric correlation. *Comprehensive Psychology*, 3, 1–9. https://doi.org/10.2466/11.IT.3.1
#' 
#' King, B. M., & Minium, E. W. (2008). *Statistical reasoning in the behavioral sciences* (5th ed.). John Wiley & Sons, Inc.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' ex1 = df2[['Teach_Motivate']]
#' order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
#' r_rank_biserial_os(ex1, levels=order)
#' 
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' r_rank_biserial_os(ex2)
#' 
#' @export
r_rank_biserial_os <- function(data, levels=NULL, mu=NULL){
  
  if (is.null(levels)){
    dataN = data}
  else{
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    dataN = as.numeric(myFieldOrd)
  }
  
  #set hypothesized median to mid range if not provided
  if (is.null(mu)) {
    mu = (min(dataN) + max(dataN)) / 2
  }
  
  #remove scores equal to hypothesized median
  dataN = dataN[dataN != mu]
  
  # Determine the differences with the hypothesized median, and the corresponding signs
  diffs = abs(dataN - mu)
  signs = sign(dataN - mu)
  
  # Use the rank function to determine the ranks
  ranks = rank(diffs)
  
  # Sum the positive and negative ranks
  rPlus = sum(ranks[signs==1])
  rNeg = sum(ranks[signs==-1])
  
  # Calculate the rank biserial correlation coefficient
  rb = (rPlus - rNeg) / (rPlus + rNeg)
  
  testResults <- data.frame(mu, rb)
  
  return(testResults)
}


