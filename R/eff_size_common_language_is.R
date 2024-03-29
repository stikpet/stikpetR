#' Common Language (CL/CLES) (Independent Samples)
#' @description 
#' the probability that a randomly selected score from the one population will be greater than a randomly sampled score from the other population.
#' 
#' @param catField A vector with the categorical data
#' @param scaleField A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' @param dmu Optional difference according to null hypothesis (default is 0)
#' @param method Optional method to use, either approximate ("appr") or use brute-force ("brute".
#' 
#' @returns 
#' A dataframe with:
#' \item{CLE cat. 1}{the effect size for the first category}
#' \item{CLE cat. 2}{the effect size for  the second category}
#' 
#' @details
#' 
#' The formula used is (McGraw & Wong, 1992, p. 361):
#' \deqn{CL = \Phi\left(z\right)}
#' 
#' With:
#' \deqn{z = \frac{\left|\bar{x}_1 - \bar{x}_2\right|}{\sqrt{s_1^2 + s_2^2}}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' \item \eqn{\Phi\left(\dots\right)} the cumulative density function of the standard normal distribution
#' }
#' 
#' The CLE for the other category is simply 1 - CLE.
#' 
#' If the brute-force method is used, the chance of a value of one category will be higher than a random score of the second category is calculated by counting all possible pairs. In case of ties a half is added.
#' 
#' @references 
#' McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. *Psychological Bulletin, 111*(2), 361–365. https://doi.org/10.1037/0033-2909.111.2.361
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['age']
#' ex1 = replace(ex1, ex1=="89 OR OLDER", "90")
#' es_common_language_is(df1['sex'], ex1)
#' 
#' #Example 2: vectors
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("nat","int","int","nat","int","int","nat","nat","int","int","int",
#' "int","int","int","nat","int",NA,"nat","int","int")
#' es_common_language_is(groups, scores)
#' 
#' 
#' @export
es_common_language_is <- function(catField, scaleField, categories=NULL, dmu=0, method="appr"){
  
  #remove rows with missing values
  df = data.frame(scaleField, catField)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  df$score = as.numeric(df$score)
  
  #the two categories
  if (!is.null(categories)){
    cat1 = categories[1]
    cat2 = categories[2]
  }
  else {
    cat1 = names(table(df$group))[1]
    cat2 = names(table(df$group))[2]
  }
  
  X1 = df$score[df$group == cat1]
  X2 = df$score[df$group == cat2]
  
  n1 = length(X1)
  n2 = length(X2)
  
  var1 = var(X1)
  var2 = var(X2)
  
  m1 = mean(X1)
  m2 = mean(X2)

  if (method=="appr"){
    z = abs(m1 - m2 - dmu)/sqrt(var1 + var2)
    c1 = pnorm(z)
    c2 = 1 - c1        
  }    
  else if (method=="brute"){
    c1 = 0
    for (i in X1){
      p = 0
      for (j in X2){
        if (i>j){
          p = p + 1}
        if (i==j){
          p = p + 0.5}
      }
      c1 = c1 + p/n2
    }
    c1 = c1 / n1
    
    c2 = 0
    for (i in X2){
      p = 0
      for (j in X1){
        if (i>j){
          p = p + 1}
        if (i==j){
          p = p + 0.5}
      }
      c2 = c2 + p/n1
    }
    c2 = c2 / n2
  }
  
  if ((m1>m2 && c1<c2) || (m2>m1 && c2 < c1)){
    c3 = c1
    c1 = c2
    c2 = c3
  }
  
  
  
  results <- data.frame(c1, c2)
  colnames(results) = c(paste("CLE", cat1), paste("CLE", cat2))
  
  return(results)
}