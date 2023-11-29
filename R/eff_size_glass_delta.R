#' Glass Delta
#' @description
#' An effect size measure when comparing two means, with a specified control group.
#' 
#' @param catField A vector with the categorical data
#' @param scaleField A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' @param dmu Optional difference according to null hypothesis (default is 0)
#' @param control Optional to indicate which category to use as control group. Default is first category found.
#' 
#' @return Glass Delata value
#' 
#' @details
#' The formula used is (Glass, 1976, p. 7):
#' \deqn{\delta = \frac{\bar{x}_1 - \bar{x}_2}{s_2}}
#' 
#' With:
#' \deqn{s_2 = \sqrt{\frac{\sum_{i=1}^{n_2} \left(x_{2,i} - \bar{x}_2\right)^2}{n_2 - 1}}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' }
#' 
#' Glass actually uses a ‘control group’ and \eqn{s_2} is then the standard deviation of the control group. 
#' 
#' @references 
#' Glass, G. V. (1976). Primary, secondary, and meta-analysis of research. *Educational Researcher, 5*(10), 3–8. https://doi.org/10.3102/0013189X005010003
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
#' es_glass_delta(df1['sex'], ex1)
#' 
#' #Example 2: vectors
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("nat","int","int","nat","int","int","nat","nat","int","int",
#' "int","int","int","int","nat","int",NA,"nat","int","int")
#' es_glass_delta(groups, scores)
#' 
#' @export
es_glass_delta <- function(catField, 
                           scaleField, 
                           categories=NULL, 
                           dmu=0, 
                           control=NULL){
  
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
  n = n1 + n2
  
  m1 = mean(X1)
  m2 = mean(X2)
  var1 = var(X1)
  var2 = var(X2)
  
  sd1 = sqrt(var1)
  sd2 = sqrt(var2)
  
  if (is.null(control) || control==cat2){
    s = sd2}
  else{
    if (control==cat1){
      s= sd1
    }
  }
  gd = (m1- m2)/s
  
  return(gd)
  
}