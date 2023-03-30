#' One-Way ANOVA table
#' 
#' @description 
#' Function to generate a one-way anova table
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' scores = c(20, 50, 80, 15, 40, 85, 30, 45, 70, 60, NA, 90, 25, 40, 70, 65, NA, 70, 98, 40, 65, 60, 35, NA, 50, 40, 75, NA, 65, 70, NA, 20, 80, 35, NA, 68, 70, 60, 70, NA, 80, 98, 10, 40, 63, 75, 80, 40, 90, 100, 33, 36, 65, 78, 50)
#' groups = c("Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Rotterdam", "Haarlem", "Diemen", "Haarlem", "Diemen", "Haarlem", "Haarlem", "Haarlem", "Haarlem", "Haarlem")
#' he_owa_table(scores, groups)
#' 
#' @export
he_owa_table <- function(scores, groups){
  exDf = na.omit(data.frame(scores, groups))
  counts <- setNames(aggregate(exDf$scores~exDf$groups, FUN=length), c("category", "n"))
  means <- setNames(aggregate(exDf$scores~exDf$groups, FUN=mean), c("category", "mean"))
  vars <- setNames(aggregate(exDf$scores~exDf$groups, FUN=var), c("category", "var"))
  myRes <- merge(counts, means, by = 'category')
  myRes <- merge(myRes, vars, by = 'category')
  
  n = sum(myRes$n)
  k = length(myRes$n)
  
  SSb = sum(vars$var*(counts$n-1))
  SSt = var(exDf$scores)*(n-1)
  SSw = SSt - SSb
  SS = c(SSb, SSw, SSt)
  df = c(k - 1, n - k, NA)
  MS = c(SSb/df[1], SSw/df[2], NA)
  F = c(MS[1]/MS[2], NA, NA)
  
  results = data.frame(SS, df, MS, F, row.names = c("between", "within", "total"))  
  
  return(results)
}

