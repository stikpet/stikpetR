#' Fisher-Freeman-Halton Exact test
#' 
#' @param var1 A vector with the data from the first variable
#' @param var2 A vector with the data from the second variable
#' @return the two-tailed p-value/sig.
#' 
#' @examples 
#' var1 <- c("female", "female","female","female","female","female","female","female",
#'           "female","female","female", "male", "male", "male", "male", "male", "male", 
#'           "male", "male", "male", "male", "male", "male", "male", "male", "male", "male",
#'           "male", "male", "male", "male", "male", "male", "male", "male", "male", "male", 
#'           "male", "male", "male", "male", "male")
#' var2 <- c("nl", "nl","nl","nl","nl","nl","nl","nl", "other", "other", "other",
#'           "nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl","nl",
#'           "other", "other", "other", "other", "other", "other", "other", "other", "other", 
#'           "other", "other", "other", "other", "other", "other")
#' ts_fisher(var1, var2)
#' 
#' @details 
#' This simply uses R's *fisher.test()* function from the stats library.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Fisher, R. A. (1922). On the Interpretation of χ2 from Contingency Tables, and the Calculation of P. *Journal of the Royal Statistical Society, 85*(1), 87–94. https://doi.org/10.2307/2340521
#' 
#' Freeman, G. H., & Halton, J. H. (1951). Note on an exact treatment of contingency, goodness of fit and other problems of significance. *Biometrika, 38*(1/2), 141–149. https://doi.org/10.2307/2332323
#'  
#' @export
ts_fisher_freeman_halton <- function(var1, var2){
  
  data = data.frame(var1, var2)
  
  #remove missing values
  data = na.omit(data)
  
  #Create a cross table first
  dataTable = table(data)
  
  res = fisher.test(dataTable)
  
  return(res$p.value)
}

