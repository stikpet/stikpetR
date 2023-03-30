#' Eta Squared (Maximum Corrected) for Cochran Q
#' 
#' @param Q the Cochran Q statistic
#' @param n the sample size (number of rows)
#' @param k the number of variables (number of columns)
#' @returns 
#' \item{es}{the effect size measure}
#' 
#' @details 
#' The formula used (Serlin et al., 1982 p. 788):
#' \deqn{\eta_Q^2 = \frac{Q}{n\times\left(k - 1\right)}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{n} the number of rows
#' \item \eqn{k} the number of columns
#' \item \eqn{Q} the Cochran Q statistic
#' }
#' 
#' The Cochran Q statistic can be obtained using *ts_cochran_q()* function.
#' The number of rows and columns of a dataframe with R's *nrow(dataframe)* and
#' *ncol(dataframe)* functions.
#' 
#' @examples 
#' Q = 0.1205357
#' n = 16
#' k = 4
#' es_eta_sq_mc(Q, n, k)
#' 
#' @references 
#' Serlin, R. C., Carr, J., & Marascuilo, L. A. (1982). A measure of association for selected nonparametric procedures. *Psychological Bulletin, 92*(3), 786–790. https://doi.org/10.1037/0033-2909.92.3.786
#' 
#' @export
es_eta_sq_mc <- function(Q, n, k){
  
  es = Q/(n*(k - 1))
  
  return(es)
}