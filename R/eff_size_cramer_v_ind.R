#' Cramer's V for Independence Test
#' 
#' @param chi2 the chi-square test statistic
#' @param n the sample size
#' @param r the number of categories in the first variable (i.e. the number of rows)
#' @param c the number of categories in the second variable (i.e. the number of columns)
#' @param cc c(NULL, "bergsma") optional to indicate correction to use (default is NULL)
#' @return Cramer's V value
#' 
#' 
#' @details
#' The formula used is:
#' \deqn{V=\sqrt\frac{\chi^{2}}{n\times \left(\text{min}\left(r, c\right) - 1\right)}}
#' 
#' *Symbols used*:
#' \itemize{
#' \item \eqn{r} the number of categories in the first variable (i.e. the number of rows)
#' \item \eqn{c} the number of categories in the second variable (i.e. the number of columns)
#' \item \eqn{n} the sample size, i.e. the sum of all frequencies
#' \item \eqn{\chi^{2}} the chi-square statistic
#' }
#' 
#' The Bergsma correction uses a different formula (Bergsma, 2013, pp. 324-325):
#' \deqn{V_B = \sqrt{\frac{\tilde{\varphi}^2}{\text{min}\left(\tilde{r}, \tilde{c}\right) - 1}}}
#' With:
#' \deqn{\tilde{\varphi}^2 = \text{max}\left(0,\varphi^2 - \frac{\left(r - 1\right)\times\left(c - 1\right)}{n - 1}\right)}
#' \deqn{\tilde{r} = r - \frac{\left(r - 1\right)^2}{n - 1}}
#' \deqn{\tilde{c} = r - \frac{\left(c - 1\right)^2}{n - 1}}
#' \deqn{\varphi^2 = \frac{\chi^{2}}{n}}
#' 
#' @references 
#' Bergsma, W. (2013). A bias-correction for Cramér’s and Tschuprow’s. Journal of the Korean Statistical Society, 42(3), 323–328. https://doi.org/10.1016/j.jkss.2012.10.002
#' 
#' Cramér, H. (1946). Mathematical methods of statistics. Princeton University Press.
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#'  
#' @examples 
#' chi2Val = 16.98975
#' n = 1941
#' nRows = 5
#' nCols = 2
#' es_cramer_v_ind(chi2Val, n, nRows, nCols)
#' es_cramer_v_ind(chi2Val, n, nRows, nCols, cc="bergsma")
#' 
#' @export
es_cramer_v_ind <- function(chi2, n, r, c, cc=NULL){
  m = min(r, c)
  
  if (!is.null(cc) && cc=="bergsma") {
    phi2 = chi2/n
    mHat = m - (m - 1)^2/(n - 1)
    df = (r - 1)*(c - 1)
    phi2 = max(0, phi2 - df/(n - 1))
    
    es = sqrt(phi2/(mHat - 1))
  }
  else{
    es = sqrt(chi2/(n*min(r-1, c-1)))  
  }
  
  return(es)
  
}