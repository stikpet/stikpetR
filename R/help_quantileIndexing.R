#' Quartile Indexing
#' 
#' Helper function for \code{\link{me_quantiles}} and \code{\link{he_quantileIndexing}} to return the index number of the quantiles.
#' 
#' @param data dataframe with scores as numbers
#' @param k : number of quantiles
#' @param method indexing method to use
#' @returns
#' a vector with the quantiles
#' 
#' @details 
#' Six alternatives for the indexing is:
#' Most basic (**SAS1**):
#' \deqn{iQ_i = n\times p_i}
#' 
#' **SAS4** method uses for indexing (SAS, 1990, p. 626; Snedecor, 1940, p. 43):
#' \deqn{iQ_i = \left(n + 1\right)\times p_i}
#' 
#' **Hog and Ledolter** use for their indexing (Hogg & Ledolter, 1992, p. 21; Hazen, 1914, p. ?):
#' 
#' \deqn{iQ_i = n\times p_i + \frac{1}{2}}
#' 
#' **MS Excel** uses for indexing (Gumbel, 1939, p. ?; Hyndman & Fan, 1996, p. 363):
#' \deqn{iQ_i = \left(n - 1\right)\times p_i + 1}
#' 
#' **Hyndman and Fan** use for their 8th version (Hyndman & Fan, 1996, p. 363):
#' 
#' \deqn{iQ_i = \left(n + \frac{1}{3}\right)\times p_i + \frac{1}{3}}
#' 
#' **Hyndman and Fan** use for their 9th version (Hyndman & Fan, 1996, p. 364):
#' \deqn{iQ_i = \left(n + \frac{1}{4}\right)\times p_i + \frac{3}{8}}
#' 
#' @references 
#' Gumbel, E. J. (1939). La Probabilité des Hypothèses. Compes Rendus de l’ Académie des Sciences, 209, 645–647.
#' 
#' Hazen, A. (1914). Storage to be provided in impounding municipal water supply. Transactions of the American Society of Civil Engineers, 77(1), 1539–1640. https://doi.org/10.1061/taceat.0002563
#' 
#' Hogg, R. V., & Ledolter, J. (1992). Applied statistics for engineers and physical scientists (2nd int.). Macmillan.
#' 
#' Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361–365. https://doi.org/10.2307/2684934
#' 
#' SAS. (1990). SAS procedures guide: Version 6 (3rd ed.). SAS Institute.
#' 
#' Snedecor, G. W. (1940). Statistical methods applied to experiments in agriculture and biology (3rd ed.). The Iowa State College Press.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
he_quantileIndexing <- function(data, 
                                k=4, 
                                method=c("sas1", "sas4", "hl", "excel", "hf8", "hf9")){
  
  if (length(method)>1){method = "sas1"}
  data = sort(data)
  n = length(data)
  props = 1/k  
  indexes = c()
  for (i in 0:k){
    if (method=="sas1"){q = n*(i*props)}
    else if (method=="sas4"){q = (n + 1)*(i*props)}
    else if (method=="hl"){q = n*(i*props) + 1/2}
    else if (method=="excel"){q = (n - 1)*(i*props) + 1}
    else if (method=="hf8"){q = (n + 1/3)*(i*props) + 1/3}
    else if (method=="hf9"){q = (n + 1/4)*(i*props)}
    
    if (q < 1) {q = 1}
    else if (q > n){q = n}
    
    indexes = append(indexes, q)
  }
  
  results <- indexes
  
  return (results)
}