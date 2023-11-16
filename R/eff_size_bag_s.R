#' Bennett-Alpert-Goldstein S 
#' @description
#' An effect size meaure, that measures the how strongly two raters or variables, agree with each other. 
#' 
#' It takes the proportions of cases that both agree, and adjusts for the number of categories. Scott's pi (see es_scott_pi()) does this as well, and improves on this measure.
#' 
#' @param field1 vector, the first categorical field
#' @param field2 vector, the first categorical field
#' @param categories vector, optional, order and/or selection for categories of field1 and field2
#' 
#' @returns
#' S, the Bennett-Alpert-Goldstein value
#' 
#' @details 
#' The formula used (Bennett et al., 1954, p. 307):
#' \deqn{S = \frac{k}{k-1}\times\left(p_0 - \frac{1}{k}\right)}
#' 
#' With:
#' \deqn{P = \sum_{i=1}^r F_{i,i}}
#' \deqn{p_0 = \frac{P}{n}}
#' 
#' *Symbols used*
#' \itemize{
#' \item \eqn{F_{i,j}}, the observed count in row i and column j.
#' \item \eqn{r}, is the number of rows (categories in the first variable)
#' \item \eqn{n}, is the total number of scores
#' }
#' 
#' @references
#' Bennett, E. M., Alpert, R., & Goldstein, A. C. (1954). Communications through limited response questioning. *Public Opinion Quarterly, 18*(3), 303. doi:10.1086/266520
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_bag_s <- function(field1, field2, categories=NULL){
  #create the cross table
  ct = tab_cross(field1, field2, categories, categories, totals="include")
  
  #basic counts
  k = nrow(ct)-1
  n = ct[k+1, k+1]
  
  p0 = 0
  for (i in 1:k){p0 = p0 + ct[i, i]}
  p0 = p0/n
  
  S = k / (k - 1) * (p0 - 1 / k)
  
  return (S)
}