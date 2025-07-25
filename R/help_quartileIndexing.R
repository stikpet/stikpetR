#' Quartile Indexing
#' 
#' Helper function for \code{\link{me_quartiles}} and \code{\link{he_quartileIndexing}} to return the index number
#'  of the first and third quartile for different methods of determining this index.
#' 
#' @param data dataframe with scores as numbers
#' @param method indexing method to use
#' @returns
#' A dataframe with:
#' \item{q1Index}{the index of the first (lower) quartile}
#' \item{q3Index}{the index of the third (upper/higher) quartile}
#' 
#' @details 
#' The **inclusive** method divides the data into two, and then includes the median in each half (if the sample size is odd). 
#' The first and third quarter are then the median of each of these two halves (Tukey, 1977, p. 32).
#' 
#' For the **inclusive** method, the index of the first quartile can be found using:
#' \deqn{iQ_1 = \begin{cases} \frac{n+2}{4} & \text{ if } n \text{ mod } 2 = 0 \\ \frac{n+3}{4} & \text{ else } \end{cases}}
#' 
#' And the third quartile:
#' \deqn{iQ_3 = \begin{cases} \frac{3\times n +2}{4} & \text{ if } n \text{ mod } 2 = 0 \\ \frac{3\times n+1}{4} & \text{ else } \end{cases}}
#' 
#' The **exclusive** method does the same as the inclusive method, but excludes the median in each half (if the sample size 
#' is odd) (Moore & McCabe, 1989, p. 33; Joarder & Firozzaman, 2001, p. 88).
#' 
#' For the **exclusive** method, the index of the first quartile can be found using:
#' \deqn{iQ_1 = \begin{cases} \frac{n+2}{4} & \text{ if } n \text{ mod } 2 = 0 \\ \frac{n + 1}{4} & \text{ else } \end{cases}}
#' 
#' And the third quartile:
#' \deqn{iQ_3 = \begin{cases} \frac{3\times n +2}{4} & \text{ if } n \text{ mod } 2 = 0 \\ \frac{3\times n + 3}{4} & \text{ else } \end{cases}}
#' 
#' Other methods use a different indexing. Six alternatives for the indexing is:
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
#' Gumbel, E. J. (1939). La Probabilite des Hypotheses. Compes Rendus de l' Academie des Sciences, 209, 645-647.
#' 
#' Hazen, A. (1914). Storage to be provided in impounding municipal water supply. Transactions of the American Society of Civil Engineers, 77(1), 1539-1640. https://doi.org/10.1061/taceat.0002563
#' 
#' Hogg, R. V., & Ledolter, J. (1992). Applied statistics for engineers and physical scientists (2nd int.). Macmillan.
#' 
#' Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. The American Statistician, 50(4), 361-365. https://doi.org/10.2307/2684934
#' 
#' Joarder, A. H., & Firozzaman, M. (2001). Quartiles for discrete data. Teaching Statistics, 23(3), 86-89. https://doi.org/10.1111/1467-9639.00063
#' 
#' Moore, D. S., & McCabe, G. P. (1989). Introduction to the practice of statistics. W.H. Freeman.
#' 
#' SAS. (1990). SAS procedures guide: Version 6 (3rd ed.). SAS Institute.
#' 
#' Snedecor, G. W. (1940). Statistical methods applied to experiments in agriculture and biology (3rd ed.). The Iowa State College Press.
#' 
#' Tukey, J. W. (1977). Exploratory data analysis. Addison-Wesley Pub. Co.
#'  
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
he_quartileIndexing <- function(data, method=c("inclusive", "exclusive", "sas1", "sas4", "hl", "excel", "hf8", "hf9")){
  if (length(method)>1){
    method = "sas1"
  }
  
  data = sort(data)
  
  n = length(data)
  
  if (method=="inclusive"){
    if ((n %% 2) == 0){
      q1Index = (n + 2)/ 4
      q3Index = (3*n + 1)/4}
    else{
      q1Index = (n + 3)/ 4
      q3Index = (3*n + 1)/4}
  }
  else if (method=="exclusive"){
    if ((n %% 2) == 0){
      q1Index = (n + 2)/ 4
      q3Index = (3*n + 1)/4}
    else{
      q1Index = (n + 1)/ 4
      q3Index = (3*n + 3)/4}
  }
  
  else if (method=="sas1"){        
    q1Index = n*1/4
    q3Index = n*3/4}
  else if (method=="sas4"){
    q1Index = (n + 1)*1/4
    q3Index = (n + 1)*3/4}
  else if (method=="hl"){  
    q1Index = n*1/4 + 1/2
    q3Index = n*3/4 + 1/2}
  else if (method=="excel"){
    q1Index = (n - 1)*1/4 + 1
    q3Index = (n - 1)*3/4 + 1}
  else if (method=="hf8"){
    q1Index = (n + 1/3)*1/4 + 1/3
    q3Index = (n + 1/3)*3/4 + 1/3}
  else if (method=="hf9"){
    q1Index = (n + 1/4)*1/4 + 3/8
    q3Index = (n + 1/4)*3/4 + 3/8}
  
  results <- data.frame(q1Index, q3Index)
  
  return (results)
}



