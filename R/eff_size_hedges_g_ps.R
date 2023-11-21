#' Hedges g (Paired Samples)
#' @description 
#' Effect size measure for paired samples. This is very similar as Hedges g for independent samples.
#' 
#' @param field1 the scores on the first variable
#' @param field2 the scores on the second variable
#' @param dmu optional the difference according to null hypothesis (default is 0)
#' @param appr approximation to use (see details), default is "none"
#' @param within boolean to use a correction for correlated pairs
#' 
#' @returns 
#' A dataframe with:
#' \item{g}{the Hedges g value}
#' \item{version}{version that was used}
#' 
#' @details 
#' The formula used is the same as for Cohen d_z
#' 
#' The same corrections can then be applied as for the independent samples version. See es_hedges_g_is() for details.
#'  
#' **Alternatives**
#' 
#' *library(effsize)*
#' 
#' datF = na.omit(data.frame(field1, field2))
#' 
#' cohen.d(datF$field1, datF$field2, paired=TRUE, within=TRUE, hedges.correction=TRUE)
#' 
#' cohen.d(datF$field1, datF$field2, paired=TRUE, within=FALSE, hedges.correction=TRUE)
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_hedges_g_ps <- function(field1, field2, dmu=0, appr=c("none", "hedges", "durlak", "xue"), within=TRUE){
  
  if (length(appr)>1) {
    appr="none"
  }
  
  datF = na.omit(data.frame(field1, field2))
  n = nrow(datF)
  
  d = datF$field1 - datF$field2
  dAvg = mean(d)
  
  s = sd(d)
  
  comment=""
  if (within) {
    s = s/sqrt(2*(1 - cor(datF$field1, datF$field2)))
    comment=", with correlation correction"
  }
  
  dz = (dAvg - dmu)/s
  df = n - 1
  m = df/2
  
  #check if exact chosen and can be calculated
  if (appr=="none" && m > 171){
    warning("Exact method could not be computed due to large sample size, Hedges approximation used instead")
    appr="hedges"    
  }
  
  
  if (appr=="none" && m <172){
    # Hedges g (exact)
    g = dz*gamma(m)/(gamma(m-0.5)*sqrt(m))
    comment = paste0("exact", comment)
  }
  else if(appr=="hedges"){
    # Hedges approximation
    g = dz*(1-3/(4*df-1))
    comment = paste0("Hedges approximation", comment)
  }
  else if(appr=="durlak"){
    # Durlak (2009, p. 927) approximation:
    g = dz*(n-3)/(n-2.25)*sqrt((n-2)/n)
    comment = paste0("Durlak approximation", comment)
  }
  else if(appr=="xue"){
    # Xue (2020, p. 3) approximation:
    g = dz*(1 - 9/df + 69/(2*df^2) - 72/(df^3) + 687/(8*df^4) - 441/(8*df^5) + 247/(16*df^6))^(1/12)
    comment = paste0("Xue approximation", comment)
  }

  #prepare results
  res <- data.frame(g, comment)
  colnames(res) <- c('g','version')
  
  return(res)
  
}