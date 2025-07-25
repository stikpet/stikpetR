#' Post-Hoc Residual Test
#' 
#' @param field1 list or dataframe with the first categorical field
#' @param field2 list or dataframe with the second categorical field
#' @param categories1 optional list with order and/or selection for categories of field1
#' @param categories2 optional list with order and/or selection for categories of field2
#' @param residual optional methdod for residual to test. Either "adjusted" (default) or "standardized".
#' 
#' @export
ph_residual <- function(field1, field2, categories1=NULL, categories2=NULL, residual="adjusted"){
  #create the cross table
  ct = tab_cross(field1, field2, categories1, categories2, totals="include")
  
  #basic counts
  nrows = nrow(ct) - 1
  ncols =  ncol(ct) - 1
  n = ct[nrows+1, ncols+1]
  
  res = data.frame()
  ipair=1
  for (i in 1:nrows){
    for (j in 1:ncols){
      res[ipair, 1] = rownames(ct)[i]
      res[ipair, 2] = colnames(ct)[j]
      res[ipair, 3] = ct[i, j]
      
      #expected count
      res[ipair, 4] = ct[nrows + 1, j] * ct[i, ncols + 1] / n
      
      #residual
      res[ipair, 5] = res[ipair, 3] - res[ipair, 4]
      
      #adjusted
      if (residual == "standardized"){
        se = res[ipair, 4]**0.5}
      else if (residual == "adjusted"){
        se = (res[ipair, 4] * (1 - ct[i, ncols + 1] / n) * (1 - ct[nrows + 1, j] / n))**0.5}
      
      res[ipair, 6] = res[ipair, 5] / se
      
      #p-value
      res[ipair, 7] = 2 * (1 - pnorm(abs(res[ipair, 6])))
      
      #adj
      res[ipair, 8] = res[ipair, 7] * nrows * ncols
      if (res[ipair, 8] > 1){
        res[ipair, 8] = 1}
      
      ipair = ipair + 1
    }
  }
  
  #column names
  if (residual == "standardized"){
    resUsed = "st. residual"}
  else if (residual == "adjusted"){
    resUsed = "adj. st. residual"}    
  colnames(res) = c("field1", "field2", "observed", "expected", "residual", resUsed, "p-value", "adj. p-value")
  
  return (res)
  
}



