#' Post-Hoc Column Proportion Test
#' 
#' @param field1 list or dataframe with the first categorical field
#' @param field2 list or dataframe with the second categorical field
#' @param categories1 optional list with order and/or selection for categories of field1
#' @param categories2 optional list with order and/or selection for categories of field2
#' @param seMethod optional methdod for standard error. Either "spss" (default) or "marascuilo".
#' 
#' @export
ph_column_proportion <- function(field1, field2, categories1=NULL, categories2=NULL, seMethod = "spss"){
  #create the cross table
  ct = tab_cross(field1, field2, categories1, categories2, totals="include")
  
  #basic counts
  nrows = nrow(ct) - 1
  ncols =  ncol(ct) - 1
  n = ct[nrows+1, ncols+1]
  
  res = data.frame()
  ipair=1
  for (i in 1:nrows){
    for (j in 1:(ncols-1)){
      for (k in (j+1):ncols){
        res[ipair, 1] = rownames(ct)[i]
        res[ipair, 2] = colnames(ct)[j]
        res[ipair, 3] = colnames(ct)[k]
        
        #column proportions
        res[ipair, 4] = ct[i, j] / ct[nrows + 1, j]
        res[ipair, 5] = ct[i, k] / ct[nrows + 1, k]
        
        #residual
        res[ipair, 6] = res[ipair, 4] - res[ipair, 5]
        
        #statistic
        if (seMethod == "spss"){
          phat = (ct[nrows + 1, j] * res[ipair, 4] + ct[nrows + 1, k] * res[ipair, 5]) / (ct[nrows + 1, j] + ct[nrows + 1, k])
          se = (phat * (1 - phat) * (1 / ct[nrows + 1, j] + 1 / ct[nrows + 1, k]))**0.5}
        else {
          se = (1 / ct[nrows + 1, j] * (res[ipair, 4] * (1 - res[ipair, 5])) + 1 / ct[nrows + 1, k] * (res[ipair, 5] * (1 - res[ipair, 4])))**0.5}
        
        res[ipair, 7] = res[ipair, 6] / se
        
        #p-value
        res[ipair, 8] = 2 * (1 - pnorm(abs(res[ipair, 7])))
        
        #adj
        res[ipair, 9] = res[ipair, 8] * ncols * (ncols - 1) / 2
        if (res[ipair, 9] > 1){
          res[ipair, 9] = 1}
        
        ipair = ipair + 1
      }
    }
  }
  
  colnames(res) = c("field1", "field2-1", "field2-2", "col. prop. 1", "col. prop. 2", "difference", "z-value", "p-value", "adj. p-value")
  
  return (res)
  
}