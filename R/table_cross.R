#' Cross Table / Contingency Table
#' 
#' @description
#' A contingency table can be defined as "tables arising when observations on a number of categorical variables are cross-classified" (Everitt, 2004, p.89).
#' 
#' There are quite a few variations on the name for this type of table. Perhaps the oldest name is actually contingency table, which was the name Pearson (1904, p. 34) gave to them. Another popular name is cross tabulation (Upton & Cook, 2002, p. 79), but also cross classification table (Zekeck, 2014, p. 71) and bivariate frequency table (Porkess, 1988, p. 48) are used. The one I used cross table which can for example be found in Newbold et al. (2013, p. 9) or Sa (2007, p. 52).
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param order1 : optional list with order for categories of field1
#' @param order2 : optional list with order for categories of field2
#' @param percent : optional which percentages to show. Either "none" (default), "all", "row", "column"
#' @param totals : optional to add margin totals. Either "exclude" (default), or "include"
#' 
#' @returns
#' dataframe : the cross table
#' 
#' @references
#' Everitt, B. (2004). *The Cambridge dictionary of statistics* (2nd ed.). Cambridge University Press.
#' Newbold, P., Carlson, W. L., & Thorne, B. (2013). *Statistics for business and economics* (8th ed). Pearson.
#' 
#' Pearson, K. (1904). *Contributions to the Mathematical Theory of Evolution*. XIII. On the theory of contingency and its relation to association and normal correlation. Dulau and Co.
#' 
#' Porkess, R. (1988). *Dictionary of statistics*. Collins.
#' 
#' Sa, J. P. M. de. (2007). *Applied statistics: Using SPSS, Statistica, MATLAB, and R* (2nd ed.). Springer.
#' 
#' Upton, G., & Cook, I. (2002). *Oxford: Dictionary of statistics*. Oxford University Press.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' tab_cross(df1[['mar1']], df1[['sex']], percent="column", totals="include")
#' 
#' orderR = c("DIVORCED", "WIDOWED", "SEPARATED", "MARRIED", "NEVER MARRIED")
#' orderC = c("MALE", "FEMALE")
#' tab_cross(df1[['mar1']], df1[['sex']], order1=orderR, order2=orderC)
#' 
#' order = c("Not scientific at all", "Not too scientific", "Pretty scientific", "Very scientific")
#' tab_cross(df1[['mar1']], df1[['accntsci']], order2=order)
#' 
#' @export
tab_cross <- function(field1, 
                      field2, 
                      order1=NULL, 
                      order2=NULL, 
                      percent=c(NULL, "all", "row", "column"), 
                      totals="exclude"){
  
  #Set default value
  if (length(percent) > 1) {percent="none"}
  
  #create basic frequency cross table
  tab = table(field1, field2)
  
  #reorder rows if order provided
  if (!is.null(order1)){
    if(!is.null(names(order1))){
      for (i in 1:length(order1)){field1[field1 == unname(order1[i])] = names(order1)[i]}
      order1 = names(order1)
      }
    tab = tab[order1, ]
  }
  
  #reorder columns if order provided
  if (!is.null(order2)){
    if(!is.null(names(order2))){
      for (i in 1:length(order2)){field2[field2 == unname(order2[i])] = names(order2)[i]}
      order2 = names(order2)
    }
    tab = tab[, order2]  
  }
  
  #add margin totals if requested
  f = 1
  if (totals!="exclude"){
    tab = addmargins(tab)
    f = 2}
  
  #convert to proper percentages if requested
  if (percent=="all"){tab = prop.table(tab)*f**2*100}
  else if (percent=="row"){tab = prop.table(tab, 1)*f*100}
  else if( percent=="column"){tab = prop.table(tab, 2)*f*100}
  
  
  return (tab)
}



