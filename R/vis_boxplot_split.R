#' Split Box Plot
#' @description 
#' Based on a categorical field the scores for each category are plotted in a separate boxplot and each of them is placed underneath each other.
#' 
#' See **vi_boxplot_single()** for more details on boxplots.
#' 
#' @param catField list or dataframe with the categories
#' @param scaleField list or dataframe with the scores
#' @param categories optional list with categories to use
#' @param ... other parameters for use in boxplot function
#' 
#' @returns
#' The split boxplot
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
vi_boxplot_split <- function(catField, scaleField, categories=NULL, ...){
  #get the original names of the variables
  if (is.null(names(catField))){
    argnames = sys.call()
    argnames = unlist(lapply(argnames[-1], as.character))
    catName=argnames[1]}
  else {
    catName=names(catField)}
  if (is.null(names(scaleField))){
    argnames = sys.call()
    argnames = unlist(lapply(argnames[-1], as.character))
    scaleName=argnames[2]}
  else {
    scaleName=names(scaleField)}
  
  #create a dataframe
  df = data.frame(catField, scaleField)
  df = na.omit(df)
  colnames(df) = c("category", "score")
  
  if (!is.null(categories)){
    df = df[df[, 1] %in% categories, ]
  }
  
  #create the split boxplot
  boxplot(score~category, data=df, xlab=scaleName, ylab=catName, horizontal=TRUE, ...)  
}