#' Stem-and-Leaf Display
#' @description
#' A stem-and-leaf display is defined as: "a method of displaying data in which each observation is split into two parts labelled the ‘stem’ and the ‘leaf’" (Everitt, 2004, p. 362). A diagram that could be used to visualize scale variables, created by Tukey (1972, p. 296).
#' 
#' In some variations of this, the cumulative frequencies are also shown, but currently this function does not provide for that.
#' 
#' @param data list with the numeric data
#' @param key optional factor to use for the stems
#' @returns
#' prints the display in console and returns a dataframe with the stems and leafs
#' 
#' @references 
#' Everitt, B. (2004). *The Cambridge dictionary of statistics* (2nd ed.). Cambridge University Press.
#' 
#' Tukey, J. W. (1972). Some graphic and semigraphic displays. In T. A. Bancroft & S. A. Brown (Eds.), *Statistical Papers in Honor of George W. Snedecor* (pp. 293–316). Iowa State University Press.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
vi_stem_and_leaf <- function(data, key=NULL){
  if (is.null(key)){
    key = 10**floor(log10(abs(max(data))))
  }
  
  n = length(data)
  
  dataSorted = sort(na.omit(data))
  
  stems = floor(dataSorted/key)
  
  unformatLeafs = dataSorted - key*stems
  format_str <- paste0("%0", nchar(as.character(key)) - 1, "d")
  leafs <- sapply(1:n, function(i) sprintf(format_str, as.integer(dataSorted[i] - key * stems[i])))
  
  display = data.frame()                
  stemIndex = 1
  stemCurrentValue = stems[1]
  currentLeaf = leafs[1]
  
  for (i in 2:n){
    if (stems[i] == stems[i-1]){
      currentLeaf = paste0(currentLeaf, " ", leafs[i])
    }
    else{
      display[stemIndex ,1] = stems[i-1]
      display[stemIndex ,2] = currentLeaf
      
      currentLeaf = leafs[i]
      stemIndex = stemIndex + 1
      
      stemCurrentValue = stemCurrentValue + 1
      
      while (stemCurrentValue < stems[i]){
        display[stemIndex ,1] = stemCurrentValue
        display[stemIndex ,2] = ""
        
        stemCurrentValue = stemCurrentValue + 1
        stemIndex = stemIndex + 1
      }
    }
  }
  display[stemIndex, 1]  = stems[n]
  display[stemIndex, 2]  = currentLeaf
  
  cat("Stem | Leaf\n") 
  
  for (i in 1:nrow(display)){
    cat(paste0(display[i,1], " | ", display[i,2], "\n"))
  }
  cat(paste0("key: ", stems[1], " | ", leafs[1], " = ", dataSorted[1]))            
  
  return (display)
}