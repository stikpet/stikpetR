#' Stem-and-Leaf Display
#' @description
#' A stem-and-leaf display is defined as: "a method of displaying data in which each observation is split into two parts labelled the ‘stem’ and the ‘leaf’" (Everitt, 2004, p. 362). A diagram that could be used to visualize scale variables, created by Tukey (1972, p. 296).
#' 
#' In some variations of this, the cumulative frequencies are also shown, but currently this function does not provide for that.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/u4XOYJz-aNE) and the visualisation is described at [PeterStatistics.com](https://peterstatistics.com/Terms/Visualisations/stemAndLeafDisplay.html)
#' 
#' @param data list with the numeric data
#' @param key optional factor to use for the stems
#' @returns
#' prints the display in console and returns a dataframe with the stems and leafs
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to create a binned frequency table
#' \code{\link{tab_frequency_bins}}, to create a binned frequency table.
#' 
#' After this you might want some descriptive measures:
#' \code{\link{me_mode_bin}}, for Mode for Binned Data.
#' \code{\link{me_mean}}, for different types of mean.
#' \code{\link{me_variation}}, for different Measures of Quantitative Variation.
#' 
#' Or a perform a test:
#' \code{\link{ts_student_t_os}}, for One-Sample Student t-Test.
#' \code{\link{ts_trimmed_mean_os}}, for One-Sample Trimmed (Yuen or Yuen-Welch) Mean Test.
#' \code{\link{ts_z_os}}, for One-Sample Z Test.
#' 
#' Alternative Visualisations:
#' \code{\link{vi_boxplot_single}}, for a Box (and Whisker) Plot.
#' \code{\link{vi_histogram}}, for a Histogram.
#' 
#' @references 
#' Everitt, B. (2004). *The Cambridge dictionary of statistics* (2nd ed.). Cambridge University Press.
#' 
#' Tukey, J. W. (1972). Some graphic and semigraphic displays. In T. A. Bancroft & S. A. Brown (Eds.), *Statistical Papers in Honor of George W. Snedecor* (pp. 293–316). Iowa State University Press.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' studentDf = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' # Example 1: dataframe
#' ex1 = studentDf[['Gen_Age']]
#' vi_stem_and_leaf(ex1);
#' 
#' # Example 2: Numeric list
#' ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
#' vi_stem_and_leaf(ex2);
#' 
#' @export
vi_stem_and_leaf <- function(data, key=NULL){
  data = na.omit(data)
  
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