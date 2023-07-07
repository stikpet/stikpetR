#' Cleveland Dot Plot
#' 
#' @description 
#' A Cleveland dot plot (Cleveland & McGill, 1987) is a bar chart where instead of bars 
#' a dot is placed at the center of the top of the bar (and then the bars removed). 
#' It is a dot plot only showing the top dot.This requires less ink. 
#' 
#' The function simply uses the dotplot() function from the lattice library.
#' 
#' A video on (Cleveland) dot plots is available [here](https://youtu.be/qs1nh0CMiIY).
#' 
#' @param data the data from which to create the plot
#' @param size the size of the dots (default is 2)
#' @return chart the Cleveland dot plot
#' 
#' @references 
#' Cleveland, W. S., & McGill, R. (1984). Graphical perception: Theory, experimentation, and application to the development of graphical methods. *Journal of the American Statistical Association, 79*(387), 531–554. https://doi.org/10.2307/2288400
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @importFrom lattice dotplot
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' vi_cleveland_dot_plot(ex1);
#' 
#' #Example 2: a list
#' ex2 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' vi_cleveland_dot_plot(ex2);
#' 
#' @export
vi_cleveland_dot_plot <- function(data, size=2){
  
  chart = lattice::dotplot(data, horizontal = FALSE, cex=size)
  
  return(chart)
}