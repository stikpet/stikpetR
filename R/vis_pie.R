#' Pie Chart
#' 
#' @description 
#' A pie-chart is a “graphic display in which a circle is cut into wedges with the area of 
#' each wedge being proportional to the percentage of cases in the category represented by 
#' that wedge” (Zedeck, 2014, p. 260). 
#' 
#' A video on pie charts is available [here](https://youtu.be/e6JtJsh-6iw).
#' 
#' @param data the data for which to create a pie-chart from
#' @param labels what to show besides the labels Either "count" (default), "percent", "none", or "both"
#' 
#' @return chart the pie chart
#' 
#' @details 
#' It is possible to either show only the labels (label="none"), the counts (label="counts"), 
#' the percentages (label="percent"), or both count and percent (label="both").
#' 
#' The function uses the basic R's graphics library *pie* function, rotated and counter clockwise.
#' 
#' The pie-chart is quite popular and often used, but actually has a few disadvantages. 
#' It can only show relative frequencies. To show other frequencies the numbers themselves 
#' have to be added. A circle has 360 degrees, equal to 100%. So by multiplying the relative 
#' frequencies with 360, the degrees for each category can be found. This means that visually 
#' the pie-chart can only show the relative frequencies.
#' 
#' Another disadvantage is when the relative frequencies are close to each other, 
#' the differences are not easily seen in a circle diagram.
#' 
#' As a third disadvantage, when there are many categories the circle diagram will 
#' look very busy and not easily to read.
#' 
#' People also have more difficulty with comparing areas and angles 
#' (what you do when looking at a pie-chart) than comparing heights (what is done with a bar-chart).
#' 
#' Also often a 3D effect is added, but this actually makes comparisons of the slices even more difficult.
#' 
#' he earliest found circle diagram is found on the inlay of a book by William Playfair (1801).
#' The name 'pie chart' might come from a misspelling of the word Pi. 
#' Pi is often associated with a circle. It might also simply come from the resemblances 
#' with a pie (as in apple-pie). 
#' However Srivastava and Rego (2011) put forward another belief that it is named after a 
#' royal French cook Pie, who served dishes in a pie-chart shape.
#' 
#' @references 
#' Playfair, W. (1801). *The statistical breviary: Shewing the resources of every state and kingdom*. T. Bensley. http://archive.org/details/statisticalbrev00playgoog
#' 
#' Srivastava, T. N., & Rego, S. (2011). *Business research methodology*. Tata McGraw-Hill.
#' 
#' Zedeck, S. (Ed.). (2014). *APA dictionary of statistics and research methods*. American Psychological Association.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' vi_pie(ex1);
#' vi_pie(ex1, labels="percent");
#' vi_pie(ex1, labels="none");
#' vi_pie(ex1, labels="both");
#' 
#' #Example 2: a list
#' ex2 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' vi_pie(ex2);
#' 
#' @export
vi_pie <- function(data, labels=c("count", "percent", "both", "none")){
  
  if (length(labels)>1) {labels="count"}
  
  f = table(data)
  dataLabels = names(f)
  
  if (labels=="count" || labels=="both") {
    dataLabels = paste(dataLabels, unname(f)) 
  }
  
  if (labels=="percent" || labels=="both") {
    pct = prop.table(f)*100
    percLabels = paste(round(pct), "%", sep = "")
    dataLabels <- paste(dataLabels, percLabels)
  }
  
  if (labels %in% c("count", "percent", "both")) {
    chart = pie(f, labels = dataLabels, init.angle = 90, clockwise = TRUE)
  }
  else{
    chart = pie(f, init.angle = 90, clockwise = TRUE)
  }
  
  return(chart)
  
}