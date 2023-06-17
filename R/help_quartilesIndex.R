#' Quartile Numeric Based on Index
#' 
#' Helper function for \code{\link{me_quartiles}} to return the quartile as a number of the first and third quartile with different methods 
#' of rounding.
#' 
#' @param data dataframe with scores as numbers
#' @param indexMethod optional to indicate which type of indexing to use
#' @param q1Frac optional to indicate what type of rounding to use for first quarter
#' @param q1Int optional to indicate the use of the integer or the midpoint method for first quarter
#' @param q3Frac optional to indicate what type of rounding to use for third quarter
#' @param q3Int optional to indicate the use of the integer or the midpoint method for third quarter
#' 
#' @returns
#' A dataframe with:
#' \item{q1}{the first (lower) quartile}
#' \item{q3}{the third (upper/higher) quartile}
#' 
#' @details 
#' If **the index is an integer** often that integer will be used to find the corresponding value in the sorted data. 
#' However, in some rare methods they argue to take the midpoint between the found index and the next one, i.e. to use:
#' \deqn{iQ_i = iQ_i + \frac{1}{2}}
#' 
#' If the index has a fractional part, we could use linear interpolation. It can be written as:
#' \deqn{X\left[\lfloor iQ_i \rfloor\right] + \frac{iQ_i - \lfloor iQ_i \rfloor}{\lceil iQ_i \rceil - \lfloor iQ_i \rfloor} \times \left(X\left[\lceil iQ_i \rceil\right] - X\left[\lfloor iQ_i \rfloor\right]\right)}
#' Where:
#' \itemize{
#' \item \eqn{X\left[x\right]} is the x-th score of the sorted scores 
#' \item  \eqn{\lfloor\dots\rfloor} is the function to always round down
#' \item  \eqn{\lceil\dots\rceil} is the function to always round up
#' }
#' 
#' Or we can use 'rounding'. But there are different versions of rounding. Besides the already mentioned round 
#' down and round up versions:
#' \itemize{
#' \item \eqn{\lfloor\dots\rceil} to indicate rounding to the nearest even integer. A value of 2.5 gets rounded to 2, while 1.5 also gets rounded to 2. This is also referred to as *bankers* method.
#' \item  \eqn{\left[\dots\right]} to indicate rounding to the nearest integer. A value that ends with .5 is then always rounded up.
#' \item  \eqn{\left< \dots\right>} to indicate to round a value ending with .5 always down
#' }
#' or even use the midpoint again i.e.:
#' \deqn{\frac{\lfloor iQ_i \rfloor + \lceil iQ_i \rceil}{2}}
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
he_quartilesIndex <-function(data, indexMethod=c("inclusive", "exclusive", "sas1", "sas4", "hl", "excel", "hf8", "hf9"), 
                             q1Frac=c("linear", "down", "up", "bankers", "nearest", "halfdown", "midpoint"), 
                             q1Int=c("int", "midpoint"), 
                             q3Frac=c("linear", "down", "up", "bankers", "nearest", "halfdown", "midpoint"), 
                             q3Int=c("int", "midpoint")){
  # Set defaults
  if (length(indexMethod)>1){indexMethod = "sas1"}
  if (length(q1Frac)>1){q1Frac = "linear"}
  if (length(q1Int)>1){q1Int = "int"}
  if (length(q3Frac)>1){q3Frac = "linear"}
  if (length(q3Int)>1){q3Int = "int"}
  
  #SORT THE DATA
  data = sort(data)
  
  #Returns quartiles based on index
  n = length(data)
  qIndexes = he_quartileIndexing(data, indexMethod)
  iq1 = qIndexes[1,1]
  iq3 = qIndexes[1,2]
  
  if (round(iq1) == iq1){
    # index is integer
    if (q1Int == "int"){
      q1 = iq1}
    else if (q1Int == "midpoint"){
      q1 = iq1 + 1/2}
  }
  else{
    # index has fraction
    if (q1Frac == "linear"){
      q1 = iq1}
    else if (q1Frac == "down"){
      q1 = floor(iq1)}
    else if (q1Frac == "up"){
      q1 = ceiling(iq1)}
    else if (q1Frac == "bankers"){
      q1 = round(iq1)}
    else if (q1Frac == "nearest"){
      q1 = as.integer(iq1 + 0.5)}
    else if (q1Frac == "halfdown"){
      if (iq1 + 0.5 == round(iq1 + 0.5)){
        q1 = floor(iq1)}
      else{
        q1 = round(iq1)}
    }
    else if (q1Frac == "midpoint"){
      q1 = (floor(iq1) + ceiling(iq1)) / 2}
  }
  
  q1i = q1
  q1iLow = floor(q1i)
  q1iHigh = ceiling(q1i)
  
  if (q1iLow==q1iHigh){
    q1 = data[q1iLow]}
  else{
    #Linear interpolation:
    q1 = data[q1iLow] + (q1i - q1iLow)/(q1iHigh - q1iLow)*(data[q1iHigh] - data[q1iLow])
  }
  
  
  if (round(iq3) == iq3){
    # index is integer
    if (q3Int == "int"){
      q3 = iq3}
    else if (q3Int == "midpoint"){
      q3 = iq3 + 1/2}
  }    
  else{
    # index has fraction
    if (q3Frac == "linear"){
      q3 = iq3}
    else if (q3Frac == "down"){
      q3 = floor(iq3)}
    else if (q3Frac == "up"){
      q3 = ceiling(iq3)}
    else if (q3Frac == "bankers"){
      q3 = round(iq3)}
    else if (q3Frac == "nearest"){
      q3 = as.integer(iq3 + 0.5)}
    else if (q3Frac == "halfdown"){
      if (iq3 + 0.5 == round(iq3 + 0.5)){
        q3 = floor(iq3)}
      else{
        q3 = round(iq3)}
    }
    else if (q3Frac == "midpoint"){
      q3 = (floor(iq3) + ceiling(iq3)) / 2}
  }
  
  q3i = q3
  q3iLow = floor(q3i)
  q3iHigh = ceiling(q3i)
  
  
  
  if (q3iLow==q3iHigh){
    q3 = data[q3iLow]}
  else{
    #Linear interpolation:
    q3 = data[q3iLow] + (q3i - q3iLow)/(q3iHigh - q3iLow)*(data[q3iHigh] - data[q3iLow])}
  
  results <- data.frame(q1, q3)
  
  return (results)
}