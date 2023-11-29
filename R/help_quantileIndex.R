#' Quantile Numeric Based on Index
#' 
#' Helper function for \code{\link{me_quartiles}} to return the quartile as a number of the first and third quartile with different methods 
#' of rounding.
#' 
#' @param data dataframe with scores as numbers
#' @param k : number of quantiles
#' @param indexMethod optional to indicate which type of indexing to use
#' @param qLfrac optional to indicate what type of rounding to use for quantiles below median
#' @param qLint optional to indicate the use of the integer or the midpoint method for quantiles below median
#' @param qHfrac optional to indicate what type of rounding to use for quantiles above median
#' @param qHint optional to indicate the use of the integer or the midpoint method for quantiles above median
#' 
#' @returns
#' A vector with the quantiles
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
he_quantilesIndex <-function(data, 
                             k=4, 
                             indexMethod=c("sas1", "sas4", "hl", "excel", "hf8", "hf9"), 
                             qLfrac=c("linear", "down", "up", "bankers", "nearest", "halfdown", "midpoint"), 
                             qLint=c("int", "midpoint"), 
                             qHfrac=c("linear", "down", "up", "bankers", "nearest", "halfdown", "midpoint"), 
                             qHint=c("int", "midpoint")){
  
  # Set defaults
  if (length(indexMethod)>1){indexMethod = "sas1"}
  if (length(qLfrac)>1){qLfrac = "linear"}
  if (length(qLint)>1){qLint = "int"}
  if (length(qHfrac)>1){qHfrac = "linear"}
  if (length(qHint)>1){qHint = "int"}
  
  #SORT THE DATA
  data = sort(data)
  
  #Returns quartiles based on index
  n = length(data)
  qIndexes = he_quantileIndexing(data, k=k, indexMethod)
  qs=c()
  for (i in 1:(k+1)){
    iq = qIndexes[i]
    if (iq < 0.5){
      if (round(iq) == iq){
        # index is integer
        if (qLint == "int"){q1 = iq}
        else if (qLint == "midpoint"){q1 = iq + 1/2}
      }
      else{
        # index has fraction
        if (qLfrac == "linear"){q1 = iq}
        else if (qLfrac == "down"){q1 = floor(iq)}
        else if (qLfrac == "up"){q1 = ceiling(iq)}
        else if (qLfrac == "bankers"){q1 = round(iq)}
        else if (qLfrac == "nearest"){q1 = as.integer(iq + 0.5)}
        else if (qLfrac == "halfdown"){
          if (iq + 0.5 == round(iq + 0.5)){
            q1 = floor(iq)}
          else{q1 = round(iq)}
        }
        else if (qLfrac == "midpoint"){q1 = (floor(iq) + ceiling(iq)) / 2}
      }
    }
    
    else{
      if (round(iq) == iq){
        # index is integer
        if (qHint == "int"){q1 = iq}
        else if (qHint == "midpoint"){q1 = iq + 1/2}
      }
      else{
        # index has fraction
        if (qHfrac == "linear"){q1 = iq}
        else if (qHfrac == "down"){q1 = floor(iq)}
        else if (qHfrac == "up"){q1 = ceiling(iq)}
        else if (qHfrac == "bankers"){q1 = round(iq)}
        else if (qHfrac == "nearest"){q1 = as.integer(iq + 0.5)}
        else if (qHfrac == "halfdown"){
          if (iq + 0.5 == round(iq + 0.5)){
            q1 = floor(iq)}
          else{q1 = round(iq)}
        }
        else if (qHfrac == "midpoint"){q1 = (floor(iq) + ceiling(iq)) / 2}
      }
    }
    
    if (q1 < 1) {q1 = 1}
    else if (q1 > n){q1 = n}
    
    q1i = q1
    q1iLow = floor(q1i)
    q1iHigh = ceiling(q1i)
    
    if (q1iLow==q1iHigh){q1 = data[q1iLow]}
    else{
      #Linear interpolation:
      q1 = data[q1iLow] + (q1i - q1iLow)/(q1iHigh - q1iLow)*(data[q1iHigh] - data[q1iLow])
    }
    
    qs = append(qs, q1)
  }
  
  return (qs)
}