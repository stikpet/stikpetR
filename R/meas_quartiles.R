#' Quartiles / Hinges
#' 
#' @description
#' The quartiles are at quarters of the data (McAlister, 1879, p. 374; Galton, 1881, p. 245). The median is at 50 percent, and the quartiles at 25 and 75 percent. Note that there are five quartiles, the minimum value is the 0-quartile, at 25 percent the first (or lower) quartile, at 50 percent the median a.k.a. the second quartile, at 75 percent the third (or upper) quartile, and the maximum as the fourth quartile.
#' 
#' Tukey (1977) also introduced the term Hinges and sorted the values in a W shape, where the bottom parts of the W are then the hinges.
#' 
#' There are quite a few different methods to determine the quartiles. This function has 19 different ones. See the details for a description.
#' 
#' @param data vector or dataframe with scores as numbers, or if text also provide levels
#' @param levels optional vector with levels in order
#' @param method optional which method to use to calculate quartiles
#' @param indexMethod optional to indicate which type of indexing to use. Either `"sas1"` (default), `"inclusive"`, `"exclusive"`, `"sas4"`, `"excel"`, `"hl"`, `"hf8"`, or `"hf9"`
#' @param q1Frac,q3Frac optional to indicate what type of rounding to use for each quartile. Either `"linear"` (default), `"down"`, `"up"`, `"bankers"`, `"nearest"`, `"halfdown"`, or `"midpoint"`
#' @param q1Int,q3Int optional to indicate the use of the integer or the midpoint method for each quartile. Either `"int"` (default), or `"midpoint"`.
#' 
#' @returns
#' A dataframe with:
#' \item{q1}{the first (lower) quartile}
#' \item{q3}{the third (upper/higher) quartile}
#' \item{q1-text}{the first (lower) quartile as text (only if levels were used)}
#' \item{q3-text}{the third (upper/higher) quartile as text (only if levels were used)}
#' 
#' @details
#' To determine the quartiles a specific indexing method can be used. See \code{\link{he_quartileIndexing}} for details on the different methods to choose from.
#' 
#' Then based on the indexes either linear interpolation or different rounding methods (bankers, nearest, down, up, half-down) can be used, or the midpoint between the two values. If the index is an integer either the integer or the mid point is used. See the \code{\link{he_quartilesIndex}} for details on this.
#' 
#' Note that the rounding method can even vary per quartile, i.e. the one used for the first quartile being different 
#' than the one for the second.
#' 
#' I've come across the following methods:
#' |method|indexing|q1 integer|q1 fractional|q3 integer|q3 fractional|
#' |------|--------|----------|-------------|----------|-------------|
#' |sas1|sas1|use int|linear|use int|linear|
#' |sas2|sas1|use int|bankers|use int|bankers|
#' |sas3|sas1|use int|up|use int|up|
#' |sas5|sas1|midpoint|up|midpoint|up|
#' |hf3b|sas1|use int|nearest|use int|halfdown|
#' |sas4|sas4|use int|linear|use int|linear|
#' |ms|sas4|use int|nearest|use int|halfdown|
#' |lohninger|sas4|use int|nearest|use int|nearest|
#' |hl2|hl|use int|linear|use int|linear|
#' |hl1|hl|use int|midpoint|use int|midpoint|
#' |excel|excel|use int|linear|use int|linear|
#' |pd2|excel|use int|down|use int|down|
#' |pd3|excel|use int|up|use int|up|
#' |pd4|excel|use int|halfdown|use int|nearest|
#' |pd5|excel|use int|midpoint|use int|midpoint|
#' |hf8|hf8|use int|linear|use int|linear|
#' |hf9|hf9|use int|linear|use int|linear|
#' 
#' The following values can be used for the *method* parameter:
#' \itemize{
#' \item inclusive = tukey =hinges = vining. (Tukey, 1977, p. 32; Siegel & Morgan, 1996, p. 77; Vining, 1998, p. 44).
#' \item exclusive = jf. (Moore & McCabe, 1989, p. 33; Joarder & Firozzaman, 2001, p. 88).
#' \item sas1 = parzen = hf4 = interpolated_inverted_cdf = maple3 = r4. (Parzen, 1979, p. 108; SAS, 1990, p. 626; Hyndman & Fan, 1996, p. 363)
#' \item sas2 = hf3 = r3. (SAS, 1990, p. 626; Hyndman & Fan, 1996, p. 362)
#' \item sas3 = hf1 = inverted_cdf = maple1 = r1 (SAS, 1990, p. 626; Hyndman & Fan, 1996, p. 362)
#' \item sas4 = hf6 = minitab = snedecor = weibull = maple5 = r6 (Hyndman & Fan, 1996, p. 363; Weibull, 1939, p. ?; Snedecor, 1940, p. 43; SAS, 1990, p. 626)
#' \item sas5 = hf2 = CDF = averaged_inverted_cdf = r2 (SAS, 1990, p. 626; Hyndman & Fan, 1996, p. 362)
#' \item hf3b = closest_observation 
#' \item ms (Mendenhall & Sincich, 1992, p. 35)
#' \item lohninger (Lohninger, n.d.)
#' \item hl1 (Hogg & Ledolter, 1992, p. 21)
#' \item hl2 = hf5 = Hazen = maple4 = r5 (Hogg & Ledolter, 1992, p. 21; Hazen, 1914, p. ?)
#' \item maple2
#' \item excel = hf7 = pd1 = linear = gumbel = maple6 = r7 (Hyndman & Fan, 1996, p. 363; Freund & Perles, 1987, p. 201; Gumbel, 1939, p. ?)
#' \item pd2 = lower
#' \item pd3 = higher
#' \item pd4 = nearest
#' \item pd5 = midpoint
#' \item hf8 = median_unbiased = maple7 = r8 (Hyndman & Fan, 1996, p. 363)
#' \item hf9 = normal_unbiased = maple8 = r9 (Hyndman & Fan, 1996, p. 363)
#' }
#' 
#' *hf* is short for Hyndman and Fan who wrote an article showcasing many different methods, 
#' *hl* is short for Hog and Ledolter, *ms* is short for Mendenhall and Sincich, *jf* is short for Joarder and Firozzaman. 
#' *sas* refers to the software package SAS, *maple* to Maple, *pd* to Python's pandas library, and *r* to R.
#' 
#' The names *linear*, *lower*, *higher*, *nearest* and *midpoint* are all used by pandas quantile function and numpy 
#' percentile function. Numpy also uses *inverted_cdf*, *averaged_inverted_cdf*, *closest_observation*, 
#' *interpolated_inverted_cdf*, *hazen*, *weibull*, *median_unbiased*, and *normal_unbiased*. 
#' 
#' 
#' @references 
#' Freund, J. E., & Perles, B. M. (1987). A new look at quartiles of ungrouped data. *The American Statistician, 41*(3), 200–203. https://doi.org/10.1080/00031305.1987.10475479
#' 
#' Galton, F. (1881). Report of the anthropometric committee. *Report of the British Association for the Advancement of Science, 51*, 225–272.
#' 
#' Gumbel, E. J. (1939). La Probabilité des Hypothèses. *Compes Rendus de l’ Académie des Sciences, 209*, 645–647.
#' 
#' Hazen, A. (1914). Storage to be provided in impounding municipal water supply. *Transactions of the American Society of Civil Engineers, 77*(1), 1539–1640. https://doi.org/10.1061/taceat.0002563
#' 
#' Hogg, R. V., & Ledolter, J. (1992). *Applied statistics for engineers and physical scientists* (2nd int.). Macmillan.
#' 
#' Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. *The American Statistician, 50*(4), 361–365. https://doi.org/10.2307/2684934
#' 
#' Joarder, A. H., & Firozzaman, M. (2001). Quartiles for discrete data. *Teaching Statistics, 23*(3), 86–89. https://doi.org/10.1111/1467-9639.00063
#' 
#' Langford, E. (2006). Quartiles in elementary statistics. *Journal of Statistics Education, 14*(3), 1–17. https://doi.org/10.1080/10691898.2006.11910589
#' 
#' Lohninger, H. (n.d.). Quartile. Fundamentals of Statistics. Retrieved April 7, 2023, from http://www.statistics4u.com/fundstat_eng/cc_quartile.html
#' 
#' McAlister, D. (1879). The law of the geometric mean. *Proceedings of the Royal Society of London, 29*(196–199), 367–376. https://doi.org/10.1098/rspl.1879.0061
#' 
#' Mendenhall, W., & Sincich, T. (1992). *Statistics for engineering and the sciences* (3rd ed.). Dellen Publishing Company.
#' 
#' Moore, D. S., & McCabe, G. P. (1989). *Introduction to the practice of statistics*. W.H. Freeman.
#' 
#' Parzen, E. (1979). Nonparametric statistical data modeling. *Journal of the American Statistical Association, 74*(365), 105–121. https://doi.org/10.1080/01621459.1979.10481621
#' 
#' SAS. (1990). SAS procedures guide: Version 6 (3rd ed.). SAS Institute.
#' 
#' Siegel, A. F., & Morgan, C. J. (1996). *Statistics and data analysis: An introduction* (2nd ed.). J. Wiley.
#' 
#' Snedecor, G. W. (1940). *Statistical methods applied to experiments in agriculture and biology* (3rd ed.). The Iowa State College Press.
#' 
#' Tukey, J. W. (1977). *Exploratory data analysis*. Addison-Wesley Pub. Co.
#' 
#' Vining, G. G. (1998). *Statistical methods for engineers*. Duxbury Press.
#' 
#' Weibull, W. (1939).* The phenomenon of rupture in solids*. Ingeniörs Vetenskaps Akademien, 153, 1–55.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @examples 
#' ex8 = c(1, 2, 3, 4, 5, 6, 7, 8)
#' me_quartiles(ex8)
#' test = c("a", "b", "c", "d", "e", "f")
#' testCoding = c("a", "b", "c", "d", "e", "f")
#' me_quartiles(test, levels=testCoding, method="excel")
#' 
#' 
#' @export
me_quartiles <- function(data, levels=NULL, 
                         method="own", 
                         indexMethod=c("inclusive", "exclusive", "sas1", "sas4", "hl", "excel", "hf8", "hf9"), 
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
  
  if (is.null(levels)){
    dataN = data}
  else{
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    dataN = as.numeric(myFieldOrd)
  }
  
  #Sort the data
  dataN = sort(dataN)
  
  #alternative namings
  if (method %in% c("inclusive", "tukey", "vining", "hinges")){
    method="inclusive"}
  else if (method %in% c("exclusive", "jf")){
    method ="exclusive"}
  else if (method %in% c("cdf", "sas5", "hf2", "averaged_inverted_cdf", "r2")){
    method = "sas5"}
  else if (method %in% c("sas4", "minitab", "hf6", "weibull", "maple5", "r6")){
    method = "sas4"}
  else if (method %in% c("excel", "hf7", "pd1", "linear", "gumbel", "maple6", "r7")){
    method = "excel"}
  else if (method %in% c("sas1", "parzen", "hf4", "interpolated_inverted_cdf", "maple3", "r4")){
    method = "sas1"}
  else if (method %in% c("sas2", "hf3", "r3")){
    method = "sas2"}
  else if (method %in% c("sas3", "hf1", "inverted_cdf", "maple1", "r1")){
    method = "sas3"}
  else if (method %in% c("hf3b", "closest_observation")){
    method = "hf3b"}
  else if (method %in% c("hl2", "hazen", "hf5", "maple4")){
    method = "hl2"}
  else if (method %in% c("np", "midpoint", "pd5")){
    method = "pd5"}
  else if (method %in% c("hf8", "median_unbiased", "maple7", "r8")){
    method = "hf8"}
  else if (method %in% c("hf9", "normal_unbiased", "maple8", "r9")){
    method = "hf9"}
  else if (method %in% c("pd2", "lower")){
    method = "pd2"}
  else if (method %in% c("pd3", "higher")){
    method = "pd3"}
  else if (method %in% c("pd4", "nearest")){
    method = "pd4"}
  
  #settings
  settings = c(indexMethod, q1Frac, q1Int, q3Frac, q3Int)
  if (method=="inclusive"){
    settings = c("inclusive", "linear","int","linear","int")}
  else if (method=="exclusive"){
    settings = c("exclusive", "linear","int","linear","int")}
  else if (method=="sas1"){
    settings = c("sas1","linear","int","linear","int")}
  else if (method=="sas2"){
    settings = c("sas1","bankers","int","bankers" ,"int")}
  else if (method=="sas3"){
    settings = c("sas1","up","int","up","int")}
  else if (method=="sas5"){
    settings = c("sas1","up","midpoint","up","midpoint")}
  else if (method=="sas4"){    
    settings = c("sas4","linear", "int","linear","int")}
  else if (method=="ms"){ 
    settings = c("sas4", "nearest","int", "halfdown","int")}
  else if (method=="lohninger"){
    settings = c("sas4", "nearest", "int","nearest","int")}
  else if (method=="hl2"){
    settings = c("hl", "linear", "int","linear","int")}
  else if (method=="hl1"){
    settings = c("hl", "midpoint","int", "midpoint","int")}
  else if (method=="excel"){
    settings = c("excel", "linear","int","linear", "int")}
  else if (method=="pd2"){
    settings = c("excel", "down", "int", "down","int")}
  else if (method=="pd3"){
    settings = c("excel", "up","int","up","int")}
  else if (method=="pd4"){
    settings = c("excel", "halfdown",  "int","nearest", "int")}
  else if (method=="hf3b"){
    settings = c("sas1", "nearest","int","halfdown","int")}
  else if (method=="pd5"){
    settings = c("excel", "midpoint","int","midpoint","int")}
  else if (method=="hf8"){
    settings = c("hf8", "linear","int","linear", "int")}
  else if (method=="hf9"){
    settings = c("hf9", "linear","int","linear", "int")}
  else if (method=="maple2"){
    settings = c("hl", "down","int","down", "int")}
  
  indexes = he_quartilesIndex(dataN, settings[1], settings[2], settings[3], settings[4], settings[5])
  q1 = indexes[1,1]
  q3 = indexes[1,2]
  
  if (is.null(levels)){        
    results = data.frame(q1, q3)}
  else{
    if (q1 == round(q1)){
      q1Text = levels[q1]}
    else{
      q1Text = paste0("between ", levels[floor(q1)], " and ", levels[ceiling(q1)])}
    
    if (q3 == round(q3)){
      q3Text = levels[q3]}
    else{
      q3Text = paste0("between ", levels[floor(q3)], " and ", levels[ceiling(q3)])}
    
    results = data.frame(q1, q3, q1Text, q3Text)
  }
  
  return (results)
}
