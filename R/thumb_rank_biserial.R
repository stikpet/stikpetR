#' Rule-of-Thumb for Rank Biserial Correlation
#' 
#' @description
#' Simple function to use a rule-of-thumb for the Rank Biserial Correlation.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/31jr7jUCni4) and the measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Correlations/RankBiserialCorrelation.html)
#' 
#' 
#' @param rb the rank-biserial correlation value
#' @param version optional version of rank-biserial that was used. Either "glass" (default) or "cureton"
#' @param qual optional setting for which rule of thumb to use. Either "cohen" (default), "vd", "sawilowsky", "cohen-conv", "lovakov", "rosenthal", "brydges"
#' @param convert optional conversion to use (only for Glass version). Either "no", "cohen_d", "vda"
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' 
#' @details 
#' If a Cureton version of rank-biserial was used, the result for independent samples is the same as Goodman-Kruskal gamma, so we can also use those rules-of-thumb.
#' 
#' If a Glass version is used, then we can either use Cliff Delta, Somers d, or a conversion to either Cohen d, or Vargha-Delaney A.
#' 
#' See the separate functions on each of these for various rules-of-thumbs.
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to obtain the measure:
#' \code{\link{r_rank_biserial_is}}, to determine a the rank biserial for independent samples.
#' \code{\link{r_rank_biserial_os}}, to determine a the rank biserial for one-sample.
#' 
#' The function uses the convert function and corresponding rules of thumb:
#' \code{\link{es_convert}}, to convert this to Cohen d. 
#' 
#' @references 
#' Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). L. Erlbaum Associates.
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
th_rank_biserial <- function(rb, version="glass", qual=NULL, convert="no"){
  if (version == 'cureton'){
    if (is.null(qual)){
      qual = "blaikie"}
      results = th_gk_gamma(rb, qual)
  }
  else if (version=='glass'){
    if (convert == "cohen_d"){
      if (is.null(qual)){
        qual = "sawilowsky"}
      d = es_convert(rb, fr="rb", to="cohend")
      results = th_cohen_d(d, qual)
    }
    else if (convert == "vda" || convert == "cles"){
      if (is.null(qual)){
        qual = "vargha"}
      vda = es_convert(rb, fr="rb", to="cle")
      results = th_cle(vda, qual)
    }
    else{
      if (qual == "metsamuuronen-somers"){
        results = th_somers_d(rb, "metsamuuronen")}
      else{
        if (is.null(qual)){
          qual = "romano"}
        results = th_cliff_delta(rb, qual)}
    }
  }
  
  return(results)
  
}
