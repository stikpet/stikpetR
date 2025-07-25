#' Rule-of-Thumb for Common Language Effect Size
#' 
#' @description
#' This function will give a qualification (classification) for a Common Language Effect Size (/ Vargha-Delaney A / Probability of Superiority)
#' 
#' The measure is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/EffectSizes/CommonLanguageEffectSize.html)
#' 
#' @param cle the Vargha-Delaney A value
#' @param qual c("vd", others via convert), optional rules-of-thumb to use, currently only 'vd' for Vargha-Delaney, otherwise a converted measure
#' @param convert c("no", "rb", "cohen_d"), optional list in case to use a rule-of-thumb from a converted measure. Either "no" for no conversion, "rb" for rank-biserial, or "cohen_d" for Cohen d.
#' 
#' @returns 
#' A dataframe with:
#' \item{classification}{the qualification of the effect size}
#' \item{reference}{a reference for the rule of thumb used}
#' 
#' @details 
#' Vargha and Delaney (2000, p. 106):
#' 
#' |\|0.5 - A\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.06 | negligible |
#' |0.11 < 0.14 | small |
#' |0.28 < 0.21 | medium |
#' |0.21 or more | large |
#' 
#' The CLE can be converted to a Rank Biserial Coefficient using:
#' \deqn{r_b = 2\times CLE - 1}
#' 
#' Rules of thumb from the **th_rank_biserial()** function could then be used, by setting: *convert="rb"*, and *qual* is any of the options in th_rank_biserial()
#' 
#' This in turn can be converted to Cohen's d using (Marfo & Okyere, 2019, p.4):
#' \deqn{d = 2\times \phi^{-1}\left(-\frac{1}{r_b - 2}\right)}
#' Rules of thumb from the **th_cohen_d()** function could then be used, by setting: *convert="cohen_d"*, and *qual* is any of the options in th_cohen_d()
#' 
#' @section Before, After and Alternatives:
#' Before this you might want to obtain the measure:
#' \code{\link{es_common_language_os}}, o determine the CLE for one-sample.
#' \code{\link{es_common_language_is}}, o determine the CLE for independent samples.
#' 
#' The function uses the convert function and corresponding rules of thumb: 
#' \code{\link{es_convert}}, for the conversions.
#' \code{\link{th_rank_biserial}}, for options for rules of thumb when converting to Rank Biserial.
#' \code{\link{th_cohen_d}}, for options for rules of thumb when converting to Cohen d.
#' 
#' 
#' @references 
#' Marfo, P., & Okyere, G. A. (2019). The accuracy of effect-size estimates under normals and contaminated normals in meta-analysis. *Heliyon, 5*(6), e01838. doi:10.1016/j.heliyon.2019.e01838
#' 
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. *Journal of Educational and Behavioral Statistics, 25*(2), 101-132. doi:10.3102/10769986025002101
#' 
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @examples 
#' # Example 1: Using Vargha and Delaney rules:
#' cle = 0.23
#' th_cle(cle)
#' 
#' # Example 2: Convert to rank-biserial and use Sawilowsky rules:
#' cle = 0.23
#' th_cle(cle, qual="sawilowsky", convert="rb")
#' 
#' @export
th_cle <- function(cle, qual="vd", convert="no"){
  
  if (convert=="no"){
    es = abs(0.5 - cle)
  
    if (qual=="vd"){
      #Use Vargha and Delaney (2000, p. 106)
      ref = "Vargha and Delaney (2000, p. 106)"
      if (es < 0.06){
        qual = "negligible"}
      else if (es < 0.14){
        qual = "small"}
      else if (es < 0.21){
        qual = "medium"}
      else{
        qual = "large"}
      
      results = data.frame(qual, ref)
      colnames(results)<-c("classification", "reference")
      
    }
  }
  else if (convert=="rb"){
    rb = es_convert(cle, fr="cle", to="rb")
    results = th_rank_biserial(rb, qual=qual)}
  else if (convert=="cohen_d"){
    rb = es_convert(cle, fr="cle", to="rb")
    d = es_convert(rb ,fr="rb", to="cohend")
    results = th_cohen_d(d, qual=qual)}
    
  return(results)
  
}



