#' Post-Hoc Goodness-of-Fit Rules-of-Thumb
#' 
#' @description 
#' This function will add a classification to the results of **es_post_hoc_gof()** using a rules-of-thumb. This is frowned upon by some, and the rule-of-thumb can vary per discipline.
#' 
#' @param eff_sizes dataframe, the dataframe from es_post_hoc_gof()
#' @param convert boolean, optional. convert the effect size to use the rule-of-thumb from another, see details
#' @param ph_results dataframe, optional. the post-hoc analysis results, required for JBM-E and Fei.
#' @param ... optional. additional arguments for the specific rule-of-thumb that are passed along. Most common 'qual=...' for a specific set of rules-of-thumb.
#' 
#' @returns
#' df, dataframe with the same dataframe as the provided *eff_sizes*, but added:
#' \item{qualification}, the qualification using the rule-of-thumb.
#' \item{reference}, a reference to the source for the rule-of-thumb.
#' 
#' If a conversion was done or needed:
#' \item{conversion description}, the value of the converted measure
#' 
#' @details
#' For Johnston-Berry-Mielke E and Fei, a conversion is always done to Cramér V, when setting *convert=True* it will convert it again to Cohen w.
#' 
#' Other possible conversions are Cohen h' to Cohen h, and Cramér V to Cohen w.
#' 
#' See the separate documentation for each of the rules-of-thumb, or conversion.
#' 
#' @seealso 
#' 
#' \code{\link{th_cohen_g}}
#' \code{\link{th_cohen_h}}
#' \code{\link{th_cohen_w}}
#' \code{\link{th_cramer_v}}
#' \code{\link{th_pearson_r}}
#' \code{\link{es_convert}}@author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
th_post_hoc_gof <- function(eff_sizes, convert=FALSE, ph_results=NULL, ...){
  df = eff_sizes
  qs = c()
  rs = c()
  conv = c()
  
  if ('Cohen h' %in% colnames(df) && convert){
    conv_lbl = 'Cohen h_2 to Cohen h'}
  else if ('Cohen w' %in% colnames(df) && convert){
    (conv_lbl = 'Cohen w to Cramér V')}
  else if ('Cramér V' %in% colnames(df) && convert){
    conv_lbl = 'Cramér V to Cohen w'}
  else if ('Johnston-Berry-Mielke E' %in% colnames(df)){
    if (convert){
      conv_lbl = 'JBM-E to Cramér V'}
    else{
      conv_lbl = 'JBM-E to Cohen w'}
  }
  else if ('Fei' %in% colnames(df)){
    if (convert){
      conv_lbl = 'Fei to Cramér V'}
    else{
      conv_lbl = 'Fei to Cohen w'}
  }
  else{
    conv_lbl="none"
  }
  
  if (df[1, length(colnames(df))] == 'not possible with this post-hoc test'){
    df = 'no effect size values found'}
  else{
    for (i in 1:nrow(df)){
      if ('Cohen g' %in% colnames(df)){
        q = th_cohen_g(df[i,'Cohen g'], ...)}                
      else if ('Cohen h' %in% colnames(df)){
        if (convert){
          cohenH = es_convert(df[i, 'Cohen h'], fr= "cohenhos", to = "cohenh")
          q = th_cohen_h(cohenH, ...)
          conv_val = cohenH}
        else{                    
          q = th_cohen_h(df[i, 'Cohen h'], ...)}
      }
      else if ('Cohen w' %in% colnames(df)){
        if (convert){
          cramerV = es_convert(df[i, 'Cohen w'], fr="cohenw", to="cramervgof", ex1=2)
          q = th_cramer_v(cramerV, ...)
          conv_val = cramerV}      
        else{
          q = th_cohen_w(df[i, 'Cohen w'], ...)}
      }
      else if ('Cramér V' %in% colnames(df)){
        if (convert){
          cohenW = es_convert(df[i, 'Cramér V'], fr="cramervgof", to="cohenw", ex1=2)
          q = th_cohen_w(cohenW, ...)
          conv_val = cohenW}      
        else{
          q = th_cramer_v(df[i, 'Cramér V'], ...)}
      }
      else if ('Johnston-Berry-Mielke E' %in% colnames(df)){
        if (is.null(ph_results)){
          q = 'need post-hoc resulst'}
        else{
          if ('n1' %in% colnames(ph_results)){
            n = ph_results[i, 'n1'] + ph_results[i, 'n2']}
          else{
            n = sum(ph_results[,'obs. count'])}
          cohenW = es_convert(df[i, 'Johnston-Berry-Mielke E'], fr="jbme", to="cohenw", ex1=ph_results[i,8]/n)
          if (convert){
            cramerV = es_convert(cohenW, fr="cohenw", to="cramervgof", ex1=2)
            q = th_cramer_v(cramerV, ...)
            conv_val = cramerV}   
          else{
            q = th_cohen_w(cohenW, ...)
            conv_val = cohenW}
        }
      }
      else if ('Fei' %in% colnames(df)){
        if (is.null(ph_results)){
          q = 'need post-hoc resulst'}
        else{
          if ('n1' %in% colnames(ph_results)){
            n = ph_results[i, 'n1'] + ph_results[i, 'n2']}
          else{
            n = sum(ph_results[,'obs. count'])}
          
          cohenW = es_convert(df[i, 'Fei'], fr="fei", to="cohenw", ex1=ph_results[i,8]/n)
          
          if (convert){
            cramerV = es_convert(cohenW, fr="cohenw", to="cramervgof", ex1=2)
            q = th_cramer_v(cramerV, ...)
            conv_val = cramerV}                    
          else{
            q = th_cohen_w(cohenW, ...)
            conv_val = cohenW}
        }
      }
      else if ('Rosenthal correlation' %in% colnames(df)){
        q = th_pearson_r(df[i, 'Rosenthal correlation'], ...)}
      
      else if ('alternative ratio' %in% colnames(df)){
        q = 'no rules-of-thumb available'}
      
      if (is.data.frame(q)){
        qs = c(qs, q[1,1])
        rs = c(rs, q[1,2])
        if (conv_lbl != "none"){
          conv = c(conv, conv_val)}
      }
      else{
        qs = c(qs, q)
        rs = c(rs, 'n.a.')
        conv = c(conv, 'n.a.')
      }
    }
    
    if (conv_lbl != "none"){
      df[conv_lbl] = conv}
    df['qualification'] = qs
    df['reference'] = rs
    
  }
  
  return (df)
  
}