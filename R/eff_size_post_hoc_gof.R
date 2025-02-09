#' Effect Sizes for a Goodness-of-Fit Post-Hoc Analysis
#' 
#' @description 
#' Determines an effect size for each test (row) from the results of ph_pairwise_bin(), ph_pairwise_gof(), ph_residual_bin(), or ph_residual_gof().
#' 
#' @param post_hoc_results dataframe with the result of either ph_pairwise_bin(), ph_pairwise_gof(), ph_residual_bin(), or ph_residual_gof()
#' @param es string optional, the effect size to determine. Either 'auto', 'coheng', 'cohenh', 'ar', 'cramerv', 'cohenw', 'jbme', 'fei', 'rosenthal'
#' @param bergsma optional boolean. Use of Bergsma correction, only for Cramér V
#' 
#' @returns
#' a dataframe with for residual post-hoc:
#' \item{category}, the label of the category
#' \item{name effect size}, the effect size value
#' 
#' for pairwise post-hoc
#' \item{category 1}, the label of the first category
#' \item{category 2}, the label of the second category
#' \item{name effect size}, the effect size value
#' 
#' @details
#' 'auto' will use Cohen h for exact tests, Rosenthal correlation for z-tests and Cramér's V otherwise.
#' 
#' Cohen g ('coheng'), Cohen h ('cohenh') and Alternative Ratio ('ar') can all be used for any test.
#' 
#' Cramér V ('cramerv'), Cohen w ('cohenw'), Johnston-Berry-Mielke E ('jbme'), and Fei ('fei') can be used with chi-square tests (or likelihood ratio tests)
#' 
#' The Rosenthal Correlation ('rosenthal') can be used with a z-test (proportion/Wald/score/residual).
#' 
#' See the separate functions for each of these for details on the calculations.
#' 
#' @seealso 
#' The post-hoc tests:
#' \code{\link{ph_pairwise_gof}}
#' \code{\link{ph_pairwise_bin}}
#' \code{\link{ph_residual_gof_gof}}
#' \code{\link{ph_residual_gof_bin}}
#' 
#' The effect sizes:
#' \code{\link{es_cohen_g}}
#' \code{\link{es_cohen_h_os}}
#' \code{\link{es_alt_ratio}}
#' \code{\link{es_cramer_v_gof}}
#' \code{\link{es_cohen_w}}
#' \code{\link{es_jbm_e}}
#' \code{\link{es_fei}}
#' \code{\link{r_rosenthal}}
#' 
#' note: the effect size functions are not used themselves in this function, but the same formulas are used.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_post_hoc_gof <- function(post_hoc_results, es='auto', bergsma=FALSE){
  #rename the post-hoc results
  df = post_hoc_results
  
  #determine the number of tests in the post-hoc results
  n_tests = nrow(df)
  
  #get the description of the test used
  if ('test' %in% colnames(df)){
    test_used = df[1, 'test']
  }
  else{
    test_used = df[1, 'test used']
  }
  if (any(grepl(paste(c("binomial", "multinomial"), collapse = "|"), test_used))){
    ph_test = 'exact'}
  else if (any(grepl(paste(c('Wald', 'Score', 'adjusted', 'standardized'), collapse = "|"), test_used))){
    ph_test = 'z-test'}
  else if (any(grepl(paste(c('G test', 'likelihood'), collapse = "|"), test_used))){
    ph_test = 'likelihood-test'}
  else{
    ph_test = 'chi2-test'}
  
  #find the name of the test-statistic column
  if (ph_test!='exact'){
    if ('z-statistic' %in% colnames(df)){
      stat_col = 'z-statistic'}
    else{
      stat_col = 'statistic'}
  }
  
  #label for effect size
  es_labels <- c(
    coheng = "Cohen g",
    cohenh = "Cohen h",
    ar = "alternative ratio",
    cramerv = "Cramér V",
    cohenw = "Cohen w",
    jbme = "Johnston-Berry-Mielke E",
    fei = "Fei",
    rosenthal = "Rosenthal correlation"
  )
  
  #set the effect size measure if es='auto'
  if (es=='auto'){
    if (ph_test=='exact'){
      es = 'cohenh'}
    else if (ph_test=='z-test'){
      es = 'rosenthal'}
    else{
      es= 'cramerv'}
  }
  
  #determine if it was a pairwise or residual test
  if ('category 1'  %in% colnames(df)){
    test_type = 'pairwise'
    col_names = unname(c('category 1', 'category 2', es_labels[es]))
    res_col=3}
  else{
    test_type = 'residual'
    n = sum(as.numeric(df[['obs. count']]))
    col_names = unname(c('category', es_labels[es]))
    res_col=2}
  
  #loop over each row (test)
  results = data.frame()
  for (i in 1:n_tests){
    # find the observed and expected counts
    if (test_type == 'pairwise'){
      p_obs = as.numeric(df[i, 'obs. prop. 1'])
      p_exp = as.numeric(df[i, 'exp. prop. 1'])
      n_row = as.numeric(df[i, 'n1']) + as.numeric(df[i,'n2'])
      
      #add the two category names to the dataframe
      results[i, 1] = df[i, 1]
      results[i, 2] = df[i, 2]
    }
    else{
      n_row = as.numeric(n)
      p_obs = as.numeric(df[i, 'obs. count'])/n_row
      p_exp = as.numeric(df[i, 'exp. count'])/n_row
      
      #add the category name to the dataframe
      results[i, 1] = df[i, 1]
    }
    
    #for non-exact tests find the test-statistic value
    if (ph_test!='exact'){
      test_statistic = as.numeric(df[i, stat_col])
    }
    
    #determine the effect size
    #effect sizes for without a test-statistic needed
    if (es=="coheng"){
      es_value = p_obs - 0.5}
    else if (es=="cohenh"){
      es_value = 2*asin(p_obs**0.5) - 2*asin(p_exp**0.5)}
    else if( es=="ar"){
      es_value = p_obs/p_exp}
    
    #effect sizes with a z-statistic
    else if (es=="rosenthal"){
      if (ph_test == 'z-test'){
        es_value = as.numeric(df[i, stat_col])/(n_row**0.5)}
      else{
        es_value = 'not possible with this post-hoc test'}
    }
    
    #effect sizes with a chi-square statistic
    else if (es=='cramerv' || es=="cohenw"){
      if (ph_test == 'chi2-test' || ph_test == 'likelihood-test'){
        es_value = (test_statistic/n_row)**0.5
        if (es=='cramerv' && bergsma){
          phi2 = test_statistic/n_row
          phi2_tilde = max(0, phi2 - 1/(n_row-1))
          es_value = (phi2_tilde/(2 - 1/(n_row-1)))**0.5}
      }
      else{
        es_value = 'not possible with this post-hoc test'}
    }
    
    else if (es=="jbme" || es=="fei"){
      if (ph_test == 'chi2-test' || ph_test == 'likelihood-test'){
        if (ph_test == 'chi2-test'){
          es_value = test_statistic*as.numeric(df[i, 'minExp'])/(n_row*(n_row - as.numeric(df[i, 'minExp'])))
        }
        else{
          es_value = -1/log(as.numeric(df[i, 'minExp'])/n_row)*test_statistic/(2*n_row)}
        
        if (es=="fei"){
          es_value= es_value**0.5}
      }
      else{
        es_value = 'not possible with this post-hoc test'
      }
    }
    
    #add the effect size value to the dataframe
    results[i, res_col] = es_value 
  }
  
  #add the column names
  colnames(results) = col_names
  
  return (results)
}