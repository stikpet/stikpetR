#' Rules of Thumb for Pearson Correlation Coefficient
#' 
#' @description 
#' 
#' This function will give a qualification (classification) to a given correlation coefficient
#' 
#' @param r the correlation coefficient
#' @param qual optional the rule of thumb to be used. Either "`bartz"` (default), `"rafter"`, `"cohen"`, `"rumsey"`, `"rosenthal"`, `"agnes"`, `"disha"`, `"hopkins"`, `"funder"`, `"gignac"`, `"hemphill"`, or `"lovakov"`
#' 
#' @returns dataframe with the qualification and source
#' 
#' @details 
#' 
#' The following rules-of-thumb can be used:
#' 
#' "rafter" => Rafter et al. (2003, p. 194)
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.25 | weak |
#' |0.25 < 0.75 | moderate |
#' |0.75 or more | strong |
#' 
#' "cohen" => Cohen (1988, p. 82)
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.20 | negligible |
#' |0.20 < 0.50 | small |
#' |0.50 < 0.80 | medium |
#' |0.80 or more | large |
#' 
#' "rumsey" => Rumsey (2011, p. 284)
#' 
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.30 | negligible |
#' |0.30 < 0.50 | weak |
#' |0.50 < 0.70 | moderate |
#' |0.70 or more | strong |
#' 
#' "gignac" => Gignac and Szodorai (2016, p. 75)
#' 
#' "hemphill" => Hemphill (2003, p. 78)
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.20 | small |
#' |0.20 < 0.30 | medium |
#' |0.30 or more | large |
#' 
#' "lovakov" => Lovakov and Agadullina (2021, p. 514)
#' 
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.12 | negligible |
#' |0.12 < 0.24 | small |
#' |0.24 < 0.41 | medium |
#' |0.41 or more | large |
#' 
#' "rosenthal" => Rosenthal (1996, p. 45)
#' 
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | negligible |
#' |0.10 < 0.30 | small |
#' |0.30 < 0.50 | medium |
#' |0.50 < 0.70 | large |
#' |0.70 or more | very large |
#' 
#' "agnes" => Agnes (2011)
#' 
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.20 | negligible |
#' |0.20 < 0.40 | low |
#' |0.40 < 0.60 | moderate |
#' |0.60 < 0.80 | marked |
#' |0.80 or more | high |
#' 
#' "bartz" => Bartz (1999, p. 184, as cited in Warmbrod 2001)
#' 
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.20 | very low |
#' |0.20 < 0.40 | low |
#' |0.40 < 0.60 | moderate |
#' |0.60 < 0.80 | strong |
#' |0.80 or more | very high |
#' 
#' "disha" => Disha (2016)
#' 
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | markedly low and negligible |
#' |0.10 < 0.30 | very low |
#' |0.30 < 0.50 | low |
#' |0.50 < 0.70 | moderate |
#' |0.70 < 0.90 | high |
#' |0.90 or more | very high |
#' 
#' "hopkins" => Hopkins (1997, as cited in Warmbrod 2001)
#' 
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.10 | trivial |
#' |0.10 < 0.30 | low |
#' |0.30 < 0.50 | moderate |
#' |0.50 < 0.70 | high |
#' |0.70 < 0.90 | very large |
#' |0.90 or more | nearly perfect |
#' 
#' "funder" => Funder and Ozer (2019, p. 166)
#' 
#' |\|r\|| Interpretation|
#' |---|----------|
#' |0.00 < 0.05 | negligible |
#' |0.05 < 0.10 | very small |
#' |0.10 < 0.20 | small |
#' |0.20 < 0.30 | medium |
#' |0.30 < 0.40 | large |
#' |0.40 or more | very large |
#' 
#' @references
#' Agnes. (2011, April 16). Correlation – Correlation coefficient, r. Finance Training Course. https://financetrainingcourse.com/education/2011/04/correlation-correlation-coefficient-r/
#' 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' Disha, M. (2016, November 3). Correlation: Meaning, types and its computation. Your Article Library. https://www.yourarticlelibrary.com/statistics-2/correlation-meaning-types-and-its-computation-statistics/92001
#' 
#' Funder, D. C., & Ozer, D. J. (2019). Evaluating effect size in psychological research: Sense and nonsense. *Advances in Methods and Practices in Psychological Science, 2*(2), 156–168. https://doi.org/10.1177/2515245919847202
#' 
#' Gignac, G. E., & Szodorai, E. T. (2016). Effect size guidelines for individual differences researchers. *Personality and Individual Differences, 102*, 74–78. https://doi.org/10.1016/j.paid.2016.06.069
#' 
#' Hemphill, J. F. (2003). Interpreting the magnitudes of correlation coefficients. *American Psychologist, 58*(1), 78–79. https://doi.org/10.1037/0003-066X.58.1.78
#' 
#' Lovakov, A., & Agadullina, E. R. (2021). Empirically derived guidelines for effect size interpretation in social psychology. *European Journal of Social Psychology, 51*(3), 485–504. https://doi.org/10.1002/ejsp.2752
#' 
#' Rafter, J. A., Abell, M. L., & Braselton, J. P. (2003). *Statistics with Maple*. Academic Press.
#' 
#' Rosenthal, J. A. (1996). Qualitative descriptors of strength of association and effect size. *Journal of Social Service Research, 21*(4), 37–59. https://doi.org/10.1300/J079v21n04_02
#' 
#' Rumsey, D. J. (2011). *Statistics for dummies* (2nd ed.). Wiley.
#' 
#' Warmbrod, J. R. (2001). Calculating, interpreting, and reporting estimates of “effect size": Magnitude of an effect or the strength of a relationship. https://www.depts.ttu.edu/aged/toolbox/effect_size.pdf
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet)
#' 
#' @export
th_pearson_r <- function(r, qual="bartz"){
  #Rafter et al. (2003, p. 194).
  if (qual=="rafter"){
    ref = "Rafter et al. (2003, p. 194)"
    if (abs(r) < 0.25){qual = "weak"}
    else if (abs(r) < 0.75){qual = "moderate"}
    else{qual = "strong"}
  }
  
  #Cohen (1988, p. 82).
  else if (qual=="cohen"){
    ref = "Cohen (1988, p. 82)"
    if (abs(r) < 0.1){qual = "negligable"}
    else if (abs(r) < 0.3){qual = "small"}
    else if (abs(r) < 0.5){qual = "medium"}
    else{qual = "large"}
  }
  #Rumsey (2011, p. 284).
  else if (qual=="rumsey"){
    ref = "Rumsey (2011, p. 284)"
    if (abs(r) < 0.3){qual = "negligable"}
    else if (abs(r) < 0.5){qual = "weak"}
    else if (abs(r) < 0.7){qual = "moderate"}
    else{qual = "strong"}
  }
  #Gignac and Szodorai (2016, p. 75); Hemphill (2003, p. 78).
  else if (qual=="gignac" || qual=="hemphill"){
    ref = "Gignac and Szodorai (2016, p. 75); Hemphill (2003, p. 78)"
    if (abs(r) < 0.1){qual = "negligable"}
    else if (abs(r) < 0.2){qual = "small"}
    else if (abs(r) < 0.3){qual = "medium"}
    else{qual = "large"}
  }
  #Lovakov and Agadullina (2021, p. 514).
  else if (qual=="lovakov"){
    ref = "Lovakov and Agadullina (2021, p. 514)"
    if (abs(r) < 0.12){qual = "negligable"}
    else if (abs(r) < 0.24){qual = "small"}
    else if (abs(r) < 0.41){qual = "medium"}
    else{qual = "large"}
  }
  #Rosenthal (1996, p. 45).
  else if (qual=="rosenthal"){
    ref = "Rosenthal (1996, p. 45)"
    if (abs(r) < 0.1){qual = "negligable"}
    else if (abs(r) < 0.3){qual = "small"}
    else if (abs(r) < 0.5){qual = "medium"}
    else if (abs(r) < 0.7){qual = "large"}
    else{qual = "very large"}
  }
  #Agnes (2011).
  else if (qual=="agnes"){
    ref = "Agnes (2011)"
    if (abs(r) < 0.2){qual = "negligable"}
    else if (abs(r) < 0.4){qual = "low"}
    else if (abs(r) < 0.6){qual = "moderate"}
    else if (abs(r) < 0.8){qual = "marked"}
    else{qual = "high"}
  }
  
  # Bartz (1999, p. 184, as cited in Warmbrod 2001).
  else if (qual=="bartz"){
    ref = "Bartz (1999, p. 184, as cited in Warmbrod 2001)"
    if (abs(r) < 0.2){qual = "very low"}
    else if (abs(r) < 0.4){qual = "low"}
    else if (abs(r) < 0.6){qual = "moderate"}
    else if (abs(r) < 0.8){qual = "strong"}
    else{qual = "very high"}
  }
  #Disha (2016).
  else if (qual=="disha"){
    ref = "Disha (2016)"
    if (abs(r) < 0.1){qual = "markedly low and negligible"}
    else if (abs(r) < 0.3){qual = "very low"}
    else if (abs(r) < 0.5){qual = "low"}
    else if (abs(r) < 0.7){qual = "moderate"}
    else if (abs(r) < 0.9){qual = "high"}
    else{qual = "very high"}
  }
  
  #Hopkins (1997, as cited in Warmbrod 2001).
  else if (qual=="hopkins"){
    ref = "Hopkins (1997, as cited in Warmbrod 2001)"
    if (abs(r) < 0.1){qual = "trivial"}
    else if (abs(r) < 0.3){qual = "low"}
    else if (abs(r) < 0.5){qual = "moderate"}
    else if (abs(r) < 0.7){qual = "high"}
    else if (abs(r) < 0.9){qual = "very large"}
    else{qual = "nearly perfect"}
  }
  
  #Funder and Ozer (2019, p. 166).
  else if (qual=="funder"){
    ref = "Funder and Ozer (2019, p. 166)"
    if (abs(r) < 0.05){qual = "negligable"}
    else if (abs(r) < 0.1){qual = "very small"}
    else if (abs(r) < 0.2){qual = "small"}
    else if (abs(r) < 0.3){qual = "medium"}
    else if (abs(r) < 0.4){qual = "large"}
    else{qual = "very large"}
  }
  
  results = data.frame(qual, ref)
  colnames(results)<-c("classification", "reference")
  
  return(results)
}