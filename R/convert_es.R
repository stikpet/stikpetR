#' Convert Effect Size
#'
#' @param es the effect size value to convert
#' @param fr name of the original effect size (see details)
#' @param to name of the effect size to convert to (see details)
#' @param ex1 extra for some conversions (see details)
#' @param ex2 extra for some conversions (see details)
#' @return the converted effect size value
#'
#' @examples
#' es_convert(0.3, fr="cohenhos", to = "cohenh")
#'
#' @details
#'
#' **COHEN D**
#'
#' **Cohen d to Odds Ratio**
#'
#' fr="cohend", to="or", ex1="chinn"
#'
#' This uses (Chinn, 2000, p. 3129):
#'
#' \deqn{OR = e^{d\times 1.81}}
#'
#' fr="cohend", to="or", ex1="borenstein"
#'
#' This uses (Borenstein et. al, 2009, p. 3):
#'
#' \deqn{OR = e^{\frac{d\times\pi}{\sqrt{3}}}}
#' 
#' **Cohen d to Rank Biserial (Cliff delta)**
#' 
#' fr = "cohend", to = "rb"
#' 
#' This uses (Marfo & Okyere, 2019, p. 4):
#' \deqn{r_b = \frac{2\times\Phi\left(\frac{d}{2}\right)-1}{\Phi\left(\frac{d}{2}\right)}}
#' 
#' 
#' **COHEN D'**
#'
#' **Convert a Cohen d' to Cohen d**
#'
#' fr="cohendos" to="cohend"
#'
#' This uses (Cohen , 1988, p. 46):
#' \deqn{d = d'\times\sqrt{2}}
#'
#'
#' **COHEN F**
#'
#' **Cohen f to Eta-squared**
#'
#' fr="cohenf" to="etasq"
#'
#' This uses (Cohen, 1988, p. 284):
#' \deqn{\eta^2 = \frac{f^2}{1 + f^2}}
#'
#' **COHEN H'**
#'
#' **Cohen h' to Cohen h**
#'
#' fr = "cohenhos", to = "cohenh"
#'
#' This uses (Cohen, 1988, p. 203):
#' \deqn{h = h'\times\sqrt{2}}
#'
#' **COHEN w**
#'
#' **Cohen w to Contingency Coefficient**
#'
#' fr="cohenw", to="cc"
#' 
#' *Cohen w to Cramér V GoF*
#' 
#' fr="cohenw", to="cramervgof", ex1=k
#' 
#' This uses (Cohen, 1988, p. 223):
#' \deqn{v = \frac{w}{\sqrt{k - 1}}}
#' 
#' *Cohen w to Cramér V ind.*
#' 
#' fr="cohenw", to="cramervind", ex1=r, ex2=c
#' 
#' This uses:
#' \deqn{v = \frac{w}{\sqrt{\min\left(r - 1, c - 1\right)}}}
#' 
#' *Cohen w to Fei*
#' 
#' fr="cohenw", to="fei", ex1=minExp/n
#' 
#' This uses:
#' \deqn{Fei = \frac{w}{\sqrt{\frac{1}{p_E}-1}}}
#'
#'
#' **CRAMER V GoF**
#'
#' **Cramer's v for Goodness-of-Fit to Cohen w**
#'
#' fr="cramervgof", to = "cohenw", ex1 = k
#'
#' This uses (Cohen, 1988, p. 223):
#' \deqn{w = v\times\sqrt{df}}
#'
#' **EPSILON SQUARED**
#'
#' **Epsilon Squared to Eta Squared**
#'
#' fr="epsilonsq", to="etasq", ex1 = n, ex2 = k
#'
#' This uses:
#' \deqn{\eta^2 = 1 - \frac{\left(1 - \epsilon^2\right)\times\left(n - k\right)}{n - 1}}
#'
#' **Epsilon Squared to Omega Squared**
#'
#' fr="epsilonsq", to="omegasq", ex1=MS_w, ex2 = SS_t
#'
#' This uses:
#' \deqn{\hat{\omega}^2 = \epsilon^2\times\left(1 - \frac{MS_w}{SS_t + MS_w}\right)}
#'
#'
#' **ETA SQUARED**
#'
#' **Eta squared to Cohen f**
#' fr="etasq", to="cohenf")
#'
#' This uses:
#' \deqn{f = \sqrt{\frac{\eta^2}{1 - \eta^2}}}
#'
#' **Eta squared to Epsilon Squared**
#' fr="etasq", to="epsilonsq", ex1=n, ex2=k
#'
#' This uses:
#' \deqn{\epsilon^2 = \frac{n\times\eta^2 - k + \left(1 - \eta^2\right)}{n - l}}
#' 
#' **FEI**
#' 
#' *Fei to Cohen w*
#' 
#' fr="fei", to="cohenw", ex1=minExp/n
#' 
#' This uses:
#' \deqn{w = Fei\times\sqrt{\frac{1}{p_E}-1}}
#' 
#' *Fei to Johnston-Berry-Mielke E*
#' fr="fei", to="jbme"
#' 
#' This uses:
#' \deqn{E = Fei^2}
#'
#' **JOHNSTON-BERRY-MIELKE**
#'
#' **Johnston-Berry-Mielke E to Cohen w**
#'
#' fr="jbme", to="cohenw", ex1=minExp/n
#'
#' This uses (Johnston et al., 2006, p. 413):
#' \deqn{w = \sqrt{\frac{E\times\left(1 - \right)}{q}}}
#' 
#' *Johnston-Berry-Mielke E to Cohen w*
#' 
#' fr="jbme", to="fei"
#' 
#' This uses:
#' \deqn{Fei = \sqrt(E)}
#' 
#' **ODDS RATIO**
#'
#' **Odds Ratio to Cohen d**
#'
#' fr="or", to="cohend", ex1="chinn"
#'
#' This uses (Chinn, 2000, p. 3129):
#'
#' \deqn{d = \frac{\ln{\left(OR\right)}}{1.81}}
#'
#' fr="or", to="cohend", ex1="borenstein"
#'
#' This uses (Borenstein et. al, 2009, p. 3):
#'
#' \deqn{d = \ln\left(OR\right)\times\frac{\sqrt{3}}{\pi}}
#'
#' **Odds Ratio to Yule Q**
#'
#' fr="or", to="yuleq"
#'
#' This uses:
#'
#' \deqn{Q = \frac{OR - 1}{OR + 1}}
#'
#' **Odds Ratio to Yule Y**
#'
#' This uses
#' \deqn{Y = \frac{\sqrt{OR} - 1}{\sqrt{OR} + 1}}
#'
#' **OMEGA SQUARED**
#'
#' **Omega Squared to Epsilon Squared**
#'
#' fr="omegasq", to="epsilonsq", ex1=MS_w, ex2 = SS_t
#'
#' This uses:
#' \deqn{\epsilon^2 = \frac{\hat{\omega}^2}{1 - \frac{MS_w}{SS_t + MS_w}}}
#' 
#' **RANK BISERIAL (CLIFF DELTA)**
#' 
#' **Rank Biserial (Cliff delta) to Cohen d**
#' 
#' fr = "rb", to = "cohend"
#' 
#' This uses (Marfo & Okyere, 2019, p. 4):
#' \deqn{d =\sqrt{2} \times \Phi^{-1}\left(-\frac{1}{r_b-2}\right)}
#' 
#' **Rank Biserial (Cliff delta) to Vargha-Delaney A**
#' 
#' fr = "rb", to = "vda"
#' 
#' This uses:
#' \deqn{r_b = 2\times A - 1}
#' 
#' **VARGHA-DELANEY A**
#' 
#' **Vargha-Delaney A to Rank Biserial (Cliff delta)**
#' 
#' fr = "vda", to = "rb"
#' 
#' This uses:
#' \deqn{A = \frac{r_b + 1}{2}}
#' 
#' **YULE Q**
#'
#' **Yule Q to Odds Ratio**
#'
#' fr="yuleq", to="or"
#'
#' This uses:
#'
#' \deqn{OR = \frac{1 + Q}{1 - Q}}
#'
#' **Yule Q to Yule Y**
#'
#' fr="yuleq", to="yuley"
#'
#' This uses:
#' \deqn{Y = \frac{1 - sqrt{1 - Q^2}}{Q}}
#'
#' **YULE Y**
#'
#' **Yule Y to Yule Q**
#'
#' fr="yuley", to=="yuleq"
#'
#' This uses:
#' \deqn{Q = \frac{2\times Y}{1 + Y^2}}
#'
#' **Yule Y to Odds Ratio**
#'
#' fr="yuley", to=="or"
#'
#' This uses
#' \deqn{OR = \left(\frac{1 + Y}{1 - Y}\right)^2}
#'
#'
#'
#' @author
#' P. Stikker
#'
#' Please visit: https://PeterStatistics.com
#'
#' YouTube channel: https://www.youtube.com/stikpet
#'
#' @references
#' Borenstein, M., Hedges, L. V., Higgins, J. P. T., & Rothstein, H. R. (2009). Converting Among Effect Sizes. In *Introduction to Meta-Analysis*. John Wiley & Sons, Ltd. https://doi.org/10.1002/9780470743386
#'
#' Chinn, S. (2000). A simple method for converting an odds ratio to effect size for use in meta-analysis. *Statistics in Medicine, 19*(22), 3127–3131. https://doi.org/10.1002/1097-0258(20001130)19:22<3127::aid-sim784>3.0.co;2-m
#'
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#'
#' Johnston, J. E., Berry, K. J., & Mielke, P. W. (2006). Measures of effect size for chi-squared and likelihood-ratio goodness-of-fit tests. *Perceptual and Motor Skills, 103*(2), 412–414. https://doi.org/10.2466/pms.103.2.412-414
#'
#' @export
es_convert <- function(es, fr, to, ex1=NULL, ex2=NULL){

  #COHEN d
  #Cohen d one-sample to Cohen d
  if(fr=="cohendos" && to=="cohend") {res = es*sqrt(2)}

  #Cohen d to Odds Ratio
  else if(fr=="cohend" && to=="or") {
    #Chinn (2000, p. 3129)
    if(ex1=="chinn"){res = exp(1.81*es)}
    #Borenstein et. al (2009, p. 3)
    else{res = exp(es*pi/sqrt(3))}
  }
  
  #Cohen d to Rank Biserial (Cliff delta)
  else if (fr=="cohend" && to=="rb"){
    pf = pnorm(es/2)
    res = (2*pf - 1)/pf
  }
  
  #COHEN F
  #Cohen f to eta squared
  else if(fr=="cohenf" && to=="etasq") {res = es**2/(1 + es**2)}


  #COHEN h'
  #Cohen h' to Cohen h
  else if(fr=="cohenhos" && to=="cohenh") {res = es*sqrt(2)}

  #COHEN w
  #Cohen w to Contingency Coefficient
  else if(fr=="cohenw" && to=="cc"){res = sqrt(es^2 / (1 + es^2))}

  #Cohen w to Cramer V
  else if(fr=="cohenw" && to=="cramervgof"){res = es/sqrt(ex1 - 1)}
  else if(fr=="cohenw" && to=="cramervind"){res = es/sqrt(min(ex1 - 1, ex2 -1))}
  
  #Cohen w to Fei
  else if(fr=="cohenw" && to=="fei"){res = es/sqrt(1/ex1 - 1)}
  
  #CONTINGENCY COEFFICIENT
  #Contingency Coefficient to Cohen w
  else if(fr=="cc" && to=="cohenw"){res = sqrt(es^2 / (1 - es^2))}


  #CRAMÉR V
  #Cramer's v GoF to Cohen w
  else if(fr=="cramervgof" && to=="cohenw") {res = es*sqrt(ex1 - 1)}
  #Cramer's v Ind. to Cohen w
  else if(fr=="cramervind" && to=="cohenw") {res = es*sqrt(min(ex1 - 1, ex2 - 1))}

  #EPSILON SQUARED
  #Epsilon squared to Eta squared
  else if(fr=="epsilonsq" && to=="etasq") {res = 1 - (1 - es)*(ex1 - ex2)/(ex1 - 1)}
  #Epsilon squared to Omega squared
  else if(fr=="epsilonsq" && to=="omegasq") {res = es*(1 - ex1/(ex2 + ex1))}


  #ETA SQUARED
  #Eta squared to Cohen f
  else if(fr=="etasq" && to=="cohenf") {res = sqrt(es/(1-es))}

  #Eta squared to Epsilon Squared
  else if(fr=="etasq" && to=="epsilonsq") {res = (ex1*es - ex2 + (1 - es))/(ex1 - ex2)}
  
  #FEI    
  #Fei to Cohen w
  else if(fr=="fei" && to=="cohenw") {res = es*sqrt(1/ex1 - 1)}
  
  #Fei to Johnston-Berry-Mielke E
  else if(fr=="fei" && to=="jbme") {res = es**2}

  #JOHNSTON-BERRY-MIELKE E
  #Johnston-Berry-Mielke E to Cohen w
  else if(fr=="jbme" && to=="cohenw") {res = sqrt(es*(1 - ex1)/(ex1))}
  
  #Johnston-Berry-Mielke E to Fei
  else if(fr=="jbme" && to=="fei") {res = sqrt(es)}

  #ODDS RATIO
  #Odds Ratio to Cohen d (Chinn, 2000, p. 3129)
  else if(fr=="or" && to=="cohend") {
    if(ex1=="chinn"){res = log(es)/1.81}
    #Borenstein et. al (2009, p. 3)
    else{res = log(es)*sqrt(3)/pi}
  }
  #Odds Ratio to Yule Q
  else if (fr=="or" && to=="yuleq") {res = (es - 1)/(es + 1)}
  #Odds Ratio to Yule Y
  else if(fr=="or" && to=="yuley") {res = (sqrt(es) - 1)/(sqrt(es) + 1)}

  #OMEGA SQUARED
  #Omega squared to Epsilon squared
  else if(fr=="omegasq" && to=="epsilonsq") {res = es/(1 - ex1/(ex2 + ex1))}


  #RANK BISERIAL
  #Rank Biserial (Cliff delta) to Cohen d
  else if (fr=="rb" && to=="cohend"){res = 2*qnorm(-1/(es-2))}
  #Rank Biserial to Vargha and Delaney A
  else if(fr=="rb" && to=="cle") {res = (es + 1)/2}

  #VARGHA AND DELANEY A
  #Vargha and Delaney A to Rank Biserial
  else if(fr=="cle" && to=="rb") {res = 2*es - 1}


  #YULE Q
  #Yule Q to Odds Ratio
  else if(fr=="yuleq" && to=="or") {res = (1 + es)/(1 - es)}

  #Yule Q to Yule Y
  else if(fr=="yuleq" && to=="yuley") {res = (1 - sqrt(1 - es**2))/es}

  #YULE Y
  #Yule Y to Odds Ratio
  else if(fr=="yuley" && to=="or") {res = ((1 + es)/(1 - es))**2}

  #Yule Y to Yule Q
  else if(fr=="yuley" && to=="yuleq") {res = (2*es)/(1 + es**2)}

  return(res)

}
