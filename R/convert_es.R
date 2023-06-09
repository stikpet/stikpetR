#' Convert Effect Size
#'
#' @param es the effect size value to convert
#' @param from name of the original effect size (see details)
#' @param to name of the effect size to convert to (see details)
#' @param ex1 extra for some conversions (see details)
#' @param ex2 extra for some conversions (see details)
#' @return the converted effect size value
#'
#' @examples
#' es_convert(0.3, from="cohenhos", to = "cohenh")
#'
#' @details
#'
#' **COHEN D**
#'
#' **Cohen d to Odds Ratio**
#'
#' from="cohend", to="or", ex1="chinn"
#'
#' This uses (Chinn, 2000, p. 3129):
#'
#' \deqn{OR = e^{d\times 1.81}}
#'
#' from="cohend", to="or", ex1="borenstein"
#'
#' This uses (Borenstein et. al, 2009, p. 3):
#'
#' \deqn{OR = e^{\frac{d\times\pi}{\sqrt{3}}}}
#'

#'
#' **COHEN D'**
#'
#' **Convert a Cohen d' to Cohen d**
#'
#' from="cohendos" to="cohend"
#'
#' This uses (Cohen , 1988, p. 46):
#' \deqn{d = d'\times\sqrt{2}}
#'
#'
#' **COHEN F**
#'
#' **Cohen f to Eta-squared**
#'
#' from="cohenf" to="etasq"
#'
#' This uses (Cohen, 1988, p. 284):
#' \deqn{\eta^2 = \frac{f^2}{1 + f^2}}
#'
#' **COHEN H'**
#'
#' **Cohen h' to Cohen h**
#'
#' from = "cohenhos", to = "cohenh"
#'
#' This uses (Cohen, 1988, p. 203):
#' \deqn{h = h'\times\sqrt{2}}
#'
#' **COHEN w**
#'
#' **Cohen w to Contingency Coefficient**
#'
#' from="cohenw", to="cc"
#'
#'
#' **CRAMER V GoF**
#'
#' **Cramer's v for Goodness-of-Fit to Cohen w**
#'
#' from="cramervgof", to = "cohenw", ex1 = k
#'
#' This uses (Cohen, 1988, p. 223):
#' \deqn{w = v\times\sqrt{df - 1}}
#'
#' **EPSILON SQUARED**
#'
#' **Epsilon Squared to Eta Squared**
#'
#' from="epsilonsq", to="etasq", ex1 = n, ex2 = k
#'
#' This uses:
#' \deqn{\eta^2 = 1 - \frac{\left(1 - \epsilon^2\right)\times\left(n - k\right)}{n - 1}}
#'
#' **Epsilon Squared to Omega Squared**
#'
#' from="epsilonsq", to="omegasq", ex1=MS_w, ex2 = SS_t
#'
#' This uses:
#' \deqn{\hat{\omega}^2 = \epsilon^2\times\left(1 - \frac{MS_w}{SS_t + MS_w}\right)}
#'
#'
#' **ETA SQUARED**
#'
#' **Eta squared to Cohen f**
#' from="etasq", to="cohenf")
#'
#' This uses:
#' \deqn{f = \sqrt{\frac{\eta^2}{1 - \eta^2}}}
#'
#' **Eta squared to Epsilon Squared**
#' from="etasq", to="epsilonsq", ex1=n, ex2=k
#'
#' This uses:
#' \deqn{\epsilon^2 = \frac{n\times\eta^2 - k + \left(1 - \eta^2\right)}{n - l}}
#'
#' **JOHNSTON-BERRY-MIELKE**
#'
#' **Johnston-Berry-Mielke E to Cohen w**
#'
#' from="jbme", to="cohenw", ex1=minExp/n
#'
#' This uses (Johnston et al., 2006, p. 413):
#' \deqn{w = \sqrt{\frac{E\times\left(1 - \right)}{q}}}
#'
#' **ODDS RATIO**
#'
#' **Odds Ratio to Cohen d**
#'
#' from="or", to="cohend", ex1="chinn"
#'
#' This uses (Chinn, 2000, p. 3129):
#'
#' \deqn{d = \frac{\ln{\left(OR\right)}}{1.81}}
#'
#' from="or", to="cohend", ex1="borenstein"
#'
#' This uses (Borenstein et. al, 2009, p. 3):
#'
#' \deqn{d = \ln\left(OR\right)\times\frac{\sqrt{3}}{\pi}}
#'
#' **Odds Ratio to Yule Q**
#'
#' from="or", to="yuleq"
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
#' from="omegasq", to="epsilonsq", ex1=MS_w, ex2 = SS_t
#'
#' This uses:
#' \deqn{\epsilon^2 = \frac{\hat{\omega}^2}{1 - \frac{MS_w}{SS_t + MS_w}}}
#'
#'
#' **YULE Q**
#'
#' **Yule Q to Odds Ratio**
#'
#' from="yuleq", to="or"
#'
#' This uses:
#'
#' \deqn{OR = \frac{1 + Q}{1 - Q}}
#'
#' **Yule Q to Yule Y**
#'
#' from="yuleq", to="yuley"
#'
#' This uses:
#' \deqn{Y = \frac{1 - sqrt{1 - Q^2}}{Q}}
#'
#' **YULE Y**
#'
#' **Yule Y to Yule Q**
#'
#' from="yuley", to=="yuleq"
#'
#' This uses:
#' \deqn{Q = \frac{2\times Y}{1 + Y^2}}
#'
#' **Yule Y to Odds Ratio**
#'
#' from="yuley", to=="or"
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
es_convert <- function(es, from, to, ex1=NULL, ex2=NULL){
  #"or", "cohend" , "yuleq", "yuley", "vda", "rb"
  #from="cohenh2", to = "cohenh"
  #from="cramervgof", to = "cohenw", ex1 = df


  #COHEN d
  #Cohen d one-sample to Cohen d
  if(from=="cohendos" && to=="cohend") {
    res = es*sqrt(2)
  }

  #Cohen d to Odds Ratio
  else if(from=="cohend" && to=="or") {
    #Chinn (2000, p. 3129)
    if(ex1=="chinn"){
      res = exp(1.81*es)
    }
    else{
      #Borenstein et. al (2009, p. 3)
      res = exp(es*pi/sqrt(3))
    }

  }

  #COHEN F
  #Cohen f to eta squared
  else if(from=="cohenf" && to=="etasq") {
    res = es**2/(1 + es**2)
  }


  #COHEN h'
  #Cohen h' to Cohen h
  else if(from=="cohenhos" && to=="cohenh") {
    res = es*sqrt(2)
  }

  #COHEN w
  #Cohen w to Contingency Coefficient
  else if(from=="cohenw" && to=="cc"){
    res = sqrt(es^2 / (1 + es^2))
  }


  #CONTINGENCY COEFFICIENT
  #Contingency Coefficient to Cohen w
  else if(from=="cc" && to=="cohenw"){
    res = sqrt(es^2 / (1 - es^2))
  }


  #CRAMÉR V
  #Cramer's v GoF to Cohen w
  else if(from=="cramervgof" && to=="cohenw") {
    res = es*sqrt(ex1 - 1)
  }

  #EPSILON SQUARED
  #Epsilon squared to Eta squared
  else if(from=="epsilonsq" && to=="etasq") {
    res = 1 - (1 - es)*(ex1 - ex2)/(ex1 - 1)
  }
  #Epsilon squared to Omega squared
  else if(from=="epsilonsq" && to=="omegasq") {
    res = es*(1 - ex1/(ex2 + ex1))
  }




  #ETA SQUARED
  #Eta squared to Cohen f
  else if(from=="etasq" && to=="cohenf") {
    res = sqrt(es/(1-es))
  }

  #Eta squared to Epsilon Squared
  else if(from=="etasq" && to=="epsilonsq") {
    res = (ex1*es - ex2 + (1 - es))/(ex1 - ex2)
  }

  #JOHNSTON-BERRY-MIELKE E
  #Johnston-Berry-Mielke E to Cohen w
  else if(from=="jbme" && to=="cohenw") {
    res = sqrt(es*(1 - ex1)/(ex1))
  }

  #ODDS RATIO
  #Odds Ratio to Cohen d (Chinn, 2000, p. 3129)
  else if(from=="or" && to=="cohend") {
    if(ex1=="chinn"){
      res = log(es)/1.81
    }
    else{
      #Borenstein et. al (2009, p. 3)
      res = log(es)*sqrt(3)/pi
    }
  }
  #Odds Ratio to Yule Q
  else if (from=="or" && to=="yuleq") {
    res = (es - 1)/(es + 1)
  }
  #Odds Ratio to Yule Y
  else if(from=="or" && to=="yuley") {
    res = (sqrt(es) - 1)/(sqrt(es) + 1)
  }

  #OMEGA SQUARED
  #Omega squared to Epsilon squared
  else if(from=="omegasq" && to=="epsilonsq") {
    res = es/(1 - ex1/(ex2 + ex1))
  }


  #RANK BISERIAL
  #Rank Biserial to Vargha and Delaney A
  else if(from=="rb" && to=="vda") {
    res = (es + 1)/2
  }



  #VARGHA AND DELANEY A
  #Vargha and Delaney A to Rank Biserial
  else if(from=="vda" && to=="rb") {
    res = 2*es - 1
  }



  #YULE Q
  #Yule Q to Odds Ratio
  else if(from=="yuleq" && to=="or") {
    res = (1 + es)/(1 - es)
  }

  #Yule Q to Yule Y
  else if(from=="yuleq" && to=="yuley") {
    res = (1 - sqrt(1 - es**2))/es
  }

  #YULE Y
  #Yule Y to Odds Ratio
  else if(from=="yuley" && to=="or") {
    res = ((1 + es)/(1 - es))**2
  }

  #Yule Y to Yule Q
  else if(from=="yuley" && to=="yuleq") {
    res = (2*es)/(1 + es**2)
  }






  return(res)

}
