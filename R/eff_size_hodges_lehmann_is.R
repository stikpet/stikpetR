#' Hodges-Lehmann Estimator (Independent Samples)
#' @description 
#' The Hodges-Lehmann estimate, is the median of all the possible differences between two sets of data. The authors (Hodges & Lehmann, 1984) describe it as the location shift that is needed to align two distributions (with similar distributions) as much as possible (p. 599).
#' 
#' It is sometimes incorrectly described as the difference between the two medians, but that is incorrect. It is not uncommon to have a different Hodges-Lehmann estimate than simply taking the difference between the two medians.
#' 
#' This measure is sometimes mentioned as an effect size measure for a Mann-Whitney U / Wilcoxon Rank Sum test (van Geloven, 2018), however since it is a median of the possible differences, it is not standardized (i.e. it doesn't range between two fixed values, and depends therefor on the data).
#' 
#' @param catField A vector with the categorical data
#' @param scores A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' @param levels Optional list with the ordinal text values in order
#' 
#' @returns 
#' HL, the Hodges-Lehmann Estimator
#' \item{CLE cat. 1}{the effect size for the first category}
#' \item{CLE cat. 2}{the effect size for  the second category}
#' 
#' @details
#' The formula for the Hodges-Lehmann estimator with two samples is (Hodges & Lehmann, 1963, p. 602):
#' 
#' \deqn{HL = \text{median}\left(x_i - y_j | 1 \leq i \leq n_x, 1 \leq j \leq n_y\right)}
#' 
#' *Symbols used:*
#' \itemize{
#'  \item  \eqn{x_i} the i-th score in category x
#' \item \eqn{x_j} the j-th score in category y
#' \item  \eqn{n_i} the number of scores in category i
#' }
#' 
#' There might be a faster method to actually determine this. Algorithm 616 (Monahan, 1984), but couldn't translate the Fortran to R
#' 
#' @references 
#' Hodges, J. L., & Lehmann, E. L. (1963). Estimates of location based on rank tests. *The Annals of Mathematical Statistics, 34*(2), 598–611. doi:10.1214/aoms/1177704172
#' 
#' Monahan, J. F. (1984). Algorithm 616: Fast computation of the Hodges-Lehmann location estimator. *ACM Transactions on Mathematical Software, 10*(3), 265–270. doi:10.1145/1271.319414
#' 
#' van Geloven, N. (2018, March 13). Mann-Whitney U toets [Wiki]. Wikistatistiek. https://wikistatistiek.amc.nl/Mann-Whitney_U_toets
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_hodges_lehmann_is <- function(catField, scores, categories=NULL, levels=NULL){
  
  #remove rows with missing values
  df = data.frame(scores, catField)
  df = na.omit(df)
  colnames(df) = c("score", "group")
  
  #replace the ordinal values if levels is provided
  if (!is.null(levels)){
    df$score = factor(df$score, ordered = TRUE, levels = levels)        
  }
  df$score = as.numeric(df$score)
  
  #the two categories
  if (!is.null(categories)){
    cat1 = categories[1]
    cat2 = categories[2]
  }
  else {
    cat1 = names(table(df$group))[1]
    cat2 = names(table(df$group))[2]
  }
  
  x1 = df$score[df$group == cat1]
  x2 = df$score[df$group == cat2]
  
  pairs <- expand.grid(x=x1, y=x2)
  difs <- pairs$x - pairs$y
  
  hl = median(difs)
  
  return(hl)
}