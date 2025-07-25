#' Common Language (CL/CLES) (Independent Samples)
#' @description 
#' The Common Language Effect Size (a.k.a. Probability of Superiority) is the probability of taking a random pair from two categories, the first is greater than the first, i.e. 
#' \deqn{P(X > Y)}
#' 
#' Note however that Wolfe and Hogg (1971) actually had this in reverse, i.e.
#' \deqn{P(X \leq Y)}
#' 
#' Some will also argue to count ties equally to each of the two categories (Grissom, 1994, p. 282), which makes the definition:
#' \deqn{P(X > Y) + \frac{P(X = Y)}{2}}
#' 
#' It was further developed by Vargha and Delaney (2000) especially in light of a Mann-Whitney U test.
#' 
#' For scale data, an approximation using the standard normal distribution is also available.
#' 
#' The term Common Language Effect Size can be found in McGraw and Wong (1992), the term Probability of Superiority is found in Grissom (1994), and the term Stochastic Superiority in Vargha and Delaney (2000)
#' 
#' @param catField A vector with the categorical data
#' @param scores A vector with the scores
#' @param categories Optional to indicate which two categories of catField to use, otherwise first two found will be used.
#' @param levels Optional list with the ordinal text values in order
#' @param dmu Optional difference according to null hypothesis (default is 0)
#' @param method Optional optional method to use. "brute" will use a brute force " that will split ties evenly, "brute-it" is the same as brute but ignores ties, "vda" will use the calculation from Vargha-Delany, and "appr" a normal approximation from McGraw-Wong
#' 
#' @returns 
#' A dataframe with:
#' \item{CLE cat. 1}{the effect size for the first category}
#' \item{CLE cat. 2}{the effect size for  the second category}
#' 
#' @details
#' For "brute" simply all possible pairs are determined and half of the ties are added, i.e. (Grissom, 1994, p. 282):
#' \deqn{P(X > Y) + \frac{P(X = Y)}{2}}
#' 
#' With "brute-it" the ties are ignored (it = ignore ties):
#' \deqn{P(X > Y)}
#' 
#' The "appr" uses the approximation from McGraw and Wong (1992, p. 361):
#' \deqn{CL = \Phi\left(z\right)}
#' With:
#' \deqn{z = \frac{\left|\bar{x}_1 - \bar{x}_2\right|}{\sqrt{s_1^2 + s_2^2}}}
#' \deqn{s_i^2 = \frac{\sum_{j=1}^{n_i} \left(x_{i,j} - \bar{x}_i\right)^2}{n_i - 1}}
#' \deqn{\bar{x}_i = \frac{\sum_{j=1}^{n_i} x_{i,j}}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#'  \item  \eqn{x_{i,j}} the j-th score in category i
#' \item \eqn{n_i} the number of scores in category i
#' \item  \eqn{\Phi\left(\dots\right)} the cumulative density function of the standard normal distribution
#' }
#' 
#' The "vda" uses the formula used from Vargha and Delaney (2000, p. 107):
#' \deqn{A = \frac{1}{n_j}\times\left(\frac{R_i}{n_i} - \frac{n_i + 1}{2}\right)}
#' with \eqn{R_i} the sum of the ranks in category i
#' 
#' It could also be calculated from the Mann-Whitney U value:
#' \deqn{A = \frac{U}{n_1\times n_2}}
#' 
#' Note that the difference between the two options (using category 1 or category 2) will be the deviation from 0.5. If all scores in the first category are lower than the scores in the second, A will be 0 using the first category, and 1 for the second.
#' 
#' If the number of scores in the first category higher than the second, is the same as the other way around, A (no matter which category used) will be 0.5.
#' 
#' The CLE can be converted to a Rank Biserial (= Cliff delta) using the **es_convert()** function. This can then be converted to a Cohen d, and then the rules-of-thumb for Cohen d could be used (**th_cohen_d()**)
#' 
#' The CLE for the other category is simply 1 - CLE, except for the case where ties are ignored ("brute-it").
#' 
#' @seealso 
#' \code{\link{th_cle}}, to find rules-of-thumb for the CLE
#' 
#' @references 
#' Grissom, R. J. (1994). Statistical analysis of ordinal categorical status after therapies. *Journal of Consulting and Clinical Psychology, 62*(2), 281-284. doi:10.1037/0022-006X.62.2.281
#' 
#' McGraw, K. O., & Wong, S. P. (1992). A common language effect size statistic. *Psychological Bulletin, 111*(2), 361-365. doi:10.1037/0033-2909.111.2.361
#' 
#' Vargha, A., & Delaney, H. D. (2000). A critique and improvement of the CL common language effect size statistics of McGraw and Wong. *Journal of Educational and Behavioral Statistics, 25*(2), 101-132. doi:10.3102/10769986025002101
#' 
#' Wolfe, D. A., & Hogg, R. V. (1971). On constructing statistics and reporting data. *The American Statistician, 25*(4), 27-30. doi:10.1080/00031305.1971.10477278
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_common_language_is <- function(catField, scores, categories=NULL, levels=NULL, dmu=0, method="brute"){
  
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
  
  n1 = length(x1)
  n2 = length(x2)
  
  if (method=="appr"){
    var1 = var(x1)
    var2 = var(x2)
    
    m1 = mean(x1)
    m2 = mean(x2)
    
    z = (m1 - m2 - dmu)/sqrt(var1 + var2)
    
    c1 = pnorm(z)
    
    c2 = 1 - c1        
  }  
  else if (method=="vda"){
    n = n1 + n2
    full_list = c(x1, x2)
    all_ranks <- rank(full_list)
    cat1Ranks <- all_ranks[1:n1]
    cat2Ranks <- all_ranks[(n1 + 1):n]
    R1 = sum(cat1Ranks)
    R2 = sum(cat2Ranks)
    c1 = 1/n2*(R1/n1 - (n1 + 1)/2)
    c2 = 1/n1*(R2/n2 - (n2 + 1)/2)
    
  }
  
  
  else if (method=="brute" || method=="brute-it"){
    pairs <- expand.grid(x1, x2)
    #total number of pairs
    n = nrow(pairs)
    
    #pairs where first is greater than second
    xGTy <- subset(pairs, var1 > var2)
    
    if (method=="brute"){
      xEQy <- subset(pairs, var1 == var2)
      n_xGTEy = nrow(xGTy) + 1/2*nrow(xEQy)
      c1 = n_xGTEy/n
      c2 = 1 - c1}
    else {
      c1 = nrow(xGTy)/n
      c2 = nrow(subset(pairs, var1 < var2))/n
    }
    
    
  }
  
  results <- data.frame(c1, c2)
  colnames(results) = c(paste("CLE", cat1), paste("CLE", cat2))
  
  return(results)
}



