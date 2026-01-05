#' Brunner-Munzel Test
#' @description
#' The Brunner-Munzel test, tests so-called stochastic dominance. It tests if you pick a random score from the first category, will there be a higher (or lower) chance than 50% it will be higher than a random score from the other category. Note this is not the same as a median test. A very short example.
#' 
#' If we have the scores 40, 50 and 60 in group A, and the scores 30, 50, and 51 in group B, the median of each group is 50. However if we pick 40 from group A, there is a 1/3 chance it is higher than a value from group B. If we pick 50 there is also a 1/3 chance, and if we pick 60 there is a 3/3 chance. Together, there is a (1/3+1/3+3/3)/3 = 5/9 chance that if you pick a random value from group A and B, that the one in group A is higher. This is quite off from the 50%.
#' 
#' The test could therefor be used if we have a binary and an ordinal variable. The Mann-Whitney U test is probably the most famous one used in these scenarios, but actually has quite some limitations, and this Brunner-Munzel test might be better suited in many cases. The Fligner-Policello test is another test that is sometimes used.
#' 
#' Brunner and Munzel (p. 21) indicate the test-statistic that is computed follows a standard normal distribution, if each category has 50 or more data points. They also remark (p. 22) that the test is no longer accurate if sample sizes are less than 10, although in Schuurhuis et al. (2025, p. 18) 15 is listed. Neubert and Brunner (2007) propose to use a studentized permutation test in these cases. Schuurhuis et al. (2025) developed an improved version of this test as well, called a \(C^2\) test.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/9jo8CrxYndA) and the test is also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Tests/BrunnerMunzelTest.html).
#' 
#' @param catField A vector or dataframe with the group data
#' @param ordField A vector or dataframe with the scores data
#' @param categories optional list with the two categories to use from catField. If not set the first two found will be used
#' @param levels optional list with the scores in order
#' @param distribution c("t", "z") optional distribution to be used
#' 
#' @returns
#' A dataframe with:
#' \item{var. est.}{the estimated overall variance}
#' \item{min n}{the sample size of the category with the lowest count}
#' \item{test-statistic}{the Brunner-Munzel statistic}
#' \item{df}{the degrees of freedom}
#' \item{p-value}{the significance (p-value), two-tailed}
#' \item{categories}{the two categories that were used}
#' \item{test used}{description of the test used and distribution}
#' 
#' @details
#' The test-statistic is calculated using (Brunner & Munzel, 2000, p. 21, eq. 4.8):
#' \deqn{W_n^{BM} = \frac{\bar{R}_2 - \bar{R}_1}{\sqrt{N \times \hat{\sigma}_N^2}}}
#' 
#' with degrees of freedom (Brunner & Munzel, 2000, p. 22):
#' \deqn{df = \frac{\hat{\sigma}_N^4}{N^2 \times \sum_{k=1}^2 \frac{\hat{\sigma}_k^4}{\left(n_k - 1\right)\times n_k^2}}}
#' 
#' The total estimated variance (Brunner & Munzel, 2000, p. 21, eq. 4.7):
#' \deqn{\hat{\sigma}_N^2 = N\times\left(\frac{\hat{\sigma}_1^2}{n_1} + \frac{\hat{\sigma}_2^2}{n_2}\right)}
#' 
#' The estimated variance for each category (Brunner & Munzel, 2000, p. 20, eq. 4.6):
#' \deqn{\hat{\sigma}_k^2 = \frac{s_{R_k^*}^2}{\left(N - n_k\right)^2}}
#' 
#' the sample variance of placement values (Brunner & Munzel, 2000, p. 20, eq. 4.5):
#' \deqn{s_{R_k^*}^2 = \frac{\sum_{i=1}^{n_k} \left(R_{ik}^* - \bar{R}_{k}^*\right)^2}{n_k - 1}}
#' 
#' mean of the placement values:
#' \deqn{\bar{R}_{k}^* = \frac{\sum_{i=1}^{n_k} R_{ik}^*}{n_k}}
#' 
#' the placement values:
#' \deqn{R_{ik}^* = R_{ik} - R_{ik}^{(k)}}
#' 
#' mean of the pooled ranks for each category:
#' \deqn{\bar{R}_k = \frac{\sum_{i=1}^{n_k} R_{ik}}{n_k}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{ik}, the rank of the i-th score in category k, when using all combined scores
#' \item \eqn{R_{ik}^{(k)}}, the rank of the i-th score in category k, when using only scores from category k
#' \item \eqn{N}, the total sample size
#' \item \eqn{n_{k}}, the number of scores in category k
#' }
#' 
#' @section Before, After and Alternatives:
#' Before running the test you might first want to get an impression using a cross table using \code{\link{tab_cross}}, or a stacked bar chart using \code{\link{vi_bar_stacked_multiple}} for a visualisation.
#' 
#' After the test you might want an effect size measure, either \code{\link{es_common_language_is}} for the CLES, \code{\link{me_hodges_lehmann_is}} for Hodges-Lehmann, or \code{\link{r_rank_biserial_is}} for the (Glass) rank biserial (Cliff delta).
#' 
#' Independent samples tests for a binary vs ordinal variable include Brunner-Munzel (\code{\link{ts_brunner_munzel}}), C-square (\code{\link{ts_c_square}}), Cliff-Delta (\code{\link{ts_cliff_delta_is}}), Fligner-Policello (\code{\link{ts_fligner_policello}}), Mann-Whitney U (\code{\link{ts_mann_whitney}}), Mood-Median (\code{\link{ts_mood_median}})
#' 
#' 
#' @references 
#' Brunner, E., & Munzel, U. (2000). The nonparametric Behrens-Fisher problem: Asymptotic theory and a small-sample approximation. *Biometrical Journal, 42*(1), 17-25. https://doi.org/10.1002/(SICI)1521-4036(200001)42:1<17::AID-BIMJ17>3.0.CO;2-U
#' 
#' Neubert, K., & Brunner, E. (2007). A studentized permutation test for the non-parametric Behrens-Fisher problem. *Computational Statistics & Data Analysis, 51*(10), 5192-5204. https://doi.org/10.1016/j.csda.2006.05.024
#' 
#' Schuurhuis, S., Konietschke, F., & Brunner, E. (2025). A new approach to the Nonparametric Behrens-Fisher problem with compatible confidence intervals (arXiv:2504.01796). arXiv. https://doi.org/10.48550/arXiv.2504.01796
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv"
#' dfr <- read.csv(dataFile, sep=";", na.strings=c("", "NA"))
#' coding = c("Far too little", "too little", "Enough", "Too much", "Far too much")
#' ts_brunner_munzel(dfr[['Gen_Gender']], dfr[['Mix_NrAct']], levels=coding)
#' 
#' @export
ts_brunner_munzel <- function(catField, ordField, categories=NULL, levels=NULL, distribution='t'){
  
  #remove rows with missing values
  df = data.frame(ordField, catField)
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
    cat1 = names(table(df[ ,2]))[1]
    cat2 = names(table(df[ ,2]))[2]
  }
  df = na.omit(df)
  #seperate the scores for each category
  X = unname(unlist((subset(df, df[ ,2] == cat1)[1])))
  Y = unname(unlist((subset(df, df[ ,2] == cat2)[1])))
  
  n1 = length(X)
  n2 = length(Y)    
  N = n1 + n2
  minN = min(n1, n2)
  
  #combine this into one long list
  A = c(X, Y)
  
  #get the ranks
  Rik = rank(A)
  Ri1 = rank(X)
  Ri2 = rank(Y)
  
  #pooled ranks for each category
  Rp1 = Rik[1:n1]
  Rp2 = Rik[(n1+1):N]
  
  #placement scores
  Rast1 = Rp1 - Ri1
  Rast2 = Rp2 - Ri2
  
  # mean of pooled ranks
  barR1 = mean(Rp1)
  barR2 = mean(Rp2)
  
  # the variance
  unbVarRast1 = var(Rast1)
  unbVarRast2 = var(Rast2)
  var1 = unbVarRast1/(N - n1)**2
  var2 = unbVarRast2/(N - n2)**2
  var = N*(var1/n1 + var2/n2)
  
  W = (barR2 - barR1)/(N*var)**0.5
  
  if (distribution=='t'){
    if (minN < 10){
      testUsed = 'Brunner-Munzel with t distribution (not recommended)'}
    else {testUsed = 'Brunner-Munzel with t distribution'}
    
    df = var**2 / (N**2 * (var1**2/((n1 - 1)*n1**2) + var2**2/((n2 - 1)*n2**2)))
    p = 2*(1-pt(abs(W), df))
  }
  else if (distribution=='z'){
    if (minN < 50){
      testUsed = 'Brunner-Munzel with standard normal distribution (not recommended)'}
    else {testUsed = 'Brunner-Munzel with standard normal distribution'}
    df = 'n.a.'
    p = 2*(1 - pnorm(abs(W)))
  }
  
  
  cats = paste0(cat1, ', ', cat2)
  
  results <- data.frame(var, minN, W, df, p, cats, testUsed)
  colnames(results)<-c("var. est.", "min n", "test statistic", "df", "p-value", "categories", "test")
  
  return(results)  
}

#' Brunner-Munzel Studentized Permutation Test
#' @description
#' This function performs a studentized permutation test of the Brunner-Munzel test. It is recommended to use this test above the regular version if sample sizes are too small in the data (e.g. less than 10 (or 15)).
#' 
#' @param catField A vector or dataframe with the group data
#' @param ordField A vector or dataframe with the scores data
#' @param categories optional list with the two categories to use from catField. If not set the first two found will be used
#' @param levels optional list with the scores in order
#' @param n_iter optional number of iterations. Default 1000
#' 
#' @returns
#' A dataframe with:
#' \item{iters}{the number of iterations}
#' \item{observed}{the observed Brunner-Munzel statistic}
#' \item{n above}{the number of permutations that had a value above but not equal equal to the observed statistic}
#' \item{n below}{the number of permutations that had a value below but not equal equal to the observed statistic}
#' \item{p-appr.}{the significance (p-value), two-tailed}
#' 
#' @details
#' The idea for a studentized permutation test was proposed by Neubert and Brunner (2007). The regular Brunner-Munzel test statistic is calculated, then the categories get shuffled and the statistic is re-calculated. This repeats.
#' 
#' The minimum of the number of permutations above vs. below the observed statistic is multiplied by two and then divided by the number of iterations, to obtain the p-value (Schuurhuis et al., 2025, p. 7)
#' 
#' See the ts_brunner_munzel for more details.
#' 
#' @section Before, After and Alternatives:
#' Before running the test you might first want to get an impression using a cross table using \code{\link{tab_cross}}, or a stacked bar chart using \code{\link{vi_bar_stacked_multiple}} for a visualisation.
#' 
#' After the test you might want an effect size measure, either \code{\link{es_common_language_is}} for the CLES, \code{\link{es_hodges_lehmann_is}} for Hodges-Lehmann, or \code{\link{r_rank_biserial_is}} for the (Glass) rank biserial (Cliff delta).
#' 
#' @references 
#' Neubert, K., & Brunner, E. (2007). A studentized permutation test for the non-parametric Behrens-Fisher problem. *Computational Statistics & Data Analysis, 51*(10), 5192-5204. https://doi.org/10.1016/j.csda.2006.05.024
#' 
#' Schuurhuis, S., Konietschke, F., & Brunner, E. (2025). A new approach to the Nonparametric Behrens-Fisher problem with compatible confidence intervals (arXiv:2504.01796). arXiv. https://doi.org/10.48550/arXiv.2504.01796
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv"
#' dfr <- read.csv(dataFile, sep=";", na.strings=c("", "NA"))
#' coding = c("Far too little", "too little", "Enough", "Too much", "Far too much")
#' ts_brunner_munzel_perm(dfr[['Gen_Gender']], dfr[['Mix_NrAct']], levels=coding)
#' 
#' @export
ts_brunner_munzel_perm <- function(catField, ordField, categories=NULL, levels=NULL, n_iter=1000){
  
  #remove rows with missing values
  df = data.frame(ordField, catField)
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
    cat1 = names(table(df[ ,2]))[1]
    cat2 = names(table(df[ ,2]))[2]
  }
  df = na.omit(df)
  #seperate the scores for each category
  X = unname(unlist((subset(df, df[ ,2] == cat1)[1])))
  Y = unname(unlist((subset(df, df[ ,2] == cat2)[1])))
  
  n1 = length(X)
  n2 = length(Y)    
  N = n1 + n2
  minN = min(n1, n2)
  
  #combine this into one long list
  A_s = c(X, Y)
  n_below = 0
  n_above = 0
  Ws = c()
  for (i in 1:(n_iter+1)){
    
    # get the shuffled scores for each category
    X_s = A_s[1:n1]
    Y_s = A_s[(n1+1):N]
    # determine the new mid-ranks
    Rik = rank(A_s)
    Ri1 = rank(X_s)
    Ri2 = rank(Y_s)
    # repeat all other steps
    Rp1 = Rik[1:n1]
    Rp2 = Rik[(n1+1):N]
    Rast1 = Rp1 - Ri1
    Rast2 = Rp2 - Ri2
    barR1 = mean(Rp1)
    barR2 = mean(Rp2)
    unbVarRast1 = var(Rast1)
    unbVarRast2 = var(Rast2)
    var1 = unbVarRast1/(N - n1)**2
    var2 = unbVarRast2/(N - n2)**2
    var = N*(var1/n1 + var2/n2)
    W_perm = (barR2 - barR1)/(N*var)**0.5
    Ws = append(Ws, W_perm)
    #print(Rast2)
    # shuffle all scores    
    A_s = sample(A_s)
  }
  
  n_above = sum(Ws > Ws[1])
  n_below = sum(Ws < Ws[1])
  p_appr = 2*min(n_above, n_below)/(n_iter)
  
  results <- data.frame(n_iter, Ws[1], n_above, n_below, p_appr)
  colnames(results) <- c("iters", "observed", "n above", "n below", "p-appr.")
  
  
  return(results)  
}



