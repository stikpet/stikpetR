#' Independent Samples Yuen(-Welch) Mean Test
#' 
#' @param scores A vector with the numeric scores
#' @param groups A vector with the group data
#' @param trimProp optional proportion to trim from each side (so in total twice this will be trimmed)
#' @return dataframe test statistic, degrees of freedom, p-value (sig.) and name of test used
#'  
#' @details 
#' The formula used is (Yuen, 1974, p. 167):
#' \deqn{t = \frac{\bar{x}_{t,1} - \bar{x}_{t,2}}{SE}}
#' \deqn{sig = 2\times\left(1 - T\left(\left|t\right|, df\right)\right)}
#' 
#' With:
#' \deqn{SE = \sqrt{\frac{s_{w,1}^2}{m_1} + \frac{s_{w,2}^2}{m_2}}}
#' \deqn{s_{w,i}^2 = \frac{SSD_{w,i}}{m_i - 1}}
#' \deqn{df = \frac{1}{\frac{c^2}{m_1 - 1} + \frac{\left(1 - c\right)^2}{m_2 -1}}}
#' \deqn{c = \frac{\frac{s_{w,1}^2}{m_1}}{\frac{s_{w,1}^2}{m_1} + \frac{s_{w,2}^2}{m_2}}}
#' \deqn{\bar{x}_{t,i} = \frac{\sum_{j=g_i+1}^{n_i - g_i}y_{i,j}}{}}
#' \deqn{g_i = \lfloor n_i\times p_t\rfloor}
#' \deqn{m_i = n_ - 2\times g_i}
#' \deqn{SSD_{w,i} = g_i\times\left(y_{i,g_i+1} - \bar{x}_{wi}\right)^2 + g_i\times\left(y_{i,n_i-g_i} - \bar{x}_{w,i}\right)^2 + \sum_{j=g+1}^{n_i - g_i} \left(y_{i,j} - \bar{x}_{w,i}\right)^2}
#' \deqn{\bar{x}_{w,i} = \frac{\bar{x}_{t,i}\times m_i + g_i\times\left(y_{i, g_i+1} + y_{i, n_i-g_i}\right)}{n_i}}
#' 
#' *Symbols used:*
#' \itemize{
#' \item \eqn{x_{t,i}} the trimmed mean of the scores in category i
#' \item \eqn{x_{w,i}} The Winsorized mean of the scores in category i
#' \item \eqn{SSD_{w,i}} the sum of squared deviations from the Winsorized mean of category i
#' \item \eqn{m_i} the number of scores in the trimmed data set from category i
#' \item \eqn{y_{i,j}} the j-th score after the scores in category i, after they are sorted from low to high
#' \item \eqn{p_t} the proportion of trimming on each side, we can define
#' }
#' 
#' @author 
#' P. Stikker
#' 
#' Please visit: https://PeterStatistics.com
#' 
#' YouTube channel: https://www.youtube.com/stikpet
#' 
#' @references 
#' Yuen, K. K. (1974). The two-sample trimmed t for unequal population variances. *Biometrika, 61*(1), 165–170. https://doi.org/10.1093/biomet/61.1.165
#'  
#' @examples 
#' scores = c(20,50,80,15,40,85,30,45,70,60, NA, 90,25,40,70,65, NA, 70,98,40)
#' groups = c("national","international","international","national","international", "international","national","national","international","international","international","international","international","international","national", "international" ,NA,"national","international","international")
#' ts_yuen_is(scores, groups)
#' 
#' @export
ts_yuen_is <- function(scores, groups, trimProp = 0.1){
    df = data.frame(scores, groups)
    df = na.omit(df)
    
    #set trimming proportion
    trimProp = 0.1
    scores1 = df$scores[df$groups==df$groups[1]]
    n1 = length(scores1)
    n1t = n1 - 2*round(n1*trimProp)
    m1t = mean(scores1, trim=trimProp)
    #note that R base mean(..., trim) will round down
    
    scores2 = df$scores[df$groups!=df$groups[1]]
    n2 = length(scores2)
    n2t = n2 - 2*round(n2*trimProp)
    m2t = mean(scores2, trim=trimProp)
    
    #Winsorize the data
    sc1Sort = sort(scores1)
    minReplace1 = sc1Sort[round(n1*trimProp)+1]
    maxReplace1 = sc1Sort[n1 - round(n1*trimProp)]
    scores1 = replace(scores1, scores1<minReplace1, minReplace1)
    scores1 = replace(scores1, scores1>maxReplace1, maxReplace1)
    
    sc2Sort = sort(scores2)
    minReplace2 = sc2Sort[round(n2*trimProp)+1]
    maxReplace2 = sc2Sort[n2 - round(n2*trimProp)]
    scores2 = replace(scores2, scores2<minReplace2, minReplace2)
    scores2 = replace(scores2, scores2>maxReplace2, maxReplace2)
    
    #now for the test
    wssd1 = var(scores1)*(n1 - 1)
    wssd2 = var(scores2)*(n2 - 1)
    
    wvar1 = wssd1/(n1t - 1)
    wvar2 = wssd2/(n2t - 1)
    
    d1 = wvar1/n1t
    d2 = wvar2/n2t
    
    se = sqrt(d1 + d2)
    
    tValue = (m1t - m2t)/se
    
    c = d1/(d1 + d2)
    df = 1/(c^2/(n1t - 1) + (1 - c)^2/(n2t - 1))
    
    pValue = 2*(1 - pt(abs(tValue), df))
    
    statistic = tValue
    testResults <- data.frame(statistic, df, pValue)
  
    return(testResults)
    
}