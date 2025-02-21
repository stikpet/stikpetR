#' P-Value Adjustments for Multiple Testing
#' 
#' @description 
#' Various methods exist to counter a problem with multiple testing. Bonferroni, Šidák, Hommel, Holm, Holm-Šidák, and Hochberg all attempt to control the family wise error rate (FWER), while Benjamini-Hochberg and Benjamini-Yekutieli attempt to control the false discovery rate (FDR).
#' 
#' FWER methods want to minimize the chance of making at least one Type I error (incorrectly rejecting the null hypothesis), while FDR methods attempt to balance the false positive and false negatives.
#' 
#' @param p_values list with the various p-values
#' @param method {'bonferroni', 'sidak', 'hommel', 'holm', 'holm-sidak', 'hochberg', 'bh', 'by', 'hommel-original', 'none'} optional method to use for adjustment, Default is 'bonferroni'
#' @param alpha : float, optional alpha level to use, only applies to 'hommel-original'. Default is 0.05.
#' 
#' @returns
#' p_adj_val : list with the adjusted p-values
#' 
#' @details
#' **none**
#' simply returns the provided p-values
#' 
#' **Bonferroni**
#' 
#' The formula used for the Bonferroni adjustment:
#' \deqn{\tilde{p}_i = \min\left(1, p_i \times k\right)}
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i}, and \eqn{k} the number of tests.
#' 
#' Dunn (1961, p. 53) uses the Bonferroni inequality to adjust confidence intervals, which is why this is also called the Dunn-Bonferroni adjustments.
#' 
#' Bonferroni describes these inequalities in two papers (1935, 1936), but unfortunately I do not read Italian. The term 'Bonferroni inequalities' can already be found in Feller (1950, p. 75)
#' 
#' **Šidák**
#' 
#' The formula used (Šidák, 1967, p. 629):
#' \deqn{\tilde{p}_i = \min\left(1, 1 - \left(1 - p_i\right)^k\right)}
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i}, and \eqn{k} the number of tests.
#' 
#' **Hommel**
#' 
#' The algorithm used (Wright, 1992, p. 1013):
#' 
#' 1. Set \eqn{a_i = p_i} for all \eqn{i}
#' 2. For each \eqn{m = k, (k - 1), ..., 2} (i.e. in descending order) do the following:
#'     2.1. For \eqn{i > (k - m)}
#'         * 2.1.1. calculate \eqn{c_i = \frac{m\times p_i}{m + i - k}}
#'         * 2.1.2. set \eqn{c_{min} = \min\left(c_i, for all i\right)}
#'         * 2.1.3. if \eqn{a_i < c_{min}} then \eqn{a_i = c_min}
#'     2.2. For \eqn{i \leq (k - m)}
#'         * 2.2.1. let \eqn{c_i = \min\left(c_{min}, m \times p_i\right)}
#'         * 2.2.2. if \eqn{a_i < c_i} then \eqn{a_i = ci}
#' 3. \eqn{p_i^{adj} = a_i}
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i} after sorting all p-values in ascending order, and \eqn{k} the number of tests.
#' 
#' Hommel (1988) original procedure is I think slightly different and implemented in 'hommel-original'. The method from Wright seems to be the method used in the multipletests() function from the Python library statsmodels.stats.multitest, and the p.adjust() function from R's stats library. The advantage of Wright's algorithm, is that it doesn't require the alpha level to be known to adjust the p-values.
#' 
#' **Holm**
#' 
#' The formula used (SAS, n.d.):
#' \deqn{\tilde{p}_i =\begin{cases} k\times p_1 & i= 1\\ \max \left(\tilde{p}_{i-1}, p_i \times \left(k+1-i\right)\right) & i=2,\dots,k \end{cases} }
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i} after sorting all p-values in ascending order, and \eqn{k} the number of tests.
#' 
#' Holm (1979, p. 67) describes this procedure, but uses alpha level.
#' 
#' **Holm-Šidák**
#' 
#' The formula used (SAS, n.d.):
#' \deqn{\tilde{p}_i =\begin{cases} 1 - \left(1 - p_1 \right)^k & i= 1\\ \max \left(\tilde{p}_{i-1}, 1 - \left(1 - p_i\right)^{k - i + 1}\right) & i=2,\dots,k \end{cases} }
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i} after sorting all p-values in ascending order, and \eqn{k} the number of tests.
#' 
#' This uses Holm (1979, p. 67) step-down approach, but instead of using the Bonferroni adjustment, it uses Šidák.
#' 
#' **Hochberg**
#' 
#' The formula used (SAS, n.d.):
#' \deqn{\tilde{p}_i =\begin{cases} p_i & i= 1\\ \min \left(\tilde{p}_{i-1}, i \times p_i\right) & i=2,\dots,k \end{cases} }
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i} after sorting all p-values in DESCENDING order, and \eqn{k} the number of tests.
#' 
#' The procedure is described by Hochberg (1988, p. 801) using alpha levels for the criteria.
#' 
#' **Benjamini-Hochberg**
#' 
#' The algorithm used (Benjamini & Hochberg, 1995, p. 293):
#' 
#' 1. Sort the p-values in ascending order
#' 2. find \eqn{j = \max\left(i | p_i \times \frac{k}{i} \leq \alpha\right)}
#' 3. All tests \eqn{i \leq j} are considered significant.
#' 
#' To find the adjusted p-values (in reverse order):
#' 
#' \deqn{\tilde{p}_i =\begin{cases} p_k & i= k\\ \min \left(\tilde{p}_{i+1}, p_i \times \frac{k}{i}\right) & i=(k-1),\dots,1 \end{cases} }
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i} after sorting all p-values in ascending order, and \eqn{k} the number of tests.
#' 
#' **Benjamini-Yekutieli**
#' 
#' The algorithm used (Benjamini & Yekutieli, 2001, p. 1169):
#' 1. Sort the p-values in ascending order
#' 3. Determine \eqn{C(k) = \sum_{i=1}^k \frac{1}{i}}
#' 2. find \eqn{j = \max\left(i | p_i \times \frac{k\times C(k)}{i} \leq \alpha\right)}
#' 3. All tests \eqn{i \leq j} are considered significant.
#' 
#' To find the adjusted p-values (in reverse order):
#' \deqn{\tilde{p}_i =\begin{cases} p_k\times C(k) & i= k\\ \min \left(\tilde{p}_{i+1}, p_i \times \frac{k\times C(k)}{i}\right) & i=(k-1),\dots,1 \end{cases} }
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i} after sorting all p-values in ascending order, and \eqn{k} the number of tests.
#' 
#' **Hommel Original**
#' 
#' Hommel (1988, p. 384) describes the following algorithm:
#' 
#' 1. Compute \eqn{j = \max\left(i \in {1,\dots, k}| p_{k-i+j} > j\times\frac{\alpha}{i} for j=1,\dots, i\right)}
#' 2. If the maximum does not exist, reject all, otherwise reject all with \eqn{p_i \leq \frac{\alpha}{j}}
#' 
#' Where \eqn{p_i} is the p-value of test \eqn{i} after sorting all p-values in ascending order, and \eqn{k} the number of tests.
#' 
#' The function will adjust the p-values using:
#' 
#' \deqn{\tilde{p}_i =\begin{cases} \min\left(1, j\times p_i\right) & j exists \\ 1 & j does not exist \end{cases} }
#' 
#' @references
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B: Statistical Methodology, 57*(1), 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x
#' 
#' Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. *The Annals of Statistics, 29*(4), 1165–1188. doi:10.1214/aos/1013699998
#' 
#' Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. In *Studi in Onore del Professore Salvatore Ortu Carboni* (pp. 13–60).
#' 
#' Bonferroni, C. E. (1936). *Teoria statistica delle classi e calcolo delle probabilità*. Pubblicazioni Del R Istituto Superiore Di Scienze Economiche e Commerciali Di Firenze, 8, 3–62.
#' 
#' Dunn, O. J. (1961). Multiple comparisons among means. *Journal of the American Statistical Association, 56*(293), 52–64. https://doi.org/10.1080/01621459.1961.10482090
#' 
#' Feller, W. (1950). *An introduction to probability theory and its applications: Vol. One*. John Wiley & Sons.
#' 
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. *Biometrika, 75*(4), 800–802. doi:10.1093/biomet/75.4.800
#' 
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure. *Scandinavian Journal of Statistics, 6*(2), 65–70.
#' 
#' Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. *Biometrika, 75*(2), 383–386. doi:10.1093/biomet/75.2.383
#' 
#' SAS. (n.d.). PROC MULTTEST: p-Value Adjustments. SAS/STAT(R) 9.22 User’s Guide. Retrieved February 1, 2025, from https://support.sas.com/documentation/cdl/en/statug/63347/HTML/default/viewer.htm#statug_multtest_sect014.htm
#' 
#' Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. *Journal of the American Statistical Association, 62*(318), 626. doi:10.2307/2283989
#' 
#' Wright, S. P. (1992). Adjusted p-values for simultaneous inference. *Biometrics, 48*(4), 1005. doi:10.2307/2532694
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
p_adjust <- function(p_values, method='bonferroni', alpha=.05){
  k = length(p_values)
  
  if (method=="none"){
    p_adj_val = p_values
  }
  else if (method=="bonferroni"){
    p_adj_val = pmin(1, p_values*k)
  }
  else if (method=='sidak'){
    #calculate the Sidak adjustments
    p_adj_val = pmin(1, 1 - (1 - p_values)^k)
  }
  else{
    # Create a data frame with index and values
    p_data <- data.frame(index = seq(1, k), value = p_values)
    # Sort by p-value (second column)
    p_sorted <- p_data[order(p_data$value), ]
    
    if (method=='hommel'){
      ai = p_sorted
      ai[, 3] = ai[, 2]
      
      ci = rep(0, k)
      for (m in seq(k, 2, -1)) {
        c_min = 1
        for (i in (k-m+1):k) {                
          ci[i] = m*p_sorted[i, 2]/(m + i - k)
          if (ci[i] < c_min){
            c_min = ci[i]
          }     
        }
        for (i in (k-m+1):k){
          if (ai[i, 3] < c_min){
            ai[i, 3] = c_min
          }
        }
        
        if (k - m > 0){
          for (i in 1:(k - m)){
            ci[i] = min(c_min, m*p_sorted[i,2])
            if (ai[i,3] < ci[i]){
              ai[i,3] = ci[i]
            }
          }
        }
      }
    }
    else if (method=='holm'){            
      #add the Holm adjusted p-values
      ai = p_sorted
      ai[,3] = ai[,2]
      ai[1,3] = min(1, k*ai[1,2])
      for (i in 2:k){
        ai[i,3] = max(ai[i-1,3], min(1, ai[i,2]*(k + 1 - i)))
      }
    }
    else if (method=='holm-sidak'){            
      #add the Holm-Sidak adjusted p-values
      ai = p_sorted
      ai[,3] = ai[,2]
      ai[1,3] = min(1, 1 - (1 - ai[1,2])**k)
      for (i in 2:k){
        ai[i,3] = max(ai[i-1,3], min(1, 1 - (1 - ai[i,2])**(k + 1 - i)))
      }
    }
    else if (method=='hochberg'){  
      p_sorted_desc <- p_data[order(p_data$value, decreasing = TRUE), ]
      #add the Holm adjusted p-values
      ai = p_sorted_desc
      ai[,3] = ai[,2]
      ai[1,3] = min(1, ai[1,2])
      for (i in 2:k){
        ai[i,3] = min(ai[i-1,3], min(1, i*ai[i,2]))
      }            
    }
    else if (method=='bh'){            
      
      ai = p_sorted
      ai[,3] = ai[,2]
      ai[k,3] = ai[k,2]
      for (i in seq(k-1, 1, -1)){
        ai[i,3] = min(ai[i+1,3], min(1, ai[i,2]*k/i))
      }
    }
    else if (method=='by'){            
      
      ai = p_sorted
      ck = sum(1/seq(1,k))
      ai[,3] = ai[,2]
      ai[k,3] = min(1, ai[k,2]*ck)
      for (i in seq(k-1, 1, -1)){
        ai[i,3] = min(ai[i+1,3], min(1, ai[i,2]*k*ck/i))
      }
    }
    
    else if (method=='hommel-original'){            
      
      i_hommel = -1
      i = 1
      while (i <= k && i_hommel== -1){
        c_min = 1
        for (j in 1:i){
          if (p_sorted[k-i+j,2]*i/j< c_min){
            c_min = p_sorted[k-i+j,2]*i/j
          }
        }
        if (c_min > alpha){
          i_hommel = i
        }
        i = i + 1
      }
      
      ai = p_sorted
      if (i_hommel > 0){
        for (i in 1:k){
          ai[i,3] = min(1, ai[i,2]*i_hommel)
        }
      }
      else {
        for (i in 1:k){
          ai[i,3] = 1
        }
      }
    }
    
    p_sort_back <- ai[order(ai$index), ]
    p_adj_val = as.numeric(p_sort_back[, 3])
  }
  
  return (p_adj_val)
  
}