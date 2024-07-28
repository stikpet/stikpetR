#' Wilcoxon Cumulative Distribution Function
#' @description
#' This function will give the cumulative probability of a sum of ranks of T, given a sample size of n.
#' 
#' @param T int with the sum of ranks
#' @param n int with the sample size
#' @param method optional the calculation method to use. Either "shift" (default), "enumerate", "recursive".
#' @returns
#' A float with the requested probability
#' @details 
#' The enumeration method will create all possible combinations of ranks 1 to n, sum each of these, and then determines the count of each unique sum of ranks. It then uses this to determine the probability and cumulative probabilities.
#' 
#' The recursive method uses the formula from McCornack (1965, p. 864):
#' \deqn{srf\left(x,y\right)=\begin{cases} 0 & x < 0 \\ 0 & x > \binom{y+1}{2} \\ 1 & y=1 \wedge \left( x=0\vee x=1 \right) \\ srf^*\left(x,y\right) & y\geq 0 \end{cases}}
#' 
#' with:
#' \deqn{srf^*\left(x,y\right) = srf\left(x-y,y-1\right) + srf\left(x,y-1\right)}
#' 
#' The shift-algorithm from Streitberg and Röhmel (1987), and can also be found in Munzel and Brunner (2002). This works as follows.
#' \itemize{
#' \item Start with listing all values from 0 to the maximum possible sum of ranks, so 0 to (n×(n+1))/2
#' \item Create a vector with the value 1 followed by n times a 0.
#' \item Create a shifted vector by moving all values by 1.
#' \item Add the two results (the original and the shifted version)
#' \item This will be the updated vector
#' \item Shift the vector now by 2
#' \item Add the two results (the updated vector with and the two shifted version)
#' \item Repeat these steps each time shifting by one more than the previous. Stop when n-times shifting has been done.
#' }
#' @references 
#' McCornack, R. L. (1965). Extended tables of the Wilcoxon matched pair signed rank statistic. *Journal of the American Statistical Association, 60*(311), 864–871. doi:10.2307/2283253
#' 
#' Munzel, U., & Brunner, E. (2002). An exact paired rank test. *Biometrical Journal, 44*(5), 584. doi:10.1002/1521-4036(200207)44:5<584::AID-BIMJ584>3.0.CO;2-9
#' 
#' Streitberg, B., & Röhmel, J. (1987). Exakte Verteilungen für Rang-und Randomisierungstests im allgemeinen c-Stichprobenproblem. *EDV in Medizin und Biologie, 18*(1), 12–19.
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
di_wcdf <- function(T, n, method="shift"){
  
  if (method=="recursive"){
    # sum all pmf results from 0 to T
    pVal=0
    for (i in 0:T){
      pVal = pVal + srf(i, n)}
    wcdf = pVal/(2**n)        
  }
  
  else if (method=="enumerate"){
    #create all possible combinations of n ranks
    rank_dist = t(expand.grid(rep(list(0:1), n)))
    
    #list all possible ranks (1 to n)
    ranks = 1:n
    
    #sum all ranks used per column
    sum_rank_dist = colSums(rank_dist*ranks)
    
    #determine the maximum rank sum
    sr_max = n*(n+1)/2
    
    #determine the count for each unique sum of ranks
    sum_rank_freq = table(sum_rank_dist)
    
    #divide by sum of all counts
    sum_rank_prop = sum_rank_freq / (2**n)
    
    #sum all results up to and including T (add one since array start at 1 in R)
    wcdf = sum(sum_rank_prop[1:(T+1)])}
  
  else if (method=="shift"){
    
    #determine the maximum rank sum
    max_rank = n*(n+1)/2
    
    #Start with listing all values from 0 to the maximum possible sum of ranks
    possible_values = list(0:max_rank)
    
    #Create a vector with the value 1 followed by n times a 0
    freqs = c(1, rep(0, max_rank))
    
    for (i in 1:n){
      #Create a shifted vector by moving all values by i
      shift = c(tail(freqs, i), head(freqs, max_rank + 1 - i))
      
      #Add the two results 
      freqs = freqs + shift
    }
    
    n_combs = sum(freqs)
    wcdf = sum(freqs[1:(T+1)])/n_combs}
  
  return (wcdf)
  
}

#' Wilcoxon Probability Mass Function
#' @description
#' This function will give the probability of a sum of ranks of T, given a sample size of n.
#' 
#' @param T int with the sum of ranks
#' @param n int with the sample size
#' @param method optional the calculation method to use. Either "shift" (default), "enumerate", "recursive".
#' @returns
#' A float with the requested probability
#' @details 
#' The enumeration method will create all possible combinations of ranks 1 to n, sum each of these, and then determines the count of each unique sum of ranks. It then uses this to determine the probability.
#' 
#' The recursive method uses the formula from McCornack (1965, p. 864):
#' \deqn{srf\left(x,y\right)=\begin{cases} 0 & x < 0 \\ 0 & x > \binom{y+1}{2} \\ 1 & y=1 \wedge \left( x=0\vee x=1 \right) \\ srf^*\left(x,y\right) & y\geq 0 \end{cases}}
#' 
#' with:
#' \deqn{srf^*\left(x,y\right) = srf\left(x-y,y-1\right) + srf\left(x,y-1\right)}
#' 
#' The shift-algorithm from Streitberg and Röhmel (1987), and can also be found in Munzel and Brunner (2002). This works as follows.
#' \itemize{
#' \item Start with listing all values from 0 to the maximum possible sum of ranks, so 0 to (n×(n+1))/2
#' \item Create a vector with the value 1 followed by n times a 0.
#' \item Create a shifted vector by moving all values by 1.
#' \item Add the two results (the original and the shifted version)
#' \item This will be the updated vector
#' \item Shift the vector now by 2
#' \item Add the two results (the updated vector with and the two shifted version)
#' \item Repeat these steps each time shifting by one more than the previous. Stop when n-times shifting has been done.
#' }
#' @references 
#' McCornack, R. L. (1965). Extended tables of the Wilcoxon matched pair signed rank statistic. *Journal of the American Statistical Association, 60*(311), 864–871. doi:10.2307/2283253
#' 
#' Munzel, U., & Brunner, E. (2002). An exact paired rank test. *Biometrical Journal, 44*(5), 584. doi:10.1002/1521-4036(200207)44:5<584::AID-BIMJ584>3.0.CO;2-9
#' 
#' Streitberg, B., & Röhmel, J. (1987). Exakte Verteilungen für Rang-und Randomisierungstests im allgemeinen c-Stichprobenproblem. *EDV in Medizin und Biologie, 18*(1), 12–19.
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
di_wpmf <- function(T, n, method="shift"){
  
  if (method=="recursive"){
    wpmf = srf(T, n)/(2**n)        
  }
  
  else if (method=="enumerate"){
    #create all possible combinations of n ranks
    rank_dist = t(expand.grid(rep(list(0:1), n)))
    
    #list all possible ranks (1 to n)
    ranks = 1:n
    
    #sum all ranks used per column
    sum_rank_dist = colSums(rank_dist*ranks)
    
    #determine the maximum rank sum
    sr_max = n*(n+1)/2
    
    #determine the count for each unique sum of ranks
    sum_rank_freq = table(sum_rank_dist)
    
    #divide by sum of all counts
    sum_rank_prop = sum_rank_freq / (2**n)
    
    #find prop for T (add one since array start at 1 in R)
    wpmf = unname(sum_rank_prop[T+1])}
  
  else if (method=="shift"){
    
    #determine the maximum rank sum
    max_rank = n*(n+1)/2
    
    #Start with listing all values from 0 to the maximum possible sum of ranks
    possible_values = list(0:max_rank)
    
    #Create a vector with the value 1 followed by n times a 0
    freqs = c(1, rep(0, max_rank))
    
    for (i in 1:n){
      #Create a shifted vector by moving all values by i
      shift = c(tail(freqs, i), head(freqs, max_rank + 1 - i))
      
      #Add the two results 
      freqs = freqs + shift
    }
    
    n_combs = sum(freqs)
    wpmf = freqs[T+1]/n_combs}
  
  return (wpmf)
  
}
#' Wilcoxon Sum-of-Ranks Frequency Function
#' @description
#' This helper function will give the count for a sum of ranks of T, given a sample size of n, using the recursive formula.
#' 
#' @param T int with the sum of ranks
#' @param n int with the sample size
#' @returns
#' A the requested count
#' @details 
#' The recursive method uses the formula from McCornack (1965, p. 864):
#' \deqn{srf\left(x,y\right)=\begin{cases} 0 & x < 0 \\ 0 & x > \binom{y+1}{2} \\ 1 & y=1 \wedge \left( x=0\vee x=1 \right) \\ srf^*\left(x,y\right) & y\geq 0 \end{cases}}
#' 
#' with:
#' \deqn{srf^*\left(x,y\right) = srf\left(x-y,y-1\right) + srf\left(x,y-1\right)}
#' 
#' }
#' @references 
#' McCornack, R. L. (1965). Extended tables of the Wilcoxon matched pair signed rank statistic. *Journal of the American Statistical Association, 60*(311), 864–871. doi:10.2307/2283253
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
srf <-function(x,y){
  #srf = sum rank frequency
  if (x<0){
    return (0)}
  else if (x > y*(y+1)/2){
    return (0)}
  else if (y==1 && (x==0 || x==1)){
    return (1)}
  else if (y>=0){
    return (srf(x-y, y-1) + srf(x, y-1))}
}
