#' Mann-Whitney-Wilcoxon Count 
#' @description
#' This function will return the number of possible ways to obtain a specified U value given n1 and n2 cases in each category.
#' 
#' @param u int, the U test statistic
#' @param n1 int, the sample size of the first category
#' @param n2 int the sample size of the second category
#' 
#' @returns
#' result : a list with the counts starting with the count for U=0
#' 
#' @details
#' A recursive formula is used:
#' \deqn{f_{n_1, n_2}(U) = \begin{cases} 0 & \text{ if } U < 0 \text{ or } U > n_1\times n_2 \\ 1 & \text{ if } (n_1=1 \text{ or } n_2=1) \text{ and } 0 \leq U \leq n_1\times n_2 \\ f_{n_1, n_2-1}(U) + f_{n_1-1, n_2}(U - n_2) & \text{ else } \end{cases}}
#' 
#' This formula is found in Mann and Whitney (1947, p. 51), Dinneen and Blakesley 1973, p. 269) and described also in Festinger (1946).
#' 
#' To convert a W statistic to a U statistic use:
#' \deqn{U = W - \frac{n_1\times\left(n_1 + 1\right)}{2}}
#' 
#' @references
#' Dinneen, L. C., & Blakesley, B. C. (1973). Algorithm AS 62: A generator for the sampling distribution of the Mann- Whitney U statistic. *Journal of the Royal Statistical Society. Series C (Applied Statistics), 22*(2), 269–273. doi:10.2307/2346934
#' 
#' Festinger, L. (1946). The significance of difference between means without reference to the frequency distribution function. *Psychometrika, 11*(2), 97–105. doi:10.1007/BF02288926
#' 
#' Mann, H. B., & Whitney, D. R. (1947). On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other. *The Annals of Mathematical Statistics, 18*(1), 50–60. doi:10.1214/aoms/1177730491
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
di_mwwf <- function(u, i, j, memo = list()) {
  # Check if result is already computed
  if (exists(paste(u, i, j, sep = "_"), where = memo)) {
    return(memo[[paste(u, i, j, sep = "_")]])
  }
  
  # Base cases
  if (u < 0 || u > i * j) {
    result <- 0
  } else if (i == 1 || j == 1) {
    result <- 1
  } else {
    # Recursive case
    result <- di_mwwf(u, i, j - 1, memo) + di_mwwf(u - j, i - 1, j, memo)
  }
  
  # Memoize the result
  memo[[paste(u, i, j, sep = "_")]] <- result
  
  return(result)
}

#' Mann-Whitney-Wilcoxon Distribution 
#' @description
#' This distribution is also referred to as a permutation distribution.
#' 
#' It is used in the Mann-Whitney U and Wilcoxon Rank Sum test.
#' 
#' In this version the U-statistic is used as input, and the sample sizes of each of the two categories. This function will return the counts (frequency) of each possible U value from 0 till the provided u value. If all possible values need to be shown, simply set u = i*j.
#' 
#' @param u int, the U test statistic
#' @param n1 int, the sample size of the first category
#' @param n2 int the sample size of the second category
#' 
#' @returns
#' result : a list with the counts starting with the count for U=0
#' 
#' @details
#' A recursive formula is used:
#' \deqn{f_{n_1, n_2}(U) = \begin{cases} 0 & \text{ if } U < 0 \text{ or } U > n_1\times n_2 \\ 1 & \text{ if } (n_1=1 \text{ or } n_2=1) \text{ and } 0 \leq U \leq n_1\times n_2 \\ f_{n_1, n_2-1}(U) + f_{n_1-1, n_2}(U - n_2) & \text{ else } \end{cases}}
#' 
#' This formula is found in Mann and Whitney (1947, p. 51), Dinneen and Blakesley 1973, p. 269) and described also in Festinger (1946).
#' 
#' To convert a W statistic to a U statistic use:
#' \deqn{U = W - \frac{n_1\times\left(n_1 + 1\right)}{2}}
#' 
#' @references
#' Dinneen, L. C., & Blakesley, B. C. (1973). Algorithm AS 62: A generator for the sampling distribution of the Mann- Whitney U statistic. *Journal of the Royal Statistical Society. Series C (Applied Statistics), 22*(2), 269–273. doi:10.2307/2346934
#' 
#' Festinger, L. (1946). The significance of difference between means without reference to the frequency distribution function. *Psychometrika, 11*(2), 97–105. doi:10.1007/BF02288926
#' 
#' Mann, H. B., & Whitney, D. R. (1947). On a Test of Whether one of Two Random Variables is Stochastically Larger than the Other. *The Annals of Mathematical Statistics, 18*(1), 50–60. doi:10.1214/aoms/1177730491
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
di_mwwd <- function(u, i, j) {
  # Create a 3D DP table with dimensions (u+1) x (i+1) x (j+1)
  dp <- array(0, dim = c(u + 1, i + 1, j + 1))
  
  # Initialize base cases
  for (y in 0:i) {
    for (z in 0:j) {
      dp[1, y + 1, z + 1] <- 1
    }
  }
  
  # Fill the DP table based on the recurrence relation
  for (x in 0:u) {
    for (y in 1:i) {
      for (z in 1:j) {
        if (x < 0 || x > y * z) {
          dp[x + 1, y + 1, z + 1] <- 0
        } else if (y == 1 || z == 1) {
          dp[x + 1, y + 1, z + 1] <- 1
        } else {
          # Fill the DP table based on the recurrence relation
          if (x <= z) {
            dp[x + 1, y + 1, z + 1] <- dp[x + 1, y + 1, z] 
            if (x - z >= 0) {
              dp[x + 1, y + 1, z + 1] <- dp[x + 1, y + 1, z + 1] + dp[x - z + 1, y, z + 1]
            }
          } else {
            dp[x + 1, y + 1, z + 1] <- dp[x + 1, y + 1, z]
            if (x - z >= 0) {
              dp[x + 1, y + 1, z + 1] <- dp[x + 1, y + 1, z + 1] + dp[x - z + 1, y, z + 1]
            }
          }
        }
      }
    }
  }
  
  # Return the distribution of results from 0 to u
  return(sapply(0:u, function(x) dp[x + 1, i + 1, j + 1]))
}

#' Mann-Whitney-Wilcoxon Probability Mass Function 
#' @description
#' This function returns the probability for the specified U statistic, given n1 and n2 cases in each category.
#' 
#' It first uses the di_mwwf function to determine the count for the u value, and divides it by the total number of possible arrangements.
#' 
#' The dwilcox() function from R's stats library does the same, and is probably more optimized than this function.
#' 
#' @param u int, the U test statistic
#' @param n1 int, the sample size of the first category
#' @param n2 int the sample size of the second category
#' 
#' @returns
#' result : a list with the counts starting with the count for U=0
#' 
#' @details
#' See the details in di_mwwf() on how the frequency is determined. This is then divided by the total number of possibilities, which is the number of ways we can choose $n_1$ items out of $n$, without replacement. This is the binomial coefficient, or number of combinations:
#' \deqn{C(n, n_1) = nCr(n, n_1) = \binom{n}{n_1} = \frac{n!}{n_1!\times\left(n - n_1\right)!}}
#' 
#' To convert a W statistic to a U statistic use:
#' \deqn{U = W - \frac{n_1\times\left(n_1 + 1\right)}{2}}
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
di_mwwpmf <- function(u, n1, n2) {
  # Calculate the PMF using di_mwwf and choose() for combinations
  return(di_mwwf(u, n1, n2) / choose(n1 + n2, n1))
}

#' Mann-Whitney-Wilcoxon Cumulative Distribution Function
#' @description
#' This function returns the cumulative probability for the specified U statistic, given n1 and n2 cases in each category.
#' 
#' It first uses the di_mwwd function to determine the distribution up to the specified u value, sums these results and divides it by the total number of possible arrangements.
#' 
#' The pwilcox() function from R's stats library does the same, and is probably more optimized than this function.
#' 
#' @param u int, the U test statistic
#' @param n1 int, the sample size of the first category
#' @param n2 int the sample size of the second category
#' 
#' @returns
#' p : the cumulative probability
#' 
#' @details
#' ee the details in di_mwwd() on how the frequency distribution is determined. The sum of these is then divided by the total number of possibilities, which is the number of ways we can choose $n_1$ items out of $n$, without replacement. This is the binomial coefficient, or number of combinations:
#' \deqn{C(n, n_1) = nCr(n, n_1) = \binom{n}{n_1} = \frac{n!}{n_1!\times\left(n - n_1\right)!}}
#' 
#' To convert a W statistic to a U statistic use:
#' \deqn{U = W - \frac{n_1\times\left(n_1 + 1\right)}{2}}
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#'  
#' @export
di_mwwcdf <- function(u, n1, n2) {
  # Determine the entire distribution up to u
  dist <- di_mwwd(u, n1, n2)
  # Calculate the sum of the distribution divided by the total number of possibilities
  p <- sum(dist) / choose(n1 + n2, n2)
  return(p)
}