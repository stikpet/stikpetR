#' Exact Multinomial Test of Goodness-of-Fit
#' 
#' @description 
#' A test that can be used with a single nominal variable, to test if the probabilities in all the categories 
#' are equal (the null hypothesis). If the test has a p-value below a pre-defined threshold (usually 0.05) the
#' assumption they are all equal in the population will be rejected. 
#' 
#' There are quite a few tests that can do this. Perhaps the most commonly used is a Pearson chi-square test, 
#' but also a G-test, Freeman-Tukey, Neyman, Mod-Log Likelihood and Cressie-Read test are possible.
#' 
#' McDonald (2014, p. 82) suggests to always use this exact test as long as the sample size is less than 1000 
#' (which was just picked as a nice round number, when n is very large the exact test becomes 
#' computational heavy even for computers).
#' 
#' @param data A vector with the data
#' @param expCounts Optional dataframe with the categories and expected counts 
#' 
#' @returns 
#' Dataframe with:
#' \item{pObs}{probability of the observed data}
#' \item{ncomb}{number of combinations used}
#' \item{pValue}{two-sided p-value}
#' \item{testUsed}{a description of the test used}
#'  
#' @details 
#' The exact multinomial test of goodness of fit is done in four steps
#' 
#' Step 1: Determine the probability of the observed counts using 
#' the probability mass function of the multinomial distribution
#' 
#' Step 2: Determine all possible permutations with repetition 
#' that create a sum equal to the sample size over the k-categories.
#' 
#' Step 3: Determine the probability of each of these permutations using the probability mass function of the multinomial distribution.
#' 
#' Step 4: Sum all probabilities found in step 3 that are equal or less than the one found in step 1.
#' 
#' @section Alternatives:
#' 
#' The *EMT* library has a similar function: *multinomial.test()*
#' 
#' The *XNomial* library has a similar function: *xmulti()*
#' 
#' @references 
#' McDonald, J. H. (2014). *Handbook of biological statistics* (3rd ed.). Sparky House Publishing.
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1[1:20, 'mar1']
#' ts_multinomial_gof(ex1)
#' 
#' #Example 2: dataframe with various settings
#' ex2 = df1[1:20, 'mar1']
#' eCounts = data.frame(c("MARRIED", "DIVORCED", "NEVER MARRIED", "SEPARATED"), c(5,5,5,5))
#' ts_multinomial_gof(ex2, expCounts=eCounts)
#' 
#' #Example 3: a list
#' ex3 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", 
#' "DIVORCED", "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", 
#' "DIVORCED", "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' ts_multinomial_gof(ex3)
#' 
#' @export
ts_multinomial_gof <- function(data, expCounts=NULL) {

  if (is.list(data)) {
    data <- unlist(data)
  }
  
  data = na.omit(data)
  
  # Determine the observed counts
  if (is.null(expCounts)) {
    # Generate frequency table
    freq <- as.data.frame(table(data))
    colnames(freq) <- c("category", "count")
    n <- sum(freq$count)
    
    # Number of categories to use (k)
    k <- nrow(freq)
    
    # Number of expected counts is simply sample size
    nE <- n
  } else {
    # If expected counts are given
    k <- nrow(expCounts)
    
    freq <- data.frame(category = character(), count = integer(), stringsAsFactors = FALSE)
    for (i in 1:k) {
      nk <- sum(data == expCounts[i, 1])
      lk <- expCounts[i, 1]
      freq <- rbind(freq, data.frame(category = lk, count = nk))
    }
    nE <- sum(expCounts[, 2])
  }
  
  n <- sum(freq$count)
  
  # The true expected counts
  if (is.null(expCounts)) {
    # Assume all to be equal
    exp_prop <- rep(1 / k, k)
  } else {
    # Check if categories match
    exp_prop <- numeric(k)
    for (i in 1:k) {
      exp_prop[i] <- expCounts[i, 2] / nE
    }
  }
  
  observed <- freq$count
  p_obs <- dmultinom(sort(observed), size = n, prob = exp_prop)
  counts <- 0:n
  
  # Generate all permutations of counts
  all_perm <- as.matrix(expand.grid(rep(list(counts), k)))
  sum_perm <- all_perm[rowSums(all_perm) == n, ]
  ncomb <- nrow(sum_perm)
  
  # Compute the p-value
  p_val <- 0
  for (i in 1:ncomb) {
    p_perm <- dmultinom(sum_perm[i, ], size = n, prob = exp_prop)
    if (p_perm <= p_obs) {
      p_val <- p_val + p_perm
    }
  }
  
  testUsed <- "one-sample multinomial exact goodness-of-fit test"
  testResults <- data.frame(
    `p obs.` = p_obs,
    `n combs.` = ncomb,
    `p-value` = p_val,
    `test` = testUsed
  )
  
  return(testResults)
}