#' Pairwise Binary-Ordinal Effect Sizes
#' @description 
#' This function determines the effect size for each comparison in a post-hoc analysis of a nominal vs. ordinal variable (e.g. a Kruskal-Wallis test).
#' 
#' @param catField vector with categories
#' @param ordField vector with the scores
#' @param categories vector, optional. the categories to use from catField
#' @param levels vector, optional. the levels or order used in ordField.
#' @param es string, optional. the effect size to determine.Either "cle" (default), "rb" or "rosenthal"
#' 
#' @returns
#' dataframe with 
#' \item{cat1}{label of first category in pair}
#' \item{cat2}{label of second category in pair}
#' \item{effect size}{the value of the effect size}
#' 
#' @details
#' The function simply goes over each possible pair of categories from the *catField* (adjusted with *categories* if used). It then runs for only the scores of those two categories the Common Language Effect Size (Vargha-Delaney A) or (Glass) Rank Biserial (Cliff delta). If the Rosenthal correlation is requested, it will perform the post-hoc Dunn test to obtain the z-statistic.
#' 
#' @seealso 
#' \code{\link{es_common_language_is}}, Common Language Effect size
#' \code{\link{r_rank_biserial_is}}, rank biserial for independent samples
#' \code{\link{ph_dunn}}, post-hoc Dunn test, used to obtain z-value for Rosenthal correlation
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
es_pairwise_bin_ord  <-
  function(catField,
           ordField,
           categories = NULL,
           levels = NULL,
           es = "cle") {
    cat1 <- c()
    cat2 <- c()
    p_res <- c()
    
    if (es %in% c('cle', 'rb')) {
      selDf = na.omit(data.frame("groups" = catField, "scores" = ordField))
      
      #only keep the ones in given list
      if (!is.null(categories)) {
        selDf <- selDf[selDf$groups %in% categories,]
      }
      
      #get the unique categories
      cats <- unique(selDf$groups)
      
      #number of categories
      k = length(cats)
      p_res2 <- c()
      for (i in 1:(k - 1)) {
        for (j in (i + 1):k) {
          cat1 <- append(cat1, cats[i])
          cat2 <- append(cat2, cats[j])
          sel_cats <- c(cats[i], cats[j])
          
          if (es == 'cle') {
            es_res = es_common_language_is(selDf["groups"],
                                         selDf["scores"],
                                         categories = sel_cats,
                                         levels = myCoding, 
                                         method="vda")
            p_res <- append(p_res, c(as.matrix(es_res[1,]))[1])
            p_res2 <- append(p_res2, c(as.matrix(es_res[1,]))[2])
          }
          else{
            es_res = r_rank_biserial_is(selDf["groups"],
                                        selDf["scores"],
                                        categories = sel_cats,
                                        levels = levels)
            p_res <- append(p_res, es_res)
          }
        }
      }
      
      if (es == 'cle') {
        results = data.frame(cat1, cat2, p_res, p_res2)
        colnames(results) = c('cat. 1', 'cat. 2', 'CLE 1', 'CLE 2')
      }
      else{
        results = data.frame(cat1, cat2, p_res)
        colnames(results) = c('cat. 1', 'cat. 2', 'rb')
      }
    }
    else {
      ph_dunn_res = ph_dunn(catField,
                            ordField,
                            categories = categories,
                            levels = levels)
      k = dim(ph_dunn_res)[1]
      ph_dunn_res['n'] = ph_dunn_res['n1'] + ph_dunn_res['n2']
      p_res <-
        sapply(1:k, function(i) {
          r_rosenthal(ph_dunn_res$statistic[i], ph_dunn_res$n[i])
        })
      cat1 = ph_dunn_res['cat. 1']
      cat2 = ph_dunn_res['cat. 2']
      
      results = data.frame(cat1, cat2, p_res)
      colnames(results) = c('cat. 1', 'cat. 2', 'Rosenthal Correlation')
    }
    
    return (results)
  }
