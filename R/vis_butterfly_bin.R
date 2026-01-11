#' Butterfly Chart of Binned Data (Pyramid chart)
#' @description
#' This function creates a simple butterfly chart by binning scaled data. This is sometimes also referred to as a Pyramid chart. 
#' 
#' @param catField list or dataframe with the categories
#' @param scaleField list or dataframe with the scores
#' @param categories optional list with the two categories to use from bin_field. If not set the first two found will be used
#' @param bins optional either a list with parameters to pass to tab_frequency_bins, a dataframe with upper and lower bounds, or simply NULL (default)
#' @param ... optional additional parameters to pass to vi_butterfly_chart
#' 
#' @returns
#' butterfly chart
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @export
vi_butterfly_bin <- function(catField, scaleField, categories=NULL, 
                             bins=NULL, ...){
  # DATA PREPARATION
  # create dataframe and remove missing values
  df <- na.omit(data.frame(score=scaleField, category=catField))
  
  #the two categories
  if (!is.null(categories)){
    #use the provided categories
    cat1 = categories[1]
    cat2 = categories[2]
  }
  else {
    # use first two categories found
    cat1 = names(table(df[ ,2]))[1]
    cat2 = names(table(df[ ,2]))[2]
  }
  
  # DETERMINE BREAKS
  #split the scores across the categories    
  scoresCat1 = unname(unlist((subset(df, df[ ,2] == cat1)[1])))
  scoresCat2 = unname(unlist((subset(df, df[ ,2] == cat2)[1])))
  #combine this into one long list
  allScores = c(scoresCat1, scoresCat2)
  
  #determine bins overall
  freq_table <- switch(
    class(bins)[1],   # use the first class if bins has multiple
    "NULL" = tab_frequency_bins(allScores),
    "list" = do.call(tab_frequency_bins, c(list(data = allScores), bins)),
    "data.frame" = tab_frequency_bins(allScores, bins = bins),
    stop("Invalid value for 'bins'")
  )
  
  bins = paste0('[', freq_table[['lower bound']], ', ', freq_table[['upper bound']], ')' )
  
  #determine counts of bins for each category
  freq_table_cat1 = tab_frequency_bins(scoresCat1, bins=data.frame(freq_table['lower bound'], freq_table['upper bound']))
  freq_table_cat2 = tab_frequency_bins(scoresCat2, bins=data.frame(freq_table['lower bound'], freq_table['upper bound']))
  
  
  X_ord = rep(paste0('[', freq_table_cat1[['lower bound']], ', ', freq_table_cat1[['upper bound']], ')' ), freq_table_cat1[['frequency']])
  Y_ord = rep(paste0('[', freq_table_cat2[['lower bound']], ', ', freq_table_cat2[['upper bound']], ')' ), freq_table_cat2[['frequency']])
  
  cats = c(rep(cat1, length(X_ord)), rep(cat2, length(Y_ord)))
  score = c(X_ord, Y_ord)
  
  dots <- list(...)
  if (!('xlbl' %in% names(dots))){
    dots$xlbl = deparse(substitute(catField))
  }
  if (!('ylbl' %in% names(dots))){
    dots$ylbl = deparse(substitute(scaleField))
  }
  
  do.call(vi_butterfly_chart, c(list(field=score, bin_field=cats, variation='butterfly'), dots))
  
}