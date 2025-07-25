% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meas_quantiles.R
\name{me_quantiles}
\alias{me_quantiles}
\title{Quantiles}
\usage{
me_quantiles(
  data,
  levels = NULL,
  k = 4,
  method = "own",
  indexMethod = c("sas1", "sas4", "hl", "excel", "hf8", "hf9"),
  q1Frac = c("linear", "down", "up", "bankers", "nearest", "halfdown", "midpoint"),
  q1Int = c("int", "midpoint"),
  q3Frac = c("linear", "down", "up", "bankers", "nearest", "halfdown", "midpoint"),
  q3Int = c("int", "midpoint")
)
}
\arguments{
\item{data}{list or dataframe}

\item{levels}{list, optional coding to use}

\item{k}{number of quantiles}

\item{method}{string, optional which method to use to calculate quartiles}

\item{indexMethod}{c("sas1", "sas4", "excel", "hl", "hf8", "hf9"), optional to indicate which type of indexing to use. Default is "sas1"}

\item{q1Frac}{c("linear", "down", "up", "bankers", "nearest", "halfdown", "midpoint"), optional to indicate what type of rounding to use for quantiles below 50 percent. Default is "linear"}

\item{q1Int}{c("int", "midpoint"), optional to indicate the use of the integer or the midpoint method for first quarter. Default is "int"}

\item{q3Frac}{c("linear", "down", "up", "bankers", "nearest", "halfdown", "midpoint"), optional to indicate what type of rounding to use for quantiles equal or above 50 percent. Default is "linear"}

\item{q3Int}{c("int", "midpoint"), optional  to indicate the use of the integer or the midpoint method for quantiles equal or above 50 percent. Default is "int"

method can be set to "own" and then provide the next parameters, or any of the methods listed in the notes.}
}
\value{
results : the quantiles, or if levels are used also additionally text versions
}
\description{
Quantiles split the data into k sections, each containing n/k scores. They can be seen as a generalisation of various 'tiles'. For example 4-quantiles is the same as the quartiles, 5-quantiles the same as quintiles, 100-quantiles the same as percentiles, etc.

Quite a few different methods exist to determine these. See the notes for more information.

This function is shown in this \href{https://youtu.be/119HkHrMu0M}{YouTube video} and the measure is also described at \href{https://peterstatistics.com/Terms/Measures/Quantiles.html}{PeterStatistics.com}
}
\details{
To determine the quartiles a specific indexing method can be used. See \strong{he_quantileIndexing()} for details on the different methods to choose from.

Then based on the indexes either linear interpolation or different rounding methods (bankers, nearest, down, up, half-down) can be used, or the midpoint between the two values. If the index is an integer either the integer or the mid point is used.

See the \strong{he_quantilesIndex()} for details on this.

Note that the rounding method can even vary per quantile, i.e. the one used for the ones below the median being different than the one those equal or above.

I've come across the following methods:\tabular{llllll}{
   method \tab indexing \tab q1 integer \tab q1 fractional \tab q3 integer \tab q3 fractional \cr
   sas1 \tab sas1 \tab use int \tab linear \tab use int \tab linear \cr
   sas2 \tab sas1 \tab use int \tab bankers \tab use int \tab bankers \cr
   sas3 \tab sas1 \tab use int \tab up \tab use int \tab up \cr
   sas5 \tab sas1 \tab midpoint \tab up \tab midpoint \tab up \cr
   hf3b \tab sas1 \tab use int \tab nearest \tab use int \tab halfdown \cr
   sas4 \tab sas4 \tab use int \tab linear \tab use int \tab linear \cr
   ms \tab sas4 \tab use int \tab nearest \tab use int \tab halfdown \cr
   lohninger \tab sas4 \tab use int \tab nearest \tab use int \tab nearest \cr
   hl2 \tab hl \tab use int \tab linear \tab use int \tab linear \cr
   hl1 \tab hl \tab use int \tab midpoint \tab use int \tab midpoint \cr
   excel \tab excel \tab use int \tab linear \tab use int \tab linear \cr
   pd2 \tab excel \tab use int \tab down \tab use int \tab down \cr
   pd3 \tab excel \tab use int \tab up \tab use int \tab up \cr
   pd4 \tab excel \tab use int \tab halfdown \tab use int \tab nearest \cr
   pd5 \tab excel \tab use int \tab midpoint \tab use int \tab midpoint \cr
   hf8 \tab hf8 \tab use int \tab linear \tab use int \tab linear \cr
   hf9 \tab hf9 \tab use int \tab linear \tab use int \tab linear \cr
}


The following values can be used for the \emph{method} parameter:
\enumerate{
\item sas1 = parzen = hf4 = interpolated_inverted_cdf = maple3 = r4. (Parzen, 1979, p. 108; SAS, 1990, p. 626; Hyndman & Fan, 1996, p. 363)
\item sas2 = hf3 = r3. (SAS, 1990, p. 626; Hyndman & Fan, 1996, p. 362)
\item sas3 = hf1 = inverted_cdf = maple1 = r1 (SAS, 1990, p. 626; Hyndman & Fan, 1996, p. 362)
\item sas4 = hf6 = minitab = snedecor = weibull = maple5 = r6 (Hyndman & Fan, 1996, p. 363; Weibull, 1939, p. ?; Snedecor, 1940, p. 43; SAS, 1990, p. 626)
\item sas5 = hf2 = CDF = averaged_inverted_cdf = r2 (SAS, 1990, p. 626; Hyndman & Fan, 1996, p. 362)
\item hf3b = closest_observation
\item ms (Mendenhall & Sincich, 1992, p. 35)
\item lohninger (Lohninger, n.d.)
\item hl1 (Hogg & Ledolter, 1992, p. 21)
\item hl2 = hf5 = Hazen = maple4 = r5 (Hogg & Ledolter, 1992, p. 21; Hazen, 1914, p. ?)
\item maple2
\item excel = hf7 = pd1 = linear = gumbel = maple6 = r7 (Hyndman & Fan, 1996, p. 363; Freund & Perles, 1987, p. 201; Gumbel, 1939, p. ?)
\item pd2 = lower
\item pd3 = higher
\item pd4 = nearest
\item pd5 = midpoint
\item hf8 = median_unbiased = maple7 = r8 (Hyndman & Fan, 1996, p. 363)
\item hf9 = normal_unbiased = maple8 = r9 (Hyndman & Fan, 1996, p. 363)
}

\emph{hf} is short for Hyndman and Fan who wrote an article showcasing many different methods, \emph{hl} is short for Hog and Ledolter, \emph{ms} is short for Mendenhall and Sincich, \emph{jf} is short for Joarder and Firozzaman. \emph{sas} refers to the software package SAS, \emph{maple} to Maple, \emph{pd} to Python's pandas library, and \emph{r} to R.

The names \emph{linear}, \emph{lower}, \emph{higher}, \emph{nearest} and \emph{midpoint} are all used by pandas quantile function and numpy percentile function. Numpy also uses \emph{inverted_cdf}, \emph{averaged_inverted_cdf}, \emph{closest_observation}, \emph{interpolated_inverted_cdf}, \emph{hazen}, \emph{weibull}, \emph{median_unbiased}, and \emph{normal_unbiased}.
}
\section{Before, After and Alternatives}{

Before this measure you might want an impression using a frequency table or a visualisation:
\code{\link{tab_frequency}}, for a frequency table
\code{\link{vi_bar_stacked_single}}, or Single Stacked Bar-Chart.
\code{\link{vi_bar_dual_axis}}, for Dual-Axis Bar Chart.

After this you might want some other descriptive measures:
\code{\link{me_consensus}}, for the Consensus.
\code{\link{me_hodges_lehmann_os}}, for the Hodges-Lehmann Estimate (One-Sample).
\code{\link{me_median}}, for the Median.
\code{\link{me_quartiles}}, for Quartiles / Hinges.
\code{\link{me_quartile_range}}, for Interquartile Range, Semi-Interquartile Range and Mid-Quartile Range.

or perform a test:
\code{\link{ts_sign_os}}, for One-Sample Sign Test.
\code{\link{ts_trinomial_os}}, for One-Sample Trinomial Test.
\code{\link{ts_wilcoxon_os}}, for One-Sample Wilcoxon Signed Rank Test.
}

\examples{
# Example 1: Dataframe
file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
studentDf = read.csv(file2, sep=';', na.strings=c("", "NA"))
ex1 = studentDf[['Teach_Motivate']]
order = c("Fully Disagree", "Disagree", "Neither disagree nor agree", "Agree", "Fully agree")
me_quantiles(ex1, levels=order)

#Example 2: Numeric data
ex2 = c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5)
me_quantiles(ex2)

#Example 3: Text data with
ex3 = c("a", "b", "f", "d", "e", "c")
order = c("a", "b", "c", "d", "e", "f")
me_quantiles(ex3, levels=order)

}
\references{
Freund, J. E., & Perles, B. M. (1987). A new look at quartiles of ungrouped data. \emph{The American Statistician, 41}(3), 200-203. https://doi.org/10.1080/00031305.1987.10475479

Galton, F. (1881). Report of the anthropometric committee. \emph{Report of the British Association for the Advancement of Science, 51}, 225-272.

Gumbel, E. J. (1939). La Probabilite des Hypotheses. \emph{Compes Rendus de l' Academie des Sciences, 209}, 645-647.

Hazen, A. (1914). Storage to be provided in impounding municipal water supply. \emph{Transactions of the American Society of Civil Engineers, 77}(1), 1539-1640. https://doi.org/10.1061/taceat.0002563

Hogg, R. V., & Ledolter, J. (1992). \emph{Applied statistics for engineers and physical scientists} (2nd int.). Macmillan.

Hyndman, R. J., & Fan, Y. (1996). Sample quantiles in statistical packages. \emph{The American Statistician, 50}(4), 361-365. https://doi.org/10.2307/2684934

Langford, E. (2006). Quartiles in elementary statistics. \emph{Journal of Statistics Education, 14}(3), 1-17. https://doi.org/10.1080/10691898.2006.11910589

Lohninger, H. (n.d.). Quartile. Fundamentals of Statistics. Retrieved April 7, 2023, from http://www.statistics4u.com/fundstat_eng/cc_quartile.html

McAlister, D. (1879). The law of the geometric mean. \emph{Proceedings of the Royal Society of London, 29}(196-199), 367-376. https://doi.org/10.1098/rspl.1879.0061

Mendenhall, W., & Sincich, T. (1992). \emph{Statistics for engineering and the sciences} (3rd ed.). Dellen Publishing Company.

Parzen, E. (1979). Nonparametric statistical data modeling. \emph{Journal of the American Statistical Association, 74}(365), 105-121. https://doi.org/10.1080/01621459.1979.10481621

SAS. (1990). SAS procedures guide: Version 6 (3rd ed.). SAS Institute.

Siegel, A. F., & Morgan, C. J. (1996). \emph{Statistics and data analysis: An introduction} (2nd ed.). J. Wiley.

Snedecor, G. W. (1940). \emph{Statistical methods applied to experiments in agriculture and biology} (3rd ed.). The Iowa State College Press.

Vining, G. G. (1998). \emph{Statistical methods for engineers}. Duxbury Press.

Weibull, W. (1939).* The phenomenon of rupture in solids*. Ingeniors Vetenskaps Akademien, 153, 1-55.
}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}
}
