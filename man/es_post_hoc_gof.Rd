% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eff_size_post_hoc_gof.R
\name{es_post_hoc_gof}
\alias{es_post_hoc_gof}
\title{Effect Sizes for a Goodness-of-Fit Post-Hoc Analysis}
\usage{
es_post_hoc_gof(post_hoc_results, es = "auto", bergsma = FALSE)
}
\arguments{
\item{post_hoc_results}{dataframe with the result of either ph_pairwise_bin(), ph_pairwise_gof(), ph_residual_bin(), or ph_residual_gof()}

\item{es}{string optional, the effect size to determine. Either 'auto', 'coheng', 'cohenh', 'ar', 'cramerv', 'cohenw', 'jbme', 'fei', 'rosenthal'}

\item{bergsma}{optional boolean. Use of Bergsma correction, only for Cramer V}
}
\value{
a dataframe with for residual post-hoc:

\describe{
\item{category}{, the label of the category}
\item{name effect size}{, the effect size value}
}

for pairwise post-hoc
\describe{
\item{category 1}{, the label of the first category}
\item{category 2}{, the label of the second category}
\item{name effect size}{, the effect size value}
}
}
\description{
Determines an effect size for each test (row) from the results of ph_pairwise_bin(), ph_pairwise_gof(), ph_residual_bin(), or ph_residual_gof().

This function is shown in this \href{https://youtu.be/Oc5Ejp-W2gU}{YouTube video} and described at \href{https://peterstatistics.com/Terms/Tests/PostHocAfterGoF.html}{PeterStatistics.com}
}
\details{
'auto' will use Cohen h for exact tests, Rosenthal correlation for z-tests and Cramer's V otherwise.

Cohen g ('coheng'), Cohen h ('cohenh') and Alternative Ratio ('ar') can all be used for any test.

Cramer V ('cramerv'), Cohen w ('cohenw'), Johnston-Berry-Mielke E ('jbme'), and Fei ('fei') can be used with chi-square tests (or likelihood ratio tests)

The Rosenthal Correlation ('rosenthal') can be used with a z-test (proportion/Wald/score/residual).

See the separate functions for each of these for details on the calculations.
}
\section{Before, After and Alternatives}{

Before this a post-hoc test might be helpful:
\code{\link{ph_pairwise_gof}}, for Pairwise Goodness-of-Fit Tests.
\code{\link{ph_pairwise_bin}}, for Pairwise Binary Test.
\code{\link{ph_residual_gof_gof}}, for Residuals Tests using Binary tests.
\code{\link{ph_residual_gof_bin}}, for Residuals Using Goodness-of-Fit Tests.

After this you might want to use a rule-of-thumb for the interpretation:
\code{\link{th_post_hoc_gof}}, for various rules-of-thumb.

Effect size in this function:
\code{\link{es_cohen_g}}, for Cohen g.
\code{\link{es_cohen_h_os}}, for Cohen h'.
\code{\link{es_alt_ratio}}, for Alternative Ratio.
\code{\link{es_cramer_v_gof}}, for Cramer's V for Goodness-of-Fit.
\code{\link{es_cohen_w}}, for Cohen's w.
\code{\link{es_jbm_e}}, for Johnston-Berry-Mielke E.
\code{\link{es_fei}}, for Fei.
\code{\link{r_rosenthal}}, for Rosenthal Correlation if a z-value is available.

note: the effect size functions are not used themselves in this function, but the same formulas are used.
}

\examples{
# Get data
dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
gssDf <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
ex1 = gssDf['mar1']

# Perform a post-hoc test
post_hoc_test = ph_pairwise_bin(ex1, test='binomial')

# Determine the effect sizes
es_post_hoc_gof(post_hoc_test, es='cohenh')

}
\author{
P. Stikker. \href{https://PeterStatistics.com}{Companion Website}, \href{https://www.youtube.com/stikpet}{YouTube Channel}, \href{https://www.patreon.com/bePatron?u=19398076}{Patreon donations}
}
