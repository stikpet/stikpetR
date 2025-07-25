#' Measures of Qualitative Variation
#' 
#' @description
#' The mode is the measure of central tendancy, to indicate the center for categorical data. Similar as the arithmetic mean is for numeric data. As with numeric data, the center alone is not always so informative. If your head is in a burning oven, and your feet are in a freezer, you are on average fine.
#' 
#' This is one of the reasons, why it is often recommended to add a measure of dispersion. It gives a clearer picture of the data, and can indicate how diverse it was (how much variation).
#' 
#' For categorical data there are a lot of different measures proposed, but I don't often see them being used. The most common one is probably the Variation Ratio. This is simply the percentage of cases that were not in the modal category.
#' 
#' The specific name of the type of measure for this qualitative variation can vary quite a lot. Some talk about dominance, differentiation, evenness, entropy, equitability, diversity, and apportionment.
#' 
#' I've tried to categorise the measures a bit, based on the calculations. Below is the overview of all measures available in this function.
#' 
#' |nr.|group|measure|source|original type|
#' |---|-----|-------|------|-------------|
#' |1|mode|Freeman Variation Ratio|(Freeman, 1965)| |
#' |2|mode|Berger-Parker Index|(Berger & Parker, 1970, p. 1345)|dominance|
#' |3|mode|Wilcox MODVR|(Wilcox, 1973, p. 7)| |
#' |4|mode|Wilcox RANVR|(Wilcox, 1973, p. 8)| |
#' |5|mean|Wilcox AVDEV|(Wilcox, 1973, p. 9)| |
#' |6|mean|Gibbs-Poston M4|(Gibbs & Poston, 1975, p. 473)|differentiation|
#' |7|mean|Gibbs-Poston M5|(Gibbs & Poston, 1975, p. 474)|differentiation|
#' |8|mean|Gibbs-Poston M6|(Gibbs & Poston, 1975, p. 474)|differentiation|
#' |9|mean|Wilcox VARNC = |(Wilcox, 1973, p. 11)| |
#' |9|mean|Gibbs-Poston M2 = |(Gibbs & Poston, 1975, p. 472)|differentiation|
#' |9|mean|Smith-Wilson E1*|(Smith & Wilson, 1996, p. 71)|evenness|
#' |10|mean|Wilcox STDEV|(Wilcox, 1973, p. 14)| |
#' |11|entropy|Shannon-Weaver Entropy|(Shannon & Weaver, 1949, p. 20)|entropy|
#' |12|entropy|Renyi Entropy|(Renyi, 1961, p. 549)|entropy|
#' |13|entropy|Wilcox HREL = |(Wilcox, 1973, p. 16)| |
#' |13|entropy|Pielou J|(Pielou, 1966, p. 141)|diversity|
#' |14|entropy|Sheldon Index|(Sheldon, 1969, p. 467)|equitability = relative diversity|
#' |15|entropy|Heip Evenness|(Heip, 1974, p. 555)|evenness|
#' |16|evenness|Hill Diversity|(Hill, 1973, p. 428)|diversity|
#' |17|evenness|Hill Evenness|(Hill, 1973, p. 429)|evenness|
#' |18|evenness|Bulla E|(Bulla, 1994, pp. 168-169)|evenness|
#' |19|evenness|Bulla D|(Bulla, 1994, p. 169)|diversity|
#' |20a|evenness|Simpson D|(Simpson, 1949, p. 688)|diversity|
#' |20b|evenness|Simpson D biased|(Smith & Wilson, 1996, p. 71)| |
#' |20c|evenness|Simpson D as diversity|(Wikipedia, n.d.)| |
#' |20d|evenness|Simpson D as diversity biased =|(Berger & Parker, 1970, p. 1345)| |
#' |20d|evenness|Gibbs-Poston M1|(Gibbs & Poston, 1975, p. 471)|differentiation|
#' |21|evenness|Gibbs-Poston M3|(Gibbs & Poston, 1975, p. 472)|differentiation|
#' |22|evenness|Smith-Wilson E2|(Smith & Wilson, 1996, p. 71)|evenness|
#' |23|evenness|Smith-Wilson E3|(Smith & Wilson, 1996, p. 71)|evenness|
#' |24|evenness|Fisher alpha|(Fisher et al., 1943, p. 55)|diversity|
#' |25|other|Wilcox MNDIF|(Wilcox, 1973, p. 9)| |
#' |26|other|Kaiser b|(Kaiser, 1968, p. 211)|apportionment|
#' 
#' \* Smith-Wilson E1 is listed with the mean group, since it uses the average frequency. It could of course also be placed in the evenness group.
#' 
#' This function is shown in this [YouTube video](https://youtu.be/9DDGMa0m4t8) and the measures are also described at [PeterStatistics.com](https://peterstatistics.com/Terms/Measures/QualitativeVariation.html)
#'  
#' @param data list or dataframe
#' @param measure optional to indicate which method to use. Either "vr" (default), "modvr", "ranvr", "avdev", "mndif", "varnc", "stdev", "hrel", "b", "m1", "m2", "m3", "m4", "m5", "m6", "d1", "d2", "d3", "d4", "bpi", "hd", "he", "swe", "re", "sw1", "sw2", "sw3", "hi", "si", "j", "b", "be", "bd", "fisher"
#' @param var1 optional additional value for some measures
#' @param var2 optional additional value for some measures
#' 
#' @returns
#' Dataframe with
#' \item{value}{the value of the requested measure}
#' \item{measure}{description of the measure calculated}
#' \item{source}{source used for calculation}
#' 
#' @details
#' The following measures can be determined:
#' 
#' \itemize{
#' \item *"modvr"*, Wilcox MODVR
#' \item *"ranvr"*, Wilcox RANVR
#' \item *"avdev"*, Wilcox AVDEV
#' \item *"mndif"*, Wilcox MNDIF
#' \item *"varnc"*, Wilcox VARNC (equal to Gibbs-Poston M2 and Smith-Wilson E1)
#' \item *"stdev"*, Wilcox STDEV
#' \item *"hrel"*, Wilcox HREL (equal to Pielou J)
#' \item *"m1"*, Gibbs-Poston M1
#' \item *"m2"*, Gibbs-Poston M2 (equal to Wilcox VARNC and Smith-Wilson E1)
#' \item *"m3"*, Gibbs-Poston M3
#' \item *"m4"*, Gibbs-Poston M4
#' \item *"m5"*, Gibbs-Poston M5
#' \item *"m6"*, Gibbs-Poston M6
#' \item *"b"*, Kaiser b
#' \item *"bd"*, Bulla D
#' \item *"be"*, Bulla E
#' \item *"bpi"*, Berger-Parker index
#' \item *"d1"*, *"d2"*, *"d3"*, *"d4"*, Simpson D and variations
#' \item *"hd"*, Hill Diversity, requires a value for *var1*
#' \item *"he"*, Hill Eveness, requires a value for *var1* and *var2*
#' \item *"hi"*, Heip Index
#' \item *"j"*, Pielou J (equal to Wilcox HREL)
#' \item *"si"*, Sheldon Index
#' \item *"sw1"*, Smith & Wilson E1 (equal to Wilcox VARNC and Gibbs-Poston M2)
#' \item *"sw2"*, Smith-Wilson E2
#' \item *"sw3"*, Smith-Wilson E3
#' \item *"swe"*, Shannon-Weaver Entropy
#' \item *"re"*, Renyi entropy, requires a value for *var1*
#' \item *"vr"*, Freeman's variation ratio
#' \item *fisher*, Fisher alpha
#' }
#' 
#' **MODE BASED MEASURES**
#' 
#' Dispersion can be seen as how much variation there is, using as a norm the center. For nominal data the measure of central tendancy is the mode, and therefor some measures of qualitative variation use the mode as the starting point.
#' 
#' The frequency of the modal category is then useful. This is simply the maximum of the frequencies.
#' 
#' **Freeman Variation Ratio** ("vr")
#' 
#' Perhaps one of the most popular measures of qualitative variation uses the mode. The (Freeman) Variation Ratio. It is simply the proportion of scores that do not belong to the modal category. In formula notation (Freeman, 1965, p. 41):
#' 
#' Formula used from Freeman (1965, p. 41):
#' \deqn{v = 1 - \frac{F_{mode}}{n}}
#' 
#' This variation ratio would become 0% if all cases fitted in the modal category, and all other categories don't have any cases.
#' 
#' A 0 (0%) would mean that all cases were in the modal category. A 1 (100%) would indicate that no cases were in the modal category. However, this seems impossible to ever occur, since the modal category is the category with the highest frequency, which is impossible to be 0, unless there are no cases at all.
#' 
#' **Berger-Parker index** ("bpi")
#' 
#' The variation ratio is the opposite of the Berger-Parker Index, which is simply the proportion of scores that did fit in the modal category. In formula notation (Berger & Parker, 1970, p. 1345):
#' \deqn{BPI = \frac{F_{mode}}{n}}
#' 
#' Berger and Parker refer to this as a dominance measure, to indicate how "dominant" the modal category is.
#' 
#' A 1 (100%) would mean that all cases were in the modal category. A 0 (0%) would indicate that no cases were in the modal category. However, this seems impossible to ever occur, since the modal category is the category with the highest frequency, which is impossible to be 0, unless there are no cases at all.
#' 
#' **Wilcox MODVR** ("modvr")
#' 
#' This looks at the difference of the frequency for each category with the modal frequency. This then gets divided by \eqn{n\times \left(k -1\right)} to standardize the results to 0 to 1.
#' 
#' It is a modification of the Freeman Variation Ratio, hence the name MODVR. Wilcox noted that the Freeman VR can never reach the maximum value of 1.
#' 
#' The formula used is (Wilcox, 1973, p. 7):
#' \deqn{\text{MODVR} = \frac{\sum_{i=1}^k F_{mode} - F_i}{n\times \left(k - 1\right)} = \frac{k\times F_{mode}-n}{n\times \left(k - 1\right)}}
#' 
#' **Wilcox RANVR** ("ranvr")
#' 
#' Short for 'range variation ratio' this measure is very similar to Freeman's VR. Instead of looking simply at the mode, it looks at the range.
#' 
#' The formula used is (Wilcox, 1973, p. 8):
#' \deqn{\text{RANVR} = 1 - \frac{F_{mode} - F_{min}}{F_{mode}}}
#' 
#' **MEAN BASED MEASURES**
#' 
#' The following measures use the average count to determine the variation. i.e.
#' \deqn{\bar{F} = \frac{\sum_{i=1}^k F_i}{k} = \frac{n}{k}}
#' 
#' **Wilcox AVDEV** ("avdev")
#' 
#' This simply follows the mean absolute deviation analogue but then using frequencies. Again this is then standardized.
#' 
#' The formula used is (Wilcox, 1973, p. 9):
#' \deqn{\text{AVDEV} = 1-\frac{\sum_{i=1}^k \left|F_i-\bar{F}\right|}{2\times \frac{n}{k}\times \left(k-1\right)}= 1-\frac{k\times \sum_{i=1}^k \left|F_i-\bar{F}\right|}{2\times n \times \left(k-1\right)}} 
#' 
#' **Gibbs-Poston M4** ("m4")
#' 
#' The formula used (Gibbs & Poston, 1975, p. 473):
#' \deqn{\text{M4} = 1-\frac{\sum_{i=1}^k \left|F_i-\bar{F}\right|}{2\times n}}
#' 
#' **Gibbs-Poston M5** ("m5")
#' 
#' The problem with M4 is that it can never be 0, so to adjust for this M5 could be used but is computationally then more difficult.
#' 
#' The formula used (Gibbs & Poston, 1975, p. 474):
#' \deqn{\text{M5} = 1-\frac{\sum_{i=1}^k \left|F_i-\bar{F}\right|}{2\times\left(n-k+1-\bar{F}\right)}}
#' 
#' **Gibbs-Poston M6** ("m6")
#' 
#' The formula used (Gibbs & Poston, 1975, p. 474):
#' \deqn{\text{M6} = k\times\left(1-\frac{\sum_{i=1}^k \left|F_i-\bar{F}\right|}{2\times n}\right) = k\times\text{M4}}
#' 
#' **Wilcox VARNC** ("varnc"), **Gibbs-Poston M2** ("m2"), and **Smith & Wilson E1** ("sw1")
#' 
#' This is similar as the variance for scale variables.
#' 
#' The formula used is (Wilcox, 1973, p. 11):
#' \deqn{\text{VARNC} = 1-\frac{\sum_{i=1}^{k}\left(F_i-\bar{F}\right)^2}{\frac{n^2\times\left(k-1\right)}{k}} = \frac{k\times\left(n^2-\sum_{i=1}^k F_i^2\right)}{n^2\times\left(k-1\right)}}
#' 
#' This is the same as Gibbs and Poston's **M2** ("m2"). Their formula looks different but has the same result (Gibbs & Poston, 1975, p. 472)
#' \deqn{\text{M2} = \frac{1-\sum_{i=1}^k p_i^2}{1-\frac{1}{k}} = \frac{\text{M1}}{1-\frac{1}{k}} = \frac{k}{k-1}\times\text{M1}}
#' 
#' It is also the same as Smith and Wilson's first evenness measure ("sw1").
#' 
#' The formula used (Smith & Wilson, 1996, p. 71):
#' \deqn{E_1 = \frac{1 - D_s}{1 - \frac{1}{k}}}
#' 
#' With \eqn{D_s} being Simpson's D, but defined as:
#' \deqn{D_s = \sum_{i=1}^k\left(\frac{F_i}{n}\right)^2}
#' 
#' **Wilcox STDEV** ("stdev")
#' 
#' As with the variance for scale variables, we can take the square root to obtain the standard deviation.
#' 
#' The formula used can be from the VARNC or the MNDIF (Wilcox, 1973, p. 14):
#' \deqn{\text{STDEV} = 1-\sqrt{\frac{\sum_{i=1}^k \left(F_i-\bar{F}\right)^2}{\left(n-\bar{F}\right)^2+\left(k-1\right)\bar{F}^2}}= 1-\sqrt{\frac{\sum_{i=1}^{k-1}\sum_{j=i+1}^k \left(F_i-F_j\right)^2}{n^2\times\left(k-1\right)}}}
#' 
#' **ENTROPY**
#' 
#' Entropy is sometimes referred to as the expected value of the surprise. It tells on average how surprised we might be about the outcome, and is also used as a measure with qualitative data.
#' 
#' I enjoyed the simple explanation on entropy from StatQuest, their video is available [here](https://www.youtube.com/watch?v=YtebGVx-Fxw).
#' 
#' It deals a lot with proportions rather than the counts themselves
#' 
#' **Shannon-Weaver Entropy** ("swe")
#' 
#' The formula used (Shannon & Weaver, 1949, p. 20):
#' \deqn{H_{sw}=-\sum_{i=1}^k p_i\times\ln\left(p_i\right)}
#' 
#' **Renyi entropy** ("re")
#' 
#' This is a generalisation for Shannon entropy.
#' 
#' The formula used is (Renyi, 1961, p. 549):
#' \deqn{H_q = \frac{1}{1 - q}\times\log_2\left(\sum_{i=1}^k p_i^q\right)}
#' 
#' **Wilcox HREL** ("hrel") and **Pielou J** ("j")
#' 
#' This uses Shannon's entropy but divides it over the maximum possible uncertainty.
#' 
#' The formula used (Wilcox, 1973, p. 16):
#' \deqn{\text{HREL} = \frac{-\sum_{i=1}^k p_i \times \text{log}_2 p_i}{\text{log}_2 k}}
#' 
#' This is the same as Pielou J. ("j")
#' 
#' The formula used (Pielou, 1966, p. 141):
#' \deqn{J=\frac{H_{sw}}{\ln\left(k\right)}}
#' 
#' **Sheldon Index** ("si")
#' 
#' The formula used (Sheldon, 1969, p. 467):
#' \deqn{E = \frac{e^{H_{sw}}}{k}}    
#' 
#' **Heip Index** ("hi")
#' 
#' The formula used is (Heip, 1974, p. 555):
#' \deqn{E_h = \frac{e^{H_{sw}} - 1}{k - 1}} 
#' 
#' **EVENNESS and DIVERSITY**
#' 
#' **Hill Diversity** ("hd")
#' 
#' The formula used is (Hill, 1973, p. 428):
#' \deqn{N_a = \begin{cases}\left(\sum_{i=1}^k p_i^a\right)^{\frac{1}{1-a}} & \text{ if } a\neq 1 \\ e^{H_{sw}} & \text{ if }  =1 \end{cases}}
#' 
#' **Hill Eveness** ("he")
#' 
#' The formula used is (Hill, 1973, p. 429):
#' \deqn{E_{a,b} = \frac{N_a}{N_b}}
#' 
#' Where \eqn{N_a} and \eqn{N_b} are Hill's diversity values for a and b.
#' 
#' **Bulla E** ("be")
#' 
#' Bulla's evenness measure.
#' 
#' The formula used is (Bulla, 1994, pp. 168-169):
#' \deqn{E_b = \frac{O - \frac{1}{k} - \frac{k - 1}{n}}{1 - \frac{1}{k} - \frac{k - 1}{n}}}
#' 
#' With:
#' \deqn{O = \sum_{i=1}^k \min\left(p_i, \frac{1}{k}\right)}
#' 
#' **Bulla D** ("bd")
#' 
#' Bulla's Evenness measure converted to a diversity measure.
#' 
#' The formula used is (Bulla, 1994, p. 169):
#' \deqn{D_b = E_b\times k}
#' 
#' Where \eqn{E_b} is Bulla E value.
#' 
#' With:
#' \deqn{O = \sum_{i=1}^k \min\left(p_i, \frac{1}{k}\right)}
#' 
#' **Simpson D** ("d1", "d2", "d3", "d4" = Gibbs-Poston M1)
#' 
#' The formula used is based on Simpson (1949, p. 688):
#' \deqn{D_1 = \frac{\sum_{i=1}^k F_i\times\left(F_i-1\right)}{n\times\left(n-1\right)}}
#' 
#' Another alternative is for a population:
#' \deqn{D_2 = \sum_{i=1}^k\left(\frac{F_i}{n}\right)^2}
#' 
#' Often the result is subtracted from 1 to reverse the scale. 
#' 
#' \deqn{D_3 = 1-\frac{\sum_{i=1}^k F_i\times\left(F_i-1\right)}{n\times\left(n-1\right)}}
#' 
#' and
#' \deqn{D_4 = 1 - \sum_{i=1}^k\left(\frac{F_i}{n}\right)^2}
#' 
#' This last one is then the same as Gibb-Poston M1 (Gibbs & Poston, 1975, p. 471):
#' \deqn{\text{M1} = 1 - \sum_{i=1}^k p_i^2}
#' 
#' **Gibbs-Poston M3** ("m3")
#' 
#' The formula used (Gibbs & Poston, 1975, p. 472):
#' \deqn{\text{M3} = \frac{1-\sum_{i=1}^k p_i^2-p_{min}}{1-\frac{1}{k}-p_{min}}}
#' 
#' With \eqn{p_{min}} the lowest proportion
#' 
#' **Smith & Wilson E2** ("sw2")
#' 
#' The formula used (Smith & Wilson, 1996, p. 71):
#' \deqn{E_2 = \frac{\ln\left(D_s\right)}{\ln\left(k\right)}}
#' 
#' With \eqn{D_s} being Simpson's D, but defined as:
#' \deqn{D_s = \sum_{i=1}^k\left(\frac{F_i}{n}\right)^2}
#' 
#' **Smith & Wilson E3** ("sw3")
#' 
#' The formula used (Smith & Wilson, 1996, p. 71):
#' \deqn{E_3 = \frac{1}{D_s \times k}}
#' 
#' With \eqn{D_s} being Simpson's D, but defined as:
#' \deqn{D_s = \sum_{i=1}^k\left(\frac{F_i}{n}\right)^2}
#' 
#' **Fisher alpha** ("fisher")
#' 
#' The formula used (Fisher et al., 1943, p. 55):
#' \deqn{k = \alpha \times \ln\left(1 + \frac{n}{\alpha}\right)}
#' 
#' The function uses a simple binary search to find the value for \eqn{\alpha} such that the result of the above formula will produce the number of categories (\eqn{k}).
#' 
#' **OTHER***
#' 
#' **Wilcox MNDIF** ("mndif")
#' 
#' Analog of the mean difference measure for scale variables.
#' 
#' The formula used is (Wilcox, 1973, p. 9):
#' \deqn{\text{MNDIF} = 1-\frac{\sum_{i=1}^{k-1}\sum_{j=i+1}^k \left|F_i-F_j\right|}{n\times\left(k-1\right)}}
#' 
#' **Kaiser b**
#' 
#' The formula used (Kaiser, 1968, p. 211):
#' \deqn{B = 1 - \sqrt{1 - \left(\sqrt[k]{\prod_{i=1}^k\frac{f_i\times k}{n}}\right)^2}}
#' 
#' Kaiser also provides rules-of-thumb for interpretation (see \code{\link{th_kaiser_b}}, for these).
#' 
#' @section Before, After and Alternatives:
#' Before this an impression using a frequency table or a visualisation might be helpful:
#' \code{\link{tab_frequency}}, for a frequency table. 
#' \code{\link{vi_bar_simple}}, for Simple Bar Chart. 
#' \code{\link{vi_cleveland_dot_plot}}, for Cleveland Dot Plot.
#' \code{\link{vi_dot_plot}}, for Dot Plot.
#' \code{\link{vi_pareto_chart}}, for for Pareto Chart.
#' \code{\link{vi_pie}}, for Pie Chart.
#' 
#' After this you might want to perform a test:
#' \code{\link{ts_pearson_gof}}, for Pearson Chi-Square Goodness-of-Fit Test. 
#' \code{\link{ts_freeman_tukey_gof}}, for Freeman-Tukey Test of Goodness-of-Fit. 
#' \code{\link{ts_freeman_tukey_read}}, for Freeman-Tukey-Read Test of Goodness-of-Fit.
#' \code{\link{ts_g_gof}}, for G (Likelihood Ratio) Goodness-of-Fit Test. 
#' \code{\link{ts_mod_log_likelihood_gof}}, for Mod-Log Likelihood Test of Goodness-of-Fit. 
#' \code{\link{ts_multinomial_gof}}, for Multinomial Goodness-of-Fit Test. 
#' \code{\link{ts_neyman_gof}}, for Neyman Test of Goodness-of-Fit. 
#' \code{\link{ts_powerdivergence_gof}}, for Power Divergence GoF Test. 
#' 
#' @references
#' Berger, W. H., & Parker, F. L. (1970). Diversity of planktonic foraminifera in deep-sea sediments. *Science, 168*(3937), 1345-1347. doi:10.1126/science.168.3937.1345
#' 
#' Bulla, L. (1994). An index of evenness and its associated diversity measure. *Oikos, 70*(1), 167-171. doi:10.2307/3545713
#' 
#' Fisher, R. A., Corbet, A. S., & Williams, C. B. (1943). The relation between the number of species and the number of individuals in a random sample of an animal population. *The Journal of Animal Ecology, 12*(1), 42-58. doi:10.2307/1411
#' 
#' Freeman, L. C. (1965). *Elementary applied statistics: For students in behavioral science*. Wiley.
#' 
#' Gibbs, J. P., & Poston, D. L. (1975). The division of labor: Conceptualization and related measures. *Social Forces, 53*(3), 468. doi:10.2307/2576589
#' 
#' Heip, C. (1974). A new index measuring evenness. *Journal of the Marine Biological Association of the United Kingdom, 54*(3), 555-557. doi:10.1017/S0025315400022736
#' 
#' Hill, M. O. (1973). Diversity and evenness: A unifying notation and its consequences. *Ecology, 54*(2), 427-432. doi:10.2307/1934352
#' 
#' Kaiser, H. F. (1968). A measure of the population quality of legislative apportionment. *American Political Science Review, 62*(1), 208-215. doi:10.2307/1953335
#' 
#' Pielou, E. C. (1966). The measurement of diversity in different types of biological collections. *Journal of Theoretical Biology, 13*, 131-144. doi:10.1016/0022-5193(66)90013-0
#' 
#' Renyi, A. (1961). On measures of entropy and information. *Contributions to the Theory of Statistics, 1*, 547-562.
#' 
#' Shannon, C. E., & Weaver, W. (1949). *The mathematical theory of communication*. The university of Illinois press.
#' 
#' Sheldon, A. L. (1969). Equitability indices: Dependence on the species count. *Ecology, 50*(3), 466-467. doi:10.2307/1933900
#' 
#' Simpson, E. H. (1949). Measurement of diversity. *Nature, 163*(4148), Article 4148. doi:10.1038/163688a0
#' 
#' Smith, B., & Wilson, J. B. (1996). A consumer's guide to evenness indices. *Oikos, 76*(1), 70-82. doi:10.2307/3545749
#' 
#' Wilcox, A. R. (1973). Indices of qualitative variation and political measurement. *Political Research Quarterly, 26*(2), 325-343. doi:10.1177/106591297302600209
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @examples 
#' #Example 1: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['mar1']
#' me_qv(ex1)
#' 
#' #Example 2: a list
#' ex2 = c("MARRIED", "DIVORCED", "MARRIED", "SEPARATED", "DIVORCED", "NEVER MARRIED", "DIVORCED", 
#' "DIVORCED", "NEVER MARRIED", "MARRIED", "MARRIED", "MARRIED", "SEPARATED", "DIVORCED", 
#' "NEVER MARRIED", "NEVER MARRIED", "DIVORCED", "DIVORCED", "MARRIED")
#' me_qv(ex2, "swe")
#' 
#' @export 
me_qv <- function(data, measure="vr", var1=2, var2=1){
  
  freqs = table(data)
  
  k = length(freqs)
  n = sum(freqs)
  fm = max(freqs)
  props = freqs/n
  
  if (measure=="modvr"){
    #Modified Variation Ratio
    src = "(Wilcox, 1973, p. 7)"
    lbl = "Wilcox MODVR"
    qv = sum(fm - freqs)/(n*(k - 1))}
  else if (measure=="ranvr"){
    #Range Variation Ratio
    src = "(Wilcox, 1973, p. 8)"
    lbl = "Wilcox RANVR"
    fl = min(freqs)
    qv = 1 - (fm - fl)/fm}
  else if (measure=="avdev"){
    #Average Deviation
    src = "(Wilcox, 1973, p. 9)"
    lbl = "Wilcox AVDEV"
    qv = 1-sum(abs(freqs-n/k)) / (2*n/k*(k-1))}
  else if (measure=="mndif"){
    #MNDif
    src = "(Wilcox, 1973, p. 9)"
    lbl = "Wilcox MNDIF"
    mndif = 0
    for (i in 1:(k-1)){
      for (j in (i+1):k){
        mndif = mndif + abs(freqs[i]-freqs[j])}
    }
    qv = 1 - mndif/(n*(k-1))}
  else if (measure=="varnc"){
    #VarNC
    src = "(Wilcox, 1973, p. 11)"
    lbl = "Wilcox VARNC"
    qv = 1 - sum((freqs-n/k)**2)/(n**2*(k-1)/k)}
  else if (measure=="stdev"){
    src = "(Wilcox, 1973, p. 14)"
    lbl = "Wilcox STDEV"
    qv = 1 - (sum((freqs-n/k)**2)/((n-n/k)**2+(k-1)*(n/k)**2))**0.5}
  else if (measure=="hrel"){
    #HRel
    src = "(Wilcox, 1973, p. 16)"
    lbl = "Wilcox HREL"
    hrel = 0
    for (i in 1:k){
      hrel = hrel + props[i]*log2(props[i])}
    qv = -hrel/log2(k)}  
  else if (measure=="m1"){
    src = "(Gibbs & Poston, 1975, p. 471)"
    lbl = "Gibbs-Poston M1"
    qv = 1 - sum(props**2)}
  else if (measure=="m2"){
    #equal to varnc
    src = "(Gibbs & Poston, 1975, p. 472)"
    lbl = "Gibbs-Poston M2"
    qv = (1 - sum(props**2)) / (1-1/k)}
  else if (measure=="m3"){
    src = "(Gibbs & Poston, 1975, p. 472)"
    lbl = "Gibbs-Poston M3"
    pl = min(props)
    qv = (1 - sum(props**2) - pl) / (1-1/k - pl)}
  else if (measure=="m4"){
    src = "(Gibbs & Poston, 1975, p. 473)"
    lbl = "Gibbs-Poston M4"
    fmean = n/k
    qv = 1 - sum(abs(freqs-fmean))/(2*n)}
  else if (measure=="m5"){
    src = "(Gibbs & Poston, 1975, p. 474)"
    lbl = "Gibbs-Poston M5"
    fmean = n/k
    qv = 1 - sum(abs(freqs-fmean))/(2*(n-k+1-fmean))}
  else if (measure=="m6"){
    src = "(Gibbs & Poston, 1975, p. 474)"
    lbl = "Gibbs-Poston M6"
    fmean = n/k
    qv = k*(1 - sum(abs(freqs-fmean))/(2*n))}
  else if (measure=="b"){
    #Kaiser B index
    src = "(Kaiser, 1968, p. 211)"
    lbl = "Kaiser b"
    qv = 1 - (1 - ((prod(freqs*k/n))**(1/k))**2)**0.5} 
  else if (measure=="bd"){
    #Bulla D
    src = "(Bulla, 1994, p. 169)"
    lbl = "Bulla D"
    o = 0
    for (p in props){
      o = o + min(p, 1/k)}
    qv = k*(o - 1/k + (k - 1)/n)/(1 - 1/k + (k-1)/n)}
  else if (measure=="be"){
    #Bulla e
    src = "(Bulla, 1994, pp. 168-169)"
    lbl = "Bulla E"
    o = 0
    for (p in props){
      o = o + min(p, 1/k)}
    qv = (o - 1/k + (k - 1)/n)/(1 - 1/k + (k-1)/n)}
  else if (measure=="bpi"){
    #Berger-Parker Index
    src = "(Berger & Parker, 1970, p. 1345)"
    lbl = "Berger-Parker D"
    qv = fm/n}
  else if (measure=="d1"){
    #Simpson's D
    src = "(Simpson, 1949, p. 688)"
    lbl = "Simpson D"
    qv = sum(freqs*(freqs-1))/(n*(n-1))}
  else if (measure=="d2"){
    #Simpson's D
    src = "(Smith & Wilson, 1996, p. 71)"
    lbl = "Simpson D biased"
    qv = sum((freqs/n)**2)}
  else if (measure=="d3"){
    #Simpson's D
    src = "(Wikipedia, n.d.)"
    lbl = "Simpson D as diversity"
    qv = 1 - sum(freqs*(freqs-1))/(n*(n-1))}
  else if (measure=="d4"){
    #Simpson's D
    src = "(Berger & Parker, 1970, p. 1345)"
    lbl = "Simpson D as diversity biased"
    qv = 1 - sum((freqs/n)**2)}
  else if (measure=="hd"){
    #Hill's Diversity
    src = "(Hill, 1973, p. 428)"
    lbl = "Hill Diversity"
    if (var1 == 1){
      qv = exp(-1*sum(props*log(props)))}
    else {
      qv = (sum(props**var1)**(1/(1-var1)))}
  }
  else if (measure=="he"){
    #Hill's Evenness
    src = "(Hill, 1973, p. 429)"
    lbl = "Hill Evenness"
    qv = me_qv(data, measure="hd", var1=var1)['value']/me_qv(data, measure="hd", var1=var2)['value']}
  else if (measure=="hi"){
    #Heip Index
    src = "(Heip, 1974, p. 555)"
    lbl = "Heip Evenness"
    h = -1*sum(props*log(props))
    qv = (exp(h) - 1)/(k - 1)}
  else if (measure=="j"){
    #Pielou J
    src = "(Pielou, 1966, p. 141)"
    lbl = "Pielou J"
    h = -1*sum(props*log(props))
    qv = h/log(k)}
  else if (measure=="si"){
    #Sheldon Index
    src = "(Sheldon, 1969, p. 467)"
    lbl = "Sheldon Evenness"
    h = -1*sum(props*log(props))
    qv = exp(h)/k}
  else if (measure=="sw1"){
    #Smith and Wilson Index 1
    src = "(Smith & Wilson, 1996, p. 71)"
    lbl = "Smith-Wilson Evenness Index 1"
    d = sum(props**2)
    qv = (1 - d)/(1 - 1/k)}
  else if (measure=="sw2"){
    #Smith and Wilson Index 2
    src = "(Smith & Wilson, 1996, p. 71)"
    lbl = "Smith-Wilson Evenness Index 2"
    d = sum(props**2)
    qv = -log(d)/log(k)}
  else if (measure=="sw3"){
    #Smith and Wilson Index 3
    src = "(Smith & Wilson, 1996, p. 71)"
    lbl = "Smith-Wilson Evenness Index 3"
    d = sum(props**2)
    qv = 1/(d*k)}    
  else if (measure=="swe"){
    #Shannon-Weaver Entropy
    src = "(Shannon & Weaver, 1949, p. 20)"
    lbl = "Shannon-Weaver Entropy"
    qv = -1*sum(props*log(props))}
  else if (measure=="re"){
    #Renyi Entropy
    src = "(Renyi, 1961, p. 549)"
    lbl = "Reneyi Entropy"
    qv = 1/(1 - var1)*log2(sum(props**var1))}
  else if (measure=="vr"){
    #Variation Ratio
    src = "(Freeman, 1965)"
    lbl = "Freeman Variation Ratio"
    pm = fm/n
    qv = 1 - pm}
  else if (measure=="fisher"){
    src = "(Fisher et al., 1943, p. 55)"
    lbl = "Fisher alpha"
    maxIter = 100
    a1 = 1
    k1 = a1 * log(1 + n/a1)
    
    if (k1 != k){
      if (k1 > k){
        a2 = 0.5}
      else {
        a2 = 2}
      
      k2 = a2 * log(1 + n / a2)
      
      if (k2 != k){
        k3 = k2
        iters = 0
        
        while (iters < maxIter && k3 != k){
          iters = iters + 1
          
          if (k2 > k){
            if (k1 > k){
              a3 = a2 - abs(a2 - a1)}
            else {
              a3 = a2 - abs(a2 - a1) / 2}
          }
          
          else {
            if (k1 < k){
              a3 = a2 + abs(a2 - a1)}
            else {
              a3 = a2 + abs(a2 - a1) / 2}
          }
          
          if (a3 == 0){
            a3 = a2 - abs(a2 - a1) / 2}
          
          k3 = a3 * log(1 + n / a3)
          
          a1 = a2
          a2 = a3
          k1 = k2
          k2 = k3
        }
      }
      
      else {
        a3 = a2}
    }
    else {
      a3 = a1}
    qv = a3
    
  }
  
  results <- data.frame(unname(qv), lbl, src)
  colnames(results)<-c("value", "measure", "source")
  return (results)
}



