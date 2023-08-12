#' Effect Sizes for Binary vs. Binary
#' @description 
#' Various measures of association/agreement/similarity for binary vs. binary cases.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' @param method : the method to use. Default is odds-ratio
#' 
#' @returns 
#' the effect size measure value
#' 
#' @details
#' The method can be set to any of the following:
#' 
#' | | | | |
#' |----|----|----|----|
#' |alroy|ample|anderberg|austin-colwell|
#' |baroni-urbani-buser-1|baroni-urbani-buser-2|becker-clogg-1|becker-clogg-2|
#' |bonett-price-1|bonett-price-2|bonett-price-3|braun-blanquet|
#' |camp-1|camp-2|camp-3|chen-popovich|
#' |clement|cohen-kappa|cohen-w|cole-c1|
#' |cole-c2|cole-c3|cole-c4|cole-c5|
#' |cole-c6|cole-c7|cole-c8|contingency|
#' |czekanowski|dennis|dice-1|dice-2|
#' |dice-3|digby|doolittle|driver-kroeber-1|
#' |driver-kroeber-2|edward|eyraud|fager-mcgowan-1|
#' |fager-mcgowan-2|faith|fleiss|forbes-1|
#' |forbes-2|fossum-kaskey|gilbert|gilbert-wells|
#' |gk-lambda-1|gk-lambda-2|gleason|gower|
#' |gower-legendre|hamann|harris-lahey|hawkins-dotson|
#' |hurlbert|jaccard|johnson|kent-foster-1|
#' |kent-foster-2|kuder-richardson|kulczynski-1|kulczynski-2|
#' |loevinger|matching|maxwell-pilliner|mcconnaughey|
#' |mcewen-michael|mountford|nei-li|ochiai-1|
#' |ochiai-2|odds-ratio|otsuka|pearson|
#' |pearson-heron|pearson-q1|pearson-q2|pearson-q3|
#' |pearson-q4|pearson-q5|peirce-1|peirce-2|
#' |peirce-3|phi|rogers-tanimoto|rogot-goldberg|
#' |russell-rao|scott|simpson|sokal-michener|
#' |sokal-sneath-1|sokal-sneath-2|sokal-sneath-3|sokal-sneath-4|
#' |sokal-sneath-5|sorgenfrei|stiles|tanimoto|
#' |tarantula|tarwid|tulloss|yule-q|
#' |yule-r|yule-y|  |  |
#' 
#' 
#'  If we have a 2x2 table with the following values:
#' 
#'  |   |Column 1| Column 2| total|
#'  |---|:-------:|:--------:|:------:|
#'  |row 1| \eqn{a} | \eqn{b} | \eqn{R_1}|
#'  |row 2| \eqn{c} | \eqn{d} | \eqn{R_2}|
#'  |total| \eqn{C_1} | \eqn{C_2} | \eqn{n}|
#' 
#'  then the following formula are used:
#' 
#'  |Nr|Label|Formula|
#'  |--|-----------|---------|
#'  |1|Russell-Rao|\eqn{\frac{a}{n}}|
#'  |2|Dice-1|\eqn{\frac{a}{R_1}}|
#'  |3|Dice-2|\eqn{\frac{a}{C_1}}|
#'  |4|Braun-Blanquet|\eqn{\frac{a}{\max\left(R_1, C_1\right)}}|
#'  |5|Simpson Similarity|\eqn{\frac{a}{\min\left(R_1, C_1\right)}}|
#'  |6|Kulczynski-1|\eqn{\frac{a}{b+c}}|
#'  |7|Jaccard|\eqn{\frac{a}{a+b+c}}|
#'  |8|Sokal-Sneath-1|\eqn{\frac{a}{a+2b+2c}}|
#'  |9|Gleasson|\eqn{\frac{2a}{2a+b+c}}|
#'  |10|Mountford|\eqn{\frac{2a}{a\left(b+c\right)+2bc}}|
#'  |11|Driver-Kroeber|\eqn{\frac{a}{\sqrt{R_1 C_1}}}|
#'  |12|Sorgenfrei|\eqn{\frac{a^2}{R_1 C_1}}|
#'  |13|Johnson|\eqn{\frac{a}{R_1}+\frac{a}{C_1}}|
#'  |14|Kulczynski-2|\eqn{\frac{1}{2}\left(\frac{a}{R_1}+\frac{a}{C_1}\right)}|
#'  |15|Fager-McGowan-1|\eqn{\frac{a}{\sqrt{R_1 C_1}}-\frac{1}{2\sqrt{\max\left(R_1, C_1\right)}}}|
#'  |16|Fager-McGowan-2|\eqn{\frac{a}{\sqrt{R_1 C_1}}-\frac{\sqrt{\max\left(R_1, C_1\right)}}{2}}|
#'  |17|tarantula|\eqn{\frac{aR_2}{cR_1}}|
#'  |18|Ample|\eqn{\left| \frac{aR_2}{cR_1}\right|}|
#'  |19|Gilbert|\eqn{\frac{an-R_1 C_1}{C_1 n + R_1 n - an - R_1 C_1}}|
#'  |20|Fossum-Kaskey|\eqn{\frac{n\left(a - \frac{1}{2}\right)^2}{R_1 C_1}}|
#'  |21|Forbes - 1|\eqn{\frac{na}{R_1 C_1}}|
#'  |22|Eyraud|\eqn{\frac{a-R_1 C_1}{R_1 R_2 C_1 C_2}}|
#'  |23|Sokal-Michener|\eqn{\frac{a+d}{n}}|
#'  |24|Faith|\eqn{\frac{a+\frac{1}{2}}{n}}|
#'  |25|Sokal-Sneath-5|\eqn{\frac{a+d}{b+c}}|
#'  |26|Rogers-Tanimoto|\eqn{\frac{a+d}{a+2\left(b+c\right)+d}}|
#'  |27|Sokal-Sneath-2|\eqn{\frac{2a+2d}{2a+b+c+2d}}|
#'  |28|Gower|\eqn{\frac{a+d}{\sqrt{R_1 R_2 C_1 C_2}}}|
#'  |29|Sokal-Sneath-4|\eqn{\frac{ad}{\sqrt{R_1 R_2 C_1 C_2}}}|
#'  |30|Rogot-Goldberg|\eqn{\frac{a}{R_1 + C_1}+\frac{d}{R_2 + C_2}}|
#'  |31|Sokal-Sneath-3|\eqn{\frac{1}{4}\left(\frac{a}{R_1}+\frac{a}{C_1}+\frac{d}{R_2}+\frac{d}{C_2}\right)}|
#'  |32|Hawkin-Dotson|\eqn{\frac{1}{2}\left(\frac{a}{a+b+c}+\frac{d}{b+c+d}\right)}|
#'  |33|Clement|\eqn{\frac{aR_2}{nR_1}+\frac{dR_1}{nR_2}}|
#'  |34|Harris-Lahey|\eqn{\frac{a\left(R_2+C_2\right)}{2n\left(a+b+c\right)}+\frac{d\left(R_1+C_1\right)}{2n\left(b+c+d\right)}}|
#'  |35|Austin-Colwell|\eqn{\frac{2}{\pi}\text{arcsin}\sqrt{\frac{a+d}{n}}}|
#'  |36|Baroni-Urbani-Buser-1|\eqn{\frac{\sqrt{ad}+a}{\sqrt{ad}+a+b+c}}|
#'  |37|Peirce-1|\eqn{\frac{ad-bc}{R_1 R_2}}|
#'  |38|Peirce-2|\eqn{\frac{ad-bc}{C_1 C_2}}|
#'  |39|Cole C1|\eqn{\frac{ad-bc}{R_1 C_1}}|
#'  |40|Loevinger|\eqn{\frac{ad-bc}{\min\left(R_1 C_2, R_2 C_1\right)}}|
#'  |41|Cole C7|\eqn{\begin{cases} \frac{ad-bc}{R_1 C_2} & \text{ if } ad \geq bc \\ \frac{ad-bc}{R_1 C_1} & \text{ if } ad < bc \text{ and } a \leq d \\ \frac{ad-bc}{R_2 C_2} & \text{ if } ad < bc \text{ and } a > d\end{cases}}|
#'  |42|Dennis|\eqn{\frac{ad-bc}{\sqrt{nR_1 C_1}}}|
#'  |43|Phi|\eqn{\frac{ad-bc}{\sqrt{R_1 R_2 C_1 C_2}}}|
#'  |44|Doolittle|\eqn{\frac{\left(ad-bc\right)^2}{R_1 R_2 C_1 C_2}}|
#'  |45|Peirce-3|\eqn{\frac{ad-bc}{ab+2bc+cd}}|
#'  |46|Cohen-kappa|\eqn{\frac{2\left(ad-bc\right)}{R_1 C_2 + R_2 C_1}}|
#'  |47|McEwen-Michael|\eqn{\frac{4\left(ad-bc\right)}{\left(a+d\right)^2 + \left(b+c\right)^2}}|
#'  |48|Kuder-Richardson|\eqn{\frac{4\left(ad-bc\right)}{R_1 R_2 + C_1 C_2 + 2ad-2bc}}|
#'  |49|Scott|\eqn{\frac{4ad-\left(b+c\right)^2}{\left(R_1 + C_1\right)\left(R_2 + C_2\right)}}|
#'  |50|Maxwell-Pilliner|\eqn{\frac{2\left(ad-bc\right)}{R_1 R_2 + C_1 C_2}}|
#'  |51|Cole C5|\eqn{\frac{\sqrt{2}\left(ad-bc\right)}{\sqrt{\left(ad-bc\right)^2+R_1 R_2 C_1 C_2}}}|
#'  |52|Hamann|\eqn{\frac{\left(a+d\right)-\left(b+c\right)}{n}}|
#'  |53|Fleiss|\eqn{\frac{\left(ad-bc\right)\left(R_1 C_2 + R_2 C_1\right)}{2R_1 R_2 C_1 C_2}}|
#'  |54|Yule Q|\eqn{\frac{ad-bc}{ad+bc}}|
#'  |55|Yule Y|\eqn{\frac{\sqrt{ad}-\sqrt{bc}}{\sqrt{ad}+\sqrt{bc}}}|
#'  |56|Digby H|\eqn{\frac{\left(ad\right)^{3/4}-\left(bc\right)^{3/4}}{\left(ad\right)^{3/4}+\left(bc\right)^{3/4}}}|
#'  |57|Edward Q|\eqn{\frac{OR^{\pi/4}-1}{OR^{\pi/4}+1}}|
#'  |58|Tarwid|\eqn{\frac{na-R_1 C_1}{na+R_1 C_1}}|
#'  |59|Bonett-Price-1|\eqn{\frac{\hat{w}^x-1}{\hat{w}^x+1}}|
#'  |60|Contingency coefficient|\eqn{\sqrt{\frac{\chi^2}{n+\chi^2}}}|
#'  |61|Cohen w|\eqn{\sqrt{\frac{\chi^2}{\chi^2}}}|
#'  |62|Pearson|\eqn{\sqrt{\frac{\phi^2}{n+\phi^2}}}|
#'  |63|Hurlbert|\eqn{\frac{ad-bc}{\left\lvert ad-bc\right\rvert}\sqrt{\frac{\chi^2-\chi_{min}^2}{\chi_{max}^2-\chi_{min}^2}}}|
#'  |64|Stiles|\eqn{\log_{10}\left(\frac{n\left(\left\lvert ad-bc \right\rvert -\frac{n}{2}\right)^2}{R_1 R_2 C_1 C_2}\right)}|
#'  |65|McConnaughey|\eqn{\frac{a^2-bc}{R_1 C_1}}|
#'  |66|Baroni-Urbani-Buser-2|\eqn{\frac{a-b-c+\sqrt{ad}}{a+b+c+\sqrt{ad}}}|
#'  |67|Kent-Foster-1|\eqn{\frac{-bc}{bR_1 + cC_1 + bc}}|
#'  |68|Kent-Foster-2|\eqn{\frac{-bc}{bR_2 + cC_2 + bc}}|
#'  |69|Tulloss|\eqn{\sqrt{U\times S\times R}}|
#'  |70|Gilbert-Wells|\eqn{\ln\left(a\right)-\ln\left(n\right)-\ln\left(\frac{R_1}{n}\right)-\ln\left(\frac{C_1}{n}\right)}|
#'  |71|Yule r|\eqn{\cos\left(\frac{\pi\sqrt{bc}}{\sqrt{ad}+\sqrt{bc}}\right)}|
#'  |72|Anderberg|\eqn{\frac{\sigma-\sigma'}{2n}}|
#'  |73|Alroy F|\eqn{\frac{a\left(n' + \sqrt{n'}\right)}{a\left(n' + \sqrt{n'}\right)+\frac{3}{2}bc}}|
#'  |74|Pearson Q1|\eqn{\sin\left(\frac{\pi}{2}\times\frac{R_2 C_1}{R_1 C_2}\right)}|
#'  |75|Goodman-Kruskal Lambda-1|\eqn{\frac{\sigma-\sigma'}{2n-\sigma'}}|
#'  |76|Goodman-Kruskal Lambda-2|\eqn{\frac{2\min\left(a,d\right)-b-c}{2\min\left(a,d\right)+b+c}}|
#'  |77|Odds Ratio|\eqn{\frac{ad}{bc}}|
#'  |78|Pearson Q4|\eqn{\sin\frac{\pi}{2}\times\frac{1}{1+\frac{2bcn}{\left(ad-bc\right)\left(b+c\right)}}}|
#'  |79|Pearson Q5|\eqn{\sin\frac{\pi}{2}\times\frac{1}{\sqrt{1+\frac{4abcdn^2}{\left(ad-bc\right)^2\left(a+d\right)\left(b+c\right)}}}}|
#'  |80|Camp (3 ver.)|\eqn{\frac{m}{\sqrt{1+\Theta\times m^2}}}|
#'  |81|Becker-Clogg-1|\eqn{\frac{g-1}{g+1}}|
#'  |82|Becker-Clogg-2|\eqn{\frac{OR^{13.3/ \delta}-1}{OR^{13.3/ \delta}+1}}|
#'  |83|Bonett-Price-r|\eqn{\cos\left(\frac{\pi}{1+\omega^c}\right)}|
#'  |84|Bonett-Price-rhat|\eqn{\cos\left(\frac{\pi}{1+\hat{\omega}^{\hat{c}}}\right)}|
#'  |85|Chen-Popovich|\eqn{\frac{ad-bc}{\lambda_x \lambda_y n^2}}|
#' 
#' 
#'  **Equation 57**
#'  
#'  \deqn{OR = \frac{ad}{bc}}
#'  
#'  **Equation 59**
#'  
#'  \deqn{x = \frac{1}{2}-\left(\frac{1}{2}-p_{min}\right)^2}
#' \deqn{p_{min} = \frac{\min\left(R_1, R_2, C_1, C_2\right)}{n}}
#' \deqn{\hat{\omega} = \frac{\left(a+0.1\right)\times\left(d+0.1\right)}{\left(b+0.1\right)\times\left(c+0.1\right)}}
#' 
#' **Equations 60, 61, and 63**
#' 
#' \deqn{\chi^2 = \frac{n\left(ad-bc\right)^2}{R_1 R_2 C_1 C_2}}
#' 
#' **Equation 62**
#' 
#' \deqn{\Phi = \frac{\left\lvert ad-bc\right\rvert}{\sqrt{R_1 R_2 C_1 C_2}}}
#' Note that Choi et. al ommit the absolute value, but this would create problems with taking the square root if bc>ad.
#' 
#' 
#' **Equation 63:**
#' 
#' \deqn{\chi_{max}^2 = \begin{cases} \frac{nR_1 C_2}{R_2 C_1} & \text{ if } ad \geq bc \\  \frac{nR_1 C_1}{R_2 C_2} & \text{ if } ad < bc \text{ and } a \leq d \\  \frac{nR_2 C_2}{R_1 C_1} & \text{ if } ad < bc \text{ and } a > d \end{cases}}
#' \deqn{\chi_{min}^2 = \frac{n^3\left(\hat{a} - g\left(\hat{a}\right)\right)^2}{R_1 R_2 C_1 C_2}}
#' \deqn{\hat{a}=\frac{R_1 C_1}{n}}
#' \deqn{g\left(\hat{a}\right) = \begin{cases} \lfloor \hat{a} \rfloor & \text{ if } ad < bc \\  \lceil \hat{a}\rceil & \text{ if } ad \geq bc \end{cases}}
#' 
#' **Equation 69**
#' 
#' \deqn{U = \log_2\left(1+\frac{\min\left(b,c\right)+a}{\max\left(b,c\right)+a}\right)}
#' \deqn{S = \frac{1}{\sqrt{\log_2\left(2+\frac{\min\left(b,c\right)}{a+1}\right)}}}
#' \deqn{R = \log_2\left(1+\frac{a}{R_1}\right)\log_2\left(1+\frac{a}{RC1}\right)}
#' 
#' **Equation 71 and 75**
#' 
#' \deqn{\sigma = \max\left(a,b\right)+\max\left(c,d\right)+\max\left(a,c\right)+\max\left(b,d\right)}
#' \deqn{\sigma' = \max\left(R_1, R_2\right)+\max\left(C_1, C_2\right)}
#' 
#' **Equation 73**
#' 
#' \deqn{n' = a+b+c}
#' 
#' **Equation 80**
#' 
#' Camp (1934, pp. 309) describes the following steps for the calculation:
#' Step 1: If total of column 1 (C1) is less than column 2 (C2), swop the two columns
#' 
#' Step 2: Calculate \eqn{p = \frac{C1}{n}}, \eqn{p_1 = \frac{a}{n}}, and \eqn{p_2 = \frac{c}{C2}}
#' 
#' Step 3: Determine \eqn{z_1}, \eqn{z_2} as the normal deviate 
#' corresponding to the area \eqn{p_1}, \eqn{p_2} resp. (inverse standard normal cumulative distribution)
#' 
#' Step 4: Determine y the normal ordinate corresponding to \eqn{p} (the height of the normal distribution)
#' 
#' Step 5: Calculate \eqn{m = \frac{p\times\left(1-p\right)\times\left(z_1 + z_2\right)}{y}}
#' 
#' Step 6: Find phi in a table of phi values
#' 
#' Camp suggested for a very basic approximation to simply use \eqn{\phi=1}.
#' 
#' For a better approximation Camp made the following table:
#' 
#' | p | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 |
#' |---|-----|-----|-----|-----|-----|
#' |phi|0.637|0.63|0.62|0.60|0.56|
#' 
#' Cureton (1968, p. 241) expanded on this table and produced:
#' 
#' | p | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
#' |-----|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
#' | 0.5 | 0.637 | 0.636 | 0.636 | 0.635 | 0.635 | 0.634 | 0.634 | 0.633 | 0.633 | 0.632 | 0.631 |
#' | 0.6 | 0.631 | 0.631 | 0.630 | 0.629 | 0.628 | 0.627 | 0.626 | 0.625 | 0.624 | 0.622 | 0.621 |
#' | 0.7 | 0.621 | 0.620 | 0.618 | 0.616 | 0.614 | 0.612 | 0.610 | 0.608 | 0.606 | 0.603 | 0.600 |
#' | 0.8 | 0.600 | 0.597 | 0.594 | 0.591 | 0.587 | 0.583 | 0.579 | 0.574 | 0.569 | 0.564 | 0.559 |
#' 
#' Step 7: Calculate \eqn{r_t = \frac{m}{\sqrt{1+\phi\times m^2}}}
#' 
#' 
#' Cureton (1968) describes quite a few shortcomings with this approximation, and circumstances when it might be appropriate.
#' 
#' **Equation 81 and 82**
#' 
#' Version 81 will calculate:
#' \deqn{\rho^* = \frac{g-1}{g+1}}
#' 
#' Version 82 will calculate:
#' \deqn{\rho^{**} = \frac{OR^{13.3/\Delta} - 1}{OR^{13.3/\Delta} + 1}}
#' 
#' With:
#' \deqn{g=e^{12.4\times\phi - 24.6\times\phi^3}}
#' \deqn{\phi = \frac{\ln\left(OR\right)}{\Delta}}
#' \deqn{OR=\frac{\left(\frac{a}{c}\right)}{\left(\frac{b}{d}\right)} = \frac{a\times d}{b\times c}}
#' \deqn{\Delta = \left(\mu_{R1} - \mu_{R2}\right) \times \left(v_{C1} - v_{C2}\right)}
#' \deqn{\mu_{R1} = \frac{-e^{-\frac{t_r^2}{2}}}{p_{R1}}, \mu_{R2} = \frac{e^{-\frac{t_r^2}{2}}}{p_{R2}}}
#' \deqn{v_{C1} = \frac{-e^{-\frac{t_c^2}{2}}}{p_{C1}}, v_{C2} = \frac{e^{-\frac{t_c^2}{2}}}{p_{C2}}}
#' \deqn{t_r = \Phi^{-1}\left(p_{R1}\right), t_c = \Phi^{-1}\left(p_{C1}\right)}
#' \deqn{p_{x} = \frac{x}{n}}
#' 
#' **Equations 83 and 84**
#' 
#' Formula for version 1 is (Bonett & Price, 2005, p. 216):
#' \deqn{\rho^* = \cos\left(\frac{\pi}{1+\omega^c}\right)}
#' With:
#' \deqn{\omega = OR = \frac{a\times d}{b\times c}}
#' \deqn{c = \frac{1-\frac{\left|R_1-C_1\right|}{5\times n} - \left(\frac{1}{2}-p_{min}\right)^2}{2}}
#' \deqn{p_{min} = \frac{\text{MIN}\left(R_1, R_2, C_1, C_2\right)}{n}}
#' 
#' Formula for version 2 is  (Bonett & Price, 2005, p. 216):
#' \deqn{\hat{\rho}^* = \cos\left(\frac{\pi}{1+\hat{\omega}^{\hat{c}}}\right)}
#' with:
#' \deqn{\hat{\omega} = \frac{\left(a+\frac{1}{2}\right)\times \left(d+\frac{1}{2}\right)}{\left(b+\frac{1}{2}\right)\times \left(c+\frac{1}{2}\right)}}
#' \deqn{\hat{c} = \frac{1-\frac{\left|R_1-C_1\right|}{5\times \left(n+2\right)} - \left(\frac{1}{2}-\hat{p}_{min}\right)^2}{2}}
#' \deqn{\hat{p}_{min} = \frac{\text{MIN}\left(R_1, R_2, C_1, C_2\right)+1}{n+2}}
#' 
#' **Equation 85**
#' 
#' \deqn{\lambda_x = \Phi^{-1}\left(\frac{R_1}{n}\right)}
#' \deqn{\lambda_y = \Phi^{-1}\left(\frac{C_1}{n}\right)}
#' 
#' with \eqn{\Phi^{-1}\left(\dots\right)} being the inverse standard normal cumulative distribution
#' 
#' **Sources for formulas**
#' 
#' The formulas were obtained from the following sources.
#' The columns W-C-H show which equation corresponds to my label in:
#' 
#' * W: Warrens (2008, pp. 219–222). Equation 4 from Warrens is the chi-square value and not added.
#' * C: Choi et al. (2010, pp. 44–45). Equations not added from this source are: Eq. 4 is a ‘3w Jaccard’, could not find a source for this and not added. Equation 12 is just the intersection (a), eq. 13 the innerproduct (a+d), and 66 the dispersion. Equation 51 is the chi-square value and measures 15 to 30 and 62 are just distance measures.
#' * H: Hubálek (Hubálek, 1982, pp. 671–673)
#' 
#' If no page number is listed in the original source, the formula was taken from Warrens, Choi et al. and/or Hubálek.
#' 
#' |nr|Label|Original|W|C|H|
#' |--|-----|-------|--|--|--|
#' |1|Russell-Rao|(Russell & Rao, 1940)|15|14|14|
#' |2|Dice-1|(Dice, 1945, p. 302)|17a|||
#' |3|Dice-2|(Dice, 1945, p. 302)|17b|||
#' |4|Braun-Blanquet|(Braun-Blanquet, 1932)|12|46|1|
#' |5|Simpson Similarity|(Simpson, 1943, p. 20, 1960, p. 301)|16|45|2|
#' |6|Kulczynski-1|(Kulczynski, 1927)|11b|64|3|
#' |7|Jaccard =|(Jaccard, 1901, 1912, p. 39)|6|1|4|
#' ||Tanimoto|(Tanimoto, 1958, p. 5)||65||
#' |8|Sokal-Sneath-1|(Sokal & Sneath, 1963, p. 129)|30a|6|6|
#' |9|Gleasson =|(Gleason, 1920, p. 31)|9||5|
#' ||Dice-3 =|(Dice, 1945, p. 302)||2||
#' ||Nei-Li =|(Nei & Li, 1979, p. 5270)||5||
#' ||Czekanowski|||3||
#' |10|Mountford|(Mountford, 1962, p. 45)|28|37|15|
#' |11|Driver-Kroeber =|(Driver & Kroeber, 1932, p. 219)|13|31|11|
#' ||Ochiai-1  =|(Ochiai, 1957)||33||
#' ||Otsuka|(Otsuka, 1936)||38||
#' |12|Sorgenfrei|(Sorgenfrei, 1958)|23|36|12|
#' |13|Johnson|(Johnson, 1967)|33|43|9|
#' |14|Driver-Kroeber-2 = |(Driver & Kroeber, 1932, p. 219)|11a|42|8|
#' ||Kulczynski-2|(Kulczynski, 1927)||41|7|
#' |15|Fager-McGowan-1|(Fager & McGowan, 1963, p. 454)|29|||
#' |16|Fager-McGowan-2|(Fager & McGowan, 1963, p. 454)||47|13|
#' |17|tarantula|(Jones & Harrold, 2005)||75||
#' |18|ample|||76||
#' |19|Gilbert|(Gilbert, 1884, p. 171)||||
#' |20|Fossum-Kaskey|(Fossum & Kaskey, 1966, p. 65)||35||
#' |21|Forbes - 1|(Forbes, 1907, p. 279)|5|34|40|
#' |22|Eyraud|(Eyraud, 1936)||74|17|
#' |23|Sokal-Michener|(Sokal & Michener, 1958, p. 1417)|22|7|20|
#' |24|Faith|(Faith, 1983, p. 290)||10||
#' |25|Sokal-Sneath-5|(Sokal & Sneath, 1963, p. 129)|30e|56|19|
#' |26|Rogers-Tanimoto|(Rogers & Tanimoto, 1960)|25|9|23|
#' |27|Sokal-Sneath-2 = |(Sokal & Sneath, 1963, p. 129)|30b|8|22|
#' ||Gower-Legendre|(Gower & Legendre, 1986)||11||
#' |28|Gower|(Gower, 1971)||50||
#' |29|Sokal-Sneath-4 = |(Sokal & Sneath, 1963, p. 130)|30d|57|25|
#' ||Ochiai-2|(Ochiai, 1957)||60||
#' |30|Rogot-Goldberg|(Rogot & Goldberg, 1966, p. 997)|32|||
#' |31|Sokal-Sneath-3|(Sokal & Sneath, 1963, p. 130)|30c|49|18|
#' |32|Hawkin-Dotson|(Hawkins & Dotson, 1975, pp. 372–373)|34|||
#' |33|Clement|(Clement, 1976, p. 258)|37|||
#' |34|Harris-Lahey|(Harris & Lahey, 1978, p. 526)|40|||
#' |35|Austin-Colwell|(Austin & Colwell, 1977, p. 205)|||21|
#' |36|Baroni-Urbani-Buser-1|(Baroni-Urbani & Buser, 1976, p. 258)|38a|71|32|
#' |37|Peirce-1|(Peirce, 1884, p. 453)|1a|||
#' |38|Peirce-2|(Peirce, 1884, p. 453)|1b||26|
#' |39|Cole C1|(Cole, 1949, p. 415)||||
#' |40|Loevinger = |(Loevinger, 1947, p. 30)|18|||
#' ||Forbes 2|(Forbes, 1925)||48|42|
#' |41|Cole C7|(Cole, 1949, p. 420)|19||34|
#' |42|Dennis|(Dennis, 1965, p. 69)||44||
#' |43|Pearson Phi =|(Pearson, 1900a, p. 12)|7|54|30|
#' ||Yule Phi =|(Yule, 1912, p. 596)||||
#' ||Cole C2|(Cole, 1949, p. 415)||||
#' |44|Doolittle|(Doolittle, 1885, p. 123)|2||31|
#' |45|Peirce-3|(Peirce, 1884)||73|16|
#' |46|Cohen-kappa|(Cohen, 1960, p. 40)|24|||
#' |47|McEwen-Michael =|(Michael, 1920, p. 57)|10|68|39|
#' ||Cole C3|(Cole, 1949, p. 415)||||
#' |48|Kuder-Richardson|(Kuder & Richardson, 1937)|14|||
#' |49|Scott|(Scott, 1955, p. 324)|21|||
#' |50|Maxwell-Pilliner|(Maxwell & Pilliner, 1968)|35|||
#' |51|Cole C5|(Cole, 1949, p. 416)||58|29|
#' |52|Hamann|(Hamann, 1961)|27|67|24|
#' |53|Fleiss|(Fleiss, 1975, p. 656)|36|||
#' |54|Yule Q =|(Yule, 1900, p. 272)|3|61|36|
#' ||Cole C4 =|(Cole, 1949, p. 415)||||
#' ||Pearson Q2|(Pearson, 1900, p. 15)||||
#' |55|Yule Y|(Yule, 1912, p. 592)|8|63|37|
#' |56|Digby H|(Digby, 1983, p. 754)|41|||
#' |57|Edward Q|(Edwards, 1957; Becker & Clogg, 1988, p. 409)||||
#' |58|Tarwid|(Tarwid, 1960, p. 117)||40|43|
#' |59|Bonett-Price-1|(Bonett & Price, 2007, p. 433)||||
#' |60|Contingency|(Pearson, 1904, p. 9)||52|28|
#' |61|Cohen w|(Cohen, 1988, p. 216)||||
#' |62|Pearson|(Pearson, 1904)||53||
#' |63|Hurlbert / Cole C8|(Hurlbert, 1969, p. 1)|||35|
#' |64|Stiles|(Stiles, 1961, p. 272)|26|59||
#' |65|McConnaughey|(McConnaughey, 1964)|31|39|10|
#' |66|Baroni-Urbani-Buser-2|(Baroni-Urbani & Buser, 1976, p. 258)|38b|72|33|
#' |67|Kent-Foster-1|(Kent & Foster, 1977, p. 311)|39a|||
#' |68|Kent-Foster-2|(Kent & Foster, 1977, p. 311)|39b|||
#' |69|Tulloss|(Tulloss, 1997, p. 133)||||
#' |70|Gilbert-Wells|(Gilbert & Wells, 1966)||||
#' |71|Yule r|(Yule, 1900, p. 276)||||
#' | |Pearson-Q3|(Pearson, 1900a, p. 16)||||
#' | |Cole C6|(Cole, 1949, p. 416)||||
#' | |Pearson-Heron|(Pearson & Heron, 1913)||55|38|
#' |72|Anderberg|(Anderberg, 1973)||70||
#' |73|Alroy F|(Alroy, 2015, eq. 6)||||
#' |74|Pearson Q1|(Pearson, 1900a, p. 15)||||
#' |75|Goodman-Kruskal Lambda-1|(Goodman & Kruskal, 1954, p. 743)||69||
#' |76|Goodman-Kruskal Lambda-2|(Goodman & Kruskal, 1954)|20|||
#' |77|Odds Ratio|(Fisher, 1935, p. 50)||||
#' |78|Pearson Q4|(Pearson, 1900, p. 16)|   |   |   |
#' |79|Pearson Q5|(Pearson, 1900, p. 16)|   |   |   |
#' |80|Camp|(Camp, 1934, p. 309)|   |   |   |
#' |81|Becker-Clogg-1|(Becker & Clogg, 1988, pp. 410–412)|   |   |   |
#' |82|Becker-Clogg-2|(Becker & Clogg, 1988, pp. 410–412)|   |   |   |
#' |83|Bonett-Price-2|(Bonett & Price, 2005, p. 216)|   |   |   |
#' |84|Bonett-Price-3|(Bonett & Price, 2005, p. 216)|   |   |   |
#' |85|Ched-Popovich|(Chen & Popovich, 2002, p. 37)|   |   |   |
#' 
#' @references
#' Alroy, J. (2015). A new twist on a very old binary similarity coefficient. *Ecology, 96*(2), 575–586. doi:10.1890/14-0471.1
#' 
#' Anderberg, M. R. (1973). *Cluster analysis for applications*. New York, NY, Academic Press.
#' 
#' Austin, B., & Colwell, R. R. (1977). Evaluation of some coefficients for use in numerical taxonomy of microorganisms. *International Journal of Systematic Bacteriology, 27*(3), 204–210. doi:10.1099/00207713-27-3-204
#' 
#' Baroni-Urbani, C., & Buser, M. W. (1976). Similarity of binary data. *Systematic Zoology, 25*(3), 251–259. doi:10.2307/2412493
#' 
#' Becker, M. P., & Clogg, C. C. (1988). A note on approximating correlations from Odds Ratios. *Sociological Methods & Research, 16*(3), 407–424. doi:10.1177/0049124188016003003
#' 
#' Bonett, D. G., & Price, R. M. (2005). Inferential methods for the tetrachoric correlation coefficient. *Journal of Educational and Behavioral Statistics, 30*(2), 213–225. https://doi.org/10.3102/10769986030002213
#' 
#' Bonett, D. G., & Price, R. M. (2007). Statistical inference for generalized yule coefficients in 2 x 2 contingency tables. *Sociological Methods & Research, 35*(3), 429–446. doi:10.1177/0049124106292358
#' 
#' Braun-Blanquet, J. (1932). *Plant sociology: The study of plant communities*. McGraw Hill.
#' 
#' Camp, B. H. (1934). *Mathematical part of elementary statistics*. D.C. Heath and Company, London.
#' 
#' Chen, P. Y., & Popovich, P. M. (2002). *Correlation: Parametric and nonparametric measures*. Sage Publications.
#' 
#' Choi, S.-S., Cha, S.-H., & Tappert, C. (2010). A survey of binary similarity and distance measures. *Journal on Systemics, Cybernetics and Informatics, 8*(1), 43–48.
#' 
#' Clement, P. W. (1976). A formula for computing inter-observer agreement. *Psychological Reports, 39*(1), 257–258. doi:10.2466/pr0.1976.39.1.257
#' 
#' Cohen, J. (1960). A coefficient of agreement for nominal scales. *Educational and Psychological Measurement, 20*(1), 37–46. doi:10.1177/001316446002000104
#' 
#' Cohen, J. (1988). *Statistical power analysis for the behavioral sciences* (2nd ed.). L. Erlbaum Associates.
#' 
#' Cole, L. C. (1949). *The measurement of interspecific associaton. Ecology, 30*(4), 411–424. doi:10.2307/1932444
#' 
#' Cureton, E. E. (1968). Tetrachoric correlation by the Camp approximation. *Educational and Psychological Measurement, 28*(2), 239–244. doi:10.1177/001316446802800202
#' 
#' Dennis, S. F. (1965). The construction of a thesaurus automatically from a sample of text. In M. E. Stevens, V. E. Giuliano, & L. B. Heilprin (Eds.), *Statistical Association Methods for Mechanized Documentation* (Vol. 14, pp. 61–148). U.S. Government Printing Office.
#' 
#' Dice, L. R. (1945). Measures of the amount of ecologic association between species. *Ecology, 26*(3), 297–302. doi:10.2307/1932409
#' 
#' Digby, P. G. N. (1983). Approximating the tetrachoric correlation coefficient. *Biometrics, 39*(3), 753–757. doi:10.2307/2531104
#' 
#' Doolittle, M. H. (1885). The verification of predictions. *Bulletin of the Philosophical Society of Washington, 7*, 122–127.
#' 
#' Driver, H. E., & Kroeber, A. L. (1932). Quantitative expression of cultural relationships. *University Of California Publications in American Archeology and Ethnology, 31*(4), 211–256.
#' 
#' Edwards, J. H. (1957). A note on the practical interpretation of 2 x 2 tables. *Journal of Epidemiology & Community Health, 11*(2), 73–78. doi:10.1136/jech.11.2.73
#' 
#' Edwards, J. H., & Edwards, A. W. F. (1984). Approximating the tetrachoric correlation coefficient. *Biometrics, 40*(2), 563–563.
#' 
#' Eyraud, H. (1936). Les principes de la mesure des correlations. *Ann. Univ. Lyon, III. Ser., Sect. A, 1*(30–47), 111.
#' 
#' Fager, E. W., & McGowan, J. A. (1963). Zooplankton species groups in the north pacific: Co-occurrences of species can be used to derive groups whose members react similarly to water-mass types. *Science, 140*(3566), 453–460. doi:10.1126/science.140.3566.453
#' 
#' Faith, D. P. (1983). Asymmetric binary similarity measures. *Oecologia, 57*(3), 287–290. doi:10.1007/BF00377169
#' 
#' Fisher, R. A. (1935). The logic of inductive inference. *Journal of the Royal Statistical Society, 98*(1), 39–82. doi:10.2307/2342435
#' 
#' Fleiss, J. L. (1975). Measuring agreement between two judges on the presence or absence of a trait. *Biometrics, 31*(3), 651–659. doi:10.2307/2529549
#' 
#' Forbes, S. A. (1907). On the local distribution of certain Illinois fishes: An essay in statistical ecology. *Illinois Natural History Survey Bulletin, 7*(8), 273–303.
#' 
#' Forbes, S. A. (1925). Method of determining and measuring the associative relations of species. *Science, 61*, 524.
#' 
#' Fossum, E. G., & Kaskey, G. (1966). Optimization and standardization of information retrieval language and systems (AF49(638)1194; p. 96). Univac Division.
#' 
#' Gilbert, G. K. (1884). Finley’s tornado predictions. *American Meteorological Journal, 1*(5), 166–172.
#' 
#' Gilbert, N., & Wells, T. C. E. (1966). *Analysis of quadrat data. Journal of Ecology, 54*(3), 675–685. doi:10.2307/2257810
#' 
#' Gleason, H. A. (1920). Some applications of the quadrat method. *Bulletin of the Torrey Botanical Club, 47*(1), 21–33. doi:10.2307/2480223
#' 
#' Glen, S. (2017, August 16). Gamma Coefficient (Goodman and Kruskal’s Gamma) & Yule’s Q. Statistics How To. https://www.statisticshowto.com/gamma-coefficient-goodman-kruskal/
#' 
#' Goodman, L. A., & Kruskal, W. H. (1954). Measures of association for cross classifications. *Journal of the American Statistical Association, 49*(268), 732–764. doi:10.2307/2281536
#' 
#' Gower, J. C. (1971). A general coefficient of similarity and some of its properties. *Biometrics, 27*(4), 857. doi:10.2307/2528823
#' 
#' Gower, J. C., & Legendre, P. (1986). Metric and Euclidean properties of dissimilarity coefficients. *Journal of Classification, 3*(1), 5–48. doi:10.1007/BF01896809
#' 
#' Hamann, U. (1961). Merkmalsbestand und verwandtschaftsbeziehungen der farinosae: Ein beitrag zum system der monokotyledonen. *Willdenowia, 2*(5), 639–768.
#' 
#' Harris, F. C., & Lahey, B. B. (1978). A method for combining occurrence and nonoccurrence interobserver agreement scores. *Journal of Applied Behavior Analysis, 11*(4), 523–527. doi:10.1901/jaba.1978.11-523
#' 
#' Hawkins, R. P., & Dotson, V. A. (1975). Reliability scores that delude: An Alice in wonderland trip through the misleading characteristics of inter-observer agreement scores in interval recording. In E. Ramp & G. Semb (Eds.), *Behavior  analysis: Areas of research and application *(pp. 359–376). Prentice Hall.
#' 
#' Hubálek, Z. (1982). Coefficients of association and similarity, based on binary (presence-absence) data: An evaluation. *Biological Reviews, 57*(4), 669–689. doi:10.1111/j.1469-185X.1982.tb00376.x
#' 
#' Hurlbert, S. H. (1969). A coefficient of interspecific assciation. *Ecology, 50*(1), 1–9. doi:10.2307/1934657
#' 
#' Jaccard, P. (1901). Étude comparative de la distribution florale dans une portion des Alpes et des Jura. *Bulletin Del La Société Vaudoise Des Sciences Naturelles, 37*, 547–579.
#' 
#' Jaccard, P. (1912). The distribution of the flora in the alpine zone. *The New Phytologist, 11*(2), 37–50.
#' 
#' Johnson, S. C. (1967). Hierarchical clustering schemes. *Psychometrika, 32*(3), 241–254. doi:10.1007/BF02289588
#' 
#' Jones, J. A., & Harrold, M. J. (2005). Empirical evaluation of the tarantula automatic fault-localization technique. Proceedings of the 20th IEEE/ACM International Conference on Automated Software Engineering, 273–282. doi:10.1145/1101908.1101949
#' 
#' Kent, R. N., & Foster, S. L. (1977). Direct observational procedures: Methodological issues in naturalistic settings. In A. R. Ciminero, K. S. Calhoun, & H. E. Adams (Eds.), *Handbook of behavioral assessment* (pp. 279–328). New York, NY: Wiley. http://archive.org/details/handbookofbehavi00cimi
#' 
#' Kuder, G. F., & Richardson, M. W. (1937). The theory of the estimation of test reliability. *Psychometrika, 2*(3), 151–160. doi:10.1007/BF02288391
#' 
#' Kulczynski, S. (1927). Die Pflanzenassoziationen der Pieninen. *Bulletin International de l’Academie Polonaise Des Sciences et Des Lettres, Classe Des Sciences Mathematiques et Naturelles, B (Sciences Naturelles)*, II, 57–203.
#' 
#' Ling, M. H. T. (2010). COPADS, I: Distance coefficients between two lists or sets. *The Python Papers Source Codes, 2*(2), 1–31.
#' 
#' Loevinger, J. (1947). A systematic approach to the construction and evaluation of tests of ability. *Psychological Monographs, 61*(4), i–49. doi:10.1037/h0093565
#' 
#' Maxwell, A. E., & Pilliner, A. E. (1968). Deriving coefficients of reliability and agreement for ratings. *The British Journal of Mathematical and Statistical Psychology, 21*(1), 105–116. doi:10.1111/j.2044-8317.1968.tb00401.x
#' 
#' McConnaughey, B. H. (1964). The determination and analysis of plankton communities. *Marine Research, 7*, 1–40.
#' 
#' Michael, E. L. (1920). Marine Ecology and the coefficient of association: A plea in behalf of quantitative biology. *Journal of Ecology, 8*(1), 54–59. doi:10.2307/2255213
#' 
#' Mountford, M. D. (1962). An index of similarity and its application to classification problems. In P. W. Murphy & D. Phil (Eds.), *Progress in soil zoology* (pp. 43–50). Butterworths.
#' 
#' Nei, M., & Li, W. H. (1979). Mathematical model for studying genetic variation in terms of restriction endonucleases. *Proceedings of the National Academy of Sciences, 76*(10), 5269–5273. doi:10.1073/pnas.76.10.5269
#' 
#' Ochiai, A. (1957). Zoogeographical studies on the soleoid fishes found in Japan and its neighbouring regions-I. *Nippon Suisan Gakkaishi, 22*(9), 522–525. doi:10.2331/suisan.22.522
#' 
#' Otsuka, Y. (1936). The faunal character of the Japanese Pleistocene marine Mollusca, as evidence of the climate having become colder during the Pleistocene in Japan. *Bulletin of the Biogeographical Society of Japan, 6*(16), 165–170.
#' 
#' Pearson, K. (1900). Contributions to the mathematical theory of evolution. VII. On the correlation of characters not quantitatively measurable. *Philosophical Transactions of the Royal Society of London, 195*, 1–405. doi:10.1098/rsta.1900.0022
#' 
#' Pearson, K. (1904). *Contributions to the Mathematical Theory of Evolution. XIII. On the theory of contingency and its relation to association and normal correlation*. Dulau and Co.
#' 
#' Pearson, K., & Heron, D. (1913). On theories of association. *Biometrika, 9*(1/2), 159–315. doi:10.2307/2331805
#' 
#' Peirce, C. S. (1884). The numerical measure of the success of predictions. *Science, 4*(93), 453–454. doi:10.1126/science.ns-4.93.453-a
#' 
#' Rogers, D. J., & Tanimoto, T. T. (1960). A computer program for classifying plants. *Science, 132*(3434), 1115–1118. doi:10.1126/science.132.3434.1115
#' 
#' Rogot, E., & Goldberg, I. D. (1966). A proposed index for measuring agreement in test-retest studies. *Journal of Chronic Diseases, 19*(9), 991–1006. doi:10.1016/0021-9681(66)90032-4
#' 
#' Russell, P. F., & Rao, T. R. (1940). On habitat and association of species of anopheline larvae in south-eastern Madras. *Journal of the Malaria Institute of India, 3*(1), 153–178.
#' 
#' Scott, W. A. (1955). Reliability of content analysis: The case of nominal scale coding. *The Public Opinion Quarterly, 19*(3), 321–325.
#' 
#' Simpson, G. G. (1943). Mammals and the nature of continents. *American Journal of Science, 241*(1), 1–31. doi:10.2475/ajs.241.1.1
#' 
#' Simpson, G. G. (1960). Notes on the measurement of faunal resemblance. *American Journal of Science, 258-A*, 300–311.
#' 
#' Sokal, P. H. A., & Sneath, R. R. (1963). Principles of numerical taxonomy. W.H. Freeman and Company.
#' 
#' Sokal, R., & Michener, C. (1958). A statistical method for evaluating systematic relationships. University of Kansas Science Bulletin.
#' 
#' Sorgenfrei, T. (1958). Molluscan assemblages from the marine middle miocene of south Jutland and their Environments. Vol. I. *Danmarks Geologiske Undersøgelse II. Række*, 79, 1–355. doi:10.34194/raekke2.v79.6868
#' 
#' Stiles, H. E. (1961). The association factor in information retrieval. *Journal of the ACM, 8*(2), 271–279. doi:10.1145/321062.321074
#' 
#' Tanimoto, T. T. (1958). An elementary mathematical theory of classification and prediction (PB167360). International Business Machines Corp., New York, NY.
#' 
#' Tarwid, K. (1960). Szacowanie zbieinosci nisz ekologicznych gatunkow droga oceny prawdopodobienstwa spotykania sie ich w polowach. *Ekologia Polska, 6*, 115–130.
#' 
#' Tulloss, R. E. (1997). Assessment of similarity indices for undesirable properties and a new tripartite similarity index based on cost functions. In M. E. Palm & I. H. Chapela (Eds.), *Mycology in sustainable development: Expanding concepts, vanishing borders* (pp. 122–143). Parkway Pub.
#' 
#' Walker, H. M., & Lev, J. (1953). *Statistical inference*. Holt. 
#' 
#' Warrens, M. J. (2008). Similarity coefficients for binary data: Properties of coefficients, coefficient matrices, multi-way metrics and multivariate coefficients [Doctoral dissertation, Universiteit Leiden]. https://hdl.handle.net/1887/12987
#' 
#' Yule, G. U. (1900). On the association of attributes in statistics: With illustrations from the material of the childhood society, &c. *Philosophical Transactions of the Royal Society of London, 194*, 257–319. doi:10.1098/rsta.1900.0019
#' 
#' Yule, G. U. (1912). On the methods of measuring association between two attributes. *Journal of the Royal Statistical Society, 75*(6), 579–652. doi:10.2307/2340126
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' 
#' @export
es_bin_bin <- function(field1, 
                       field2, 
                       categories1=NULL, 
                       categories2=NULL, 
                       method = "odds-ratio"){

  #get the cross table
      tb = tab_cross(field1, field2, order1=categories1, order2=categories2)
      
  #cell values of sample cross table
  a = tb[1, 1]
  b = tb[1, 2]
  c = tb[2, 1]
  d = tb[2, 2]
  r1 = a + b
  r2 = c + d
  c1 = a + c
  c2 = b + d
  n = r1 + r2
  
  if (method == "russell-rao"){
    es = a / n
  } else if (method == "dice-1"){
    es = a / r1
  } else if (method == "dice-2"){
    es = a / c1
  } else if (method == "braun-blanquet"){
      if (r1 > c1){es = a / r1} 
      else {es = a / c1}
  } else if (method == "simpson"){
      if (r1 < c1){es = a / r1}
      else {es = a / c1}
  } else if (method == "kulczynski-1"){
      es = a / (b + c)
  } else if (method == "jaccard" || method == "tanimoto"){
      es = a / (a + b + c)
  } else if (method == "sokal-sneath-1"){
      es = a / (a + 2 * b + 2 * c)
  } else if (method == "gleason" || method == "dice-3" || method == "nei-li" || method == "czekanowski"){
      es = 2 * a / (2 * a + b + c)
  } else if (method == "mountford"){
      es = 2 * a / (a * (b + c) + 2 * b * c)
  } else if (method == "driver-kroeber-1" || method == "ochiai-1" || method == "otsuka"){
      es = a / sqrt((a + b) * (a + c))
  } else if (method == "sorgenfrei"){
      es = a ^ 2 / (r1 * c1)
  } else if (method == "johnson"){
      es = a / r1 + a / c1
  } else if (method == "kulczynski-2"||method == "driver-kroeber-2"){
      es = a * (2 * a + b + c) / (2 * r1 * c1)
  } else if (method == "fager-mcgowan-1"){
      es = a / sqrt(r1 * c1) - 1 / (2 * sqrt(max(r1, c1)))
  } else if (method == "fager-mcgowan-2"){
      es = a / sqrt(r1 * c1) - sqrt(max(r1, c1)) / 2
  } else if (method == "tarantula"){
      es = a * r2 / (c * r1)
  } else if (method == "ample"){
      es = abs(a * r2 / (c * r1))
  } else if (method == "gilbert"){
      es = (a * n - r1 * c1) / (c1 * n + r1 * n - a * n - r1 * c1)
  } else if (method == "fossum-kaskey"){
      es = n * (a - 1 / 2) ^ 2 / (r1 * c1)
  } else if (method == "eyraud"){
      es = (a - r1 * c1) / (r1 * r2 * c1 * c2)
  } else if (method == "sokal-michener"||method == "matching"){
      es = (a + d) / (a + b + c + d)
  } else if (method == "faith"){
      es = (a + 1 / 2 * d) / n
  } else if (method == "sokal-sneath-5"){
      es = (a + d) / (b + c)
  } else if (method == "rogers-tanimoto"){
      es = (a + d) / (a + 2 * b + 2 * c + d)
  } else if (method == "sokal-sneath-2" || method == "gower-legendre"){
      es = (2 * a + 2 * d) / (2 * a + b + c + 2 * d)
  } else if (method == "gower"){
        es = (a + d) / sqrt(r1 * r2 * c1 * c2)
  } else if (method == "sokal-sneath-4" || method == "ochiai-2"){
      es = a * d / sqrt((a + b) * (a + c) * (c + d) * (b + d))
  } else if (method == "rogot-goldberg"){
      es = a / (r1 + c1) + d / (r2 + c2)
  } else if (method == "sokal-sneath-3"){
      es = 1 / 4 * (a / r1 + a / c1 + d / r2 + d / c2)
  } else if (method == "hawkins-dotson"){
      es = 1 / 2 * (a / (a + b + c) + d / (b + c + d))
  } else if (method == "clement"){
      es = a * r2 / (n * r1) + d * r1 / (n * r2)
  } else if (method == "harris-lahey"){
      es = a * (r2 + c2) / (2 * n * (a + b + c)) + d * (r1 + c1) / (2 * n * (b + c + d))
  } else if (method == "austin-colwell"){
      es = 2 / pi * asin(sqrt((a + d) / n))
  } else if (method == "forbes-1"){
      es = n * a / (r1 * c1)
  } else if (method == "baroni-urbani-buser-1"){
      es = (sqrt(a * d) + a) / (sqrt(a * d) + a + b + c)
  } else if (method == "peirce-1"){
      es = (a * d - b * c) / (r1 * r2)
  } else if (method == "peirce-2"){
      es = (a * d - b * c) / (c1 * c2)
  } else if (method == "cole-c1"){
      c1 = (a * d - b * c) / ((a + b) * (a + c))
      es = c1
  } else if (method == "loevinger" || method == "forbes-2"){
      es = (a * d - b * c) / min(r1 * c2, r2 * c1)
  } else if (method == "cole-c7"){
      if (a * d >= (b * c)){
        C7 = (a * d - b * c) / ((a + b) * (b + d))
      } else if (a <= d){
          C7 = (a * d - b * c) / ((a + b) * (a + c))
      } else {
          C7 = (a * d - b * c) / ((b + d) * (c + d))}
      es = C7
  } else if (method == "dennis"){
      es = (a * d - b * c) / sqrt(n * r1 * c1)
  } else if (method == "phi" || method == "cole-c2"){
      es = (a * d - b * c) / sqrt(r1 * r2 * c1 * c2)
  } else if (method == "doolittle"){
      es = (a * d - b * c) ^ 2 / (r1 * r2 * c1 * c2)
  } else if (method == "peirce-3"){
      es = (a * d + b * c) / (a * b + 2 * b * c + c * d)
  } else if (method == "cohen-kappa"){
      es = 2 * (a * d - b * c) / (r1 * c2 + r2 * c1)
  } else if (method == "mcewen-michael" || method == "cole-c3"){
      m = 4 * (a * d - b * c) / ((a + d) ^ 2 + (b + c) ^ 2)
      es = m
  } else if (method == "kuder-richardson"){
      es = 4 * (a * d - b * c) / (r1 * r2 + c1 * c2 + 2 * a * d - 2 * b * c)
  } else if (method == "scott"){
      es = (4 * a * d - (b + c) ^ 2) / ((r1 + c1) * (r2 + c2))
  } else if (method == "maxwell-pilliner"){
      es = 2 * (a * d - b * c) / (r1 * r2 + c1 * c2)
  } else if (method == "cole-c5"){
      C5 = sqrt(2) * (a * d - b * c) / sqrt((a * d - b * c) ^ 2 + (a + b) * (c + d) * (a + c) * (b + d))
      es = C5
  } else if (method == "hamann"){
      es = (a + d - b - c) / (a + b + c + d)
  } else if (method == "fleiss"){
      es = (a * d - b * c) * (r1 * c2 + r2 * c1) / (2 * r1 * r2 * c1 * c2)
  } else if (method == "yule-q" || method == "cole-c4" || method == "pearson-q2"){
      es = (a * d - b * c) / (a * d + b * c)
  } else if (method == "yule-y"){
      y = (sqrt(a * d) - sqrt(b * c)) / (sqrt(a * d) + sqrt(b * c))
      es = y
  } else if (method == "digby"){
        es = ((a * d) ^ (3 / 4) - (b * c) ^ (3 / 4)) / ((a * d) ^ (3 / 4) + (b * c) ^ (3 / 4))
  } else if (method == "edward"){
      OdRat = a * d / (b * c)
      es = (OdRat ^ (pi / 4) - 1) / (OdRat ^ (pi / 4) + 1)
  } else if (method == "tarwid"){
      es = (n * a - r1 * c1) / (n * a + r1 * c1)
  } else if (method == "mcconnaughey"){
      es = (a ^ 2 - b * c) / ((a + b) * (a + c))
  } else if (method == "baroni-urbani-buser-2"){
      es = (sqrt(a * d) + a - b - c) / (sqrt(a * d) + a + b + c)
  } else if (method == "kent-foster-1"){
      es = -b * c / (b * r1 + c * c1 + b * c)
  } else if (method == "kent-foster-2"){
      es = -b * c / (b * r2 + c * c2 + b * c)
  } else if (method == "tulloss"){
      mn = b
      mx = c
      if (c < b){
          mn = c
          mx = b
      }
      f1 = (log(1 + (mn + a) / (mx + a)) / log(10)) / (log(2) / log(10))
      f2 = 1 / sqrt((log(2 + mn / (a + 1)) / log(10)) / (log(2) / log(10)))
      f3 = (log(1 + a / (a + b)) / log(10)) * (log(1 + a / (a + c)) / log(10)) / ((log(2) / log(10)) ^ 2)
      
      es = sqrt(f1 * f2 * f3)
  } else if (method == "gilbert-wells"){
      es = log(a) - log(n) - log(r1 / n) - log(c1 / n)
  } else if (method == "yule-r" || method == "pearson-q3" || method == "cole-c6" || method == "pearson-heron"){
      es = cos(pi * sqrt(b * c) / (sqrt(a * d) + sqrt(b * c)))
  } else if (method == "anderberg"){
      S = max(a, b) + max(c, d) + max(a, c) + max(b, d)
      s2 = max(r1, r2) + max(c1, c2)
      es = (S - s2) / (2 * n)
  } else if (method == "alroy"){
      na = a + b + c
      es = a * (na + sqrt(na)) / (a * (na + sqrt(na)) + 3 / 2 * b * c)
  } else if (method == "pearson-q1"){
      #Pearson requires for Q1 that ad>bc, a>d, and c>b
      sw = 1
      Runs = 0
      while (Runs < 2){
        #swop rows
        if (a * d < (b * c) || a < d || c < b){
        #swop the rows
          at = a
          a = c
          c = at
          bt = b
          b = d
          d = bt
          sw = -sw
        }
        #swop columns
        if (a * d < (b * c) || a < d || c < b){
        #swop the columns
          at = a
          a = b
          b = at
          ct = c
          c = d
          d = ct
          sw = -sw
        }
        #swop rows again
        if (a * d < (b * c) || a < d || c < b){
        #swop the rows
          at = a
          a = c
          c = at
          bt = b
          b = d
          d = bt
          sw = -sw
        }
        #pivot the table
        if (a * d < (b * c) || a < d || c < b){
          bt = b
          b = c
          c = bt
          sw = -sw
        }
        Runs = Runs + 1
      }
      es = sw * sin(pi / 2 * (c + d) * (a + c) / ((a + b) * (b + d)))
        
  } else if (method == "gk-lambda-1"){
      S = max(a, b) + max(c, d) + max(a, c) + max(b, d)
      s2 = max(r1, r2) + max(c1, c2)
      es = (S - s2) / (2 * n - s2)
        
  } else if (method == "gk-lambda-2"){
      m = min(a, d)
      es = (2 * m - b - c) / (2 * m + b + c)
  } else if (method == "contingency"){
      chi = n * (a * d - b * c) ^ 2 / (r1 * r2 * c1 * c2)
      es = sqrt(chi / (n + chi))
  } else if (method == "cohen-w"){
      chi = n * (a * d - b * c) ^ 2 / (r1 * r2 * c1 * c2)
      es = sqrt(chi / n)
  } else if (method == "pearson"){
      Phi = abs(a * d - b * c) / sqrt(r1 * r2 * c1 * c2)
      es = sqrt(Phi / (n + Phi))
  } else if (method == "hurlbert" || method == "cole-c8"){
      chi = n * (a * d - b * c) ^ 2 / (r1 * r2 * c1 * c2)
      ahat = r1 * c1 / n
      gahat = floor(ahat)
      if (a * d >= b * c){
        chimax = n * r1 * c2 / r2 * c1
        gahat = ceiling(ahat)
      } else if (a * d < b * c && a <= d){
        chimax = n * r1 * c1 / r2 * c2
      } else {
        chimax = n * r2 * c2 / r1 * c1
      }
      chimin = n ^ 3 * (ahat - gahat) ^ 2 / (r1 * r2 * c1 * c2)
      es = (a * d - b * c) / abs(a * d - b * c) * sqrt((chi - chimin) / (chimax - chimin))
  } else if (method == "stiles"){
      es = log((a + b + c + d) * (abs(a * d - b * c) - (a + b + c + d) / 2) ^ 2 / ((a + b) * (a + c) * (b + d) * (c + d))) / log(10)
      
  } else if (method == "bonett-price-1"){
      w = ((a + 0.1) * (d + 0.1)) / ((b + 0.1) * (c + 0.1))
      pmin = min(r1, r2, c1, c2) / n
      x = 1 / 2 - (1 / 2 - pmin) ^ 2
      es = (w ^ x - 1) / (w ^ x + 1)
  
  } else if (method == "odds-ratio"){
        es = a * d / (b * c)
  } else if (method == "pearson-q4"){
        es = sin(pi / 2 * 1 / (1 + 2 * b * c * n / ((a * d - b * c) * (b + c))))
  } else if (method == "pearson-q5"){
      k = 4 * a * b * c * d * n ^ 2 / ((a * d - b * c) ^ 2 * (a + d) * (b + c))
      es = sin(pi / 2 * 1 / ((1 + k) ^ 0.5))
  } else if (method == "camp-1" || method == "camp-2" || method == "camp-3"){
      #create duplicates of original
      C1a = c1
      C2a = c2
      ca=a
      cb=b
      cc=c
      cd=d
      
      #step 1: if C1 < C2 swop the two columns
      switch=1
      if (c1 < c2){
        a2 = ca
        ca = cb
        cb = a2
        c2 = cc
        cc = cd
        cd = c2
        C1a = ca + cc
        C2a = cb + cd
        switch=-1
      }
      
      #step 2: determine three proportions
      p1 = ca/C1a
      p2 = cd/C2a
      p = C1a/n
      
      #step 3: determine the corresponding z-values
      z1 = qnorm(p1)
      z2 = qnorm(p2)
      
      #step 4: determine the height of the normal distribution
      y = dnorm(p)
      
      #step 5: calculate m
      m = p*(1-p)*(z1+z2)/y
      
      #step 6: find phi in table of phi-values
      if (method=="camp-1"){
        phi = 1
      }
      else if (method=="camp-2"){
        phiValues = data.frame("0.5" = 0.637, "0.6" = 0.63, "0.7" = 0.62, "0.8" = 0.6, "0.9" = 0.56)
        if (p<0.5){
          phi = phiValues[paste("X" ,round((1-p)*10)/10, sep="")]
        }
        else{
          phi = phiValues[paste("X" ,round(p*10)/10, sep="")]  
        }
        
        
      }
      else {
        phiValues = data.frame("0.5" = c(0.637, 0.636, 0.636, 0.635, 0.635, 0.634, 0.634, 0.633, 0.633, 0.632, 0.631), "0.6" = c(0.631, 0.631, 0.63, 0.629, 0.628, 0.627, 0.626, 0.625, 0.624, 0.622, 0.621), "0.7" = c(0.621, 0.62, 0.618, 0.616, 0.614, 0.612, 0.61, 0.608, 0.606, 0.603, 0.6), "0.8" = c(0.6, 0.597, 0.594, 0.591, 0.587, 0.583, 0.579, 0.574, 0.569, 0.564, 0.559))
        
        if (p<0.5){
          phi = phiValues[round((1-p)*100,0)-floor((1-p)*10)*10+1, paste("X" ,floor((1-p)*10)/10, sep="")]  
        }
        else {
          phi = phiValues[round(p*100,0)-floor(p*10)*10+1, paste("X" ,floor(p*10)/10, sep="")]  
        }
      }
      
      #step 7: calculate r
      es = switch * m/(1+phi*m**2)**0.5
      
  } else if (method=="becker-clogg-1" || method=="becker-clogg-2") {
      pR1 = r1/n
      pR2 = r2/n
      pC1 = c1/n
      pC2 = c2/n
      tr = qnorm(pR1)
      tc = qnorm(pC1)
      mR1 = -exp(-tr**2/2) / pR1
      mR2 = exp(-tr**2/2) / pR2
      vC1 = -exp(-tc**2/2) / pC1
      vC2 = exp(-tc**2/2) / pC2
      delta = (mR1 - mR2)*(vC1 - vC2)
      OR = a*d/(b*c)
      if (method=="becker-clogg-2") {
        es = (OR**(13.3/delta) - 1) / (OR**(13.3/delta) + 1)}
      else {
        phiBC = log(OR) / delta
        g = exp(12.4*phiBC - 24.6*phiBC**3)
        es = (g - 1)/(g + 1)}
  } else if(method=="bonett-price-2"){
      pMin = min(r1,r2, c1,c2)/n
      cBP = (1 - abs(r1 - c1)/(5*n) - (0.5 - pMin)**2)/2
      omg = a*d/(b*c) 
      es = cos(pi / (1 + omg**cBP))
  } else if(method=="bonett-price-3"){
    pMin2 = (min(r1,r2, c1,c2) + 1)/(n+2)
    cBP2 = (1 - abs(r1 - c1)/(5*(n+2)) - (0.5 - pMin2)**2)/2
    omg2 = (a+0.5)*(d+0.5)/((b+0.5)*(c+0.5)) 
    es = cos(pi  / (1 + omg2**cBP2))
  } else if(method == "chen-popovich"){
    p1 = r1 / n
    p2 = c1 / n
    z1 = qnorm(p1)
    z2 = qnorm(p2)
    es = (a * d - b * c) / (z1 * z2 * n ^ 2)
  }
  
  return (es)
  
  }




