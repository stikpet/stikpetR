#' Number of Bins
#' 
#' @description 
#' To decide on the appropriate number of bins, many different rules can be applied. This function will 
#' determine the number of bins, based on the chosen method.
#' 
#' @param data vector or dataframe
#' @param method optional to indicate the method to use. Either "src", "sturges", "qr", "rice", "ts", "exp", "velleman", "doane", "scott", "fd", "shinshim", "stone", or "knuth"
#' @param maxBins optional for in iterations with "shinshim", "stone" and "knuth"
#' @param qmethod optional quartile method calculation to use for IQR when "fd" is used. See me_quartiles for options
#' 
#' @returns 
#' k : integer with optimum number of bins according to chosen method
#' 
#' @details
#' The first few methods are determining the number of bins (k) using the sample size (n).
#' 
#' **Square Root Choice (src)**
#' 
#' This method uses (unknown source):
#' \deqn{k = \lceil\sqrt{n}\rceil}
#' 
#' **Sturges Choice (sturges)**
#' 
#' This method uses (Sturges, 1926, p. 65):
#' \deqn{k = \lceil\log_2\left(n\right)\rceil + 1}
#' 
#' **Quartic Root (qr)**
#' 
#' This method uses (anonymous, as cited in Lohaka, 2007, p. 87):
#' \deqn{k = \lceil 2.5\times\sqrt[4]{n}\rceil}
#' 
#' **Rice Rule (rice)**
#' 
#' This method uses (Lane, n.d., p. 85):
#' \deqn{k = \lceil 2\times \sqrt[3]{n}\lceil}
#' 
#' **Terrell and Scotte (ts)**
#' 
#' This method uses (Terrell & Scott, 1985, p. 212):
#' \deqn{k = \lceil \sqrt[3]{2\times n\rceil}}
#' 
#' **Exponential (exp)**
#' 
#' This method uses (Iman & Conover, 1989, p. 54):
#' \deqn{k = \lceil \log_2\left(n\right)\rceil}
#' 
#' **Velleman (velleman)**
#' 
#' This method uses (Velleman, 1976 as cited in Lohaka, 2007, p. 89):
#' \deqn{k = \begin{cases}\lceil 2\times\sqrt{n}\rceil & \text{ if } n\leq 100 \\ \lceil 10\times\log_{10}\left(n\right)\rceil & \text{ if } n > 100\end{cases}}
#' 
#' 
#' **Doane (doane)**
#' 
#' This method uses (Doane, 1976, pp. 181-182):
#' \deqn{k = 1 + \lceil\log_2\left(n\right) + \log_2\left(1+\frac{\left|g_1\right|}{\sigma_{g_1}}\right)\rceil}
#' In the formula's \eqn{g_1} is the 3rd moment skewness:
#' \deqn{g_1 = \frac{\sum_{i=1}^n\left(x_i-\bar{x}\right)^3} {n\times\sigma^3} = \frac{1}{n}\times\sum_{i=1}^n\left(\frac{x_i-\bar{x}}{\sigma}\right)^3}
#' With:
#' \deqn{\sigma = \sqrt{\frac{\sum_{i=1}^n\left(x_i-\bar{x}\right)^2}{n}}}
#' The \eqn{\sigma_{g_1}} is defined using the formula:
#' \deqn{\sigma_{g_1}=\sqrt{\frac{6\times\left(n-2\right)}{\left(n+1\right)\left(n+3\right)}}}
#' 
#' Next are methods that determine the bin sizes (h), which can then be used to determine the number of bins (k) using:
#' \deqn{k = \lceil\frac{\text{max}\left(x\right)-\text{min}\left(x\right)}{h}\rceil}
#' 
#' **Scott (scott)**
#' 
#' This method uses (Scott, 1979, p. 608):
#' \deqn{h = \frac{3.49\times s}{\sqrt[3]{n}}}
#' 
#' Where \eqn{s} is the sample standard deviation:
#' \deqn{s = \sqrt{\frac{\sum_{i=1}^n\left(x_i-\bar{x}\right)^2}{n-1}}}
#' 
#' **Freedman and Diaconis (fd)**
#' 
#' This method uses (Freedman & Diaconis, 1981, p. 3):
#' \deqn{h = 2\times\frac{\text{IQR}\left(x\right)}{\sqrt[3]{n}}}
#' 
#' Where \eqn{\text{IQR}} is the inter-quartile range.
#' 
#' The last three methods all minimize a cost function (or maximize a profit function). They make use of the following steps:
#' 
#' 1. Divide the data into k bins and count the frequency in each bin
#' 1. Compute the cost function
#' 1. Repeat the first two steps while changing k, until a k is found that minimizes the cost function
#' 
#' **Shimazaki and Shinomoto (shinshim)**
#' 
#' This method uses as a cost function (Shimazaki & Shinomoto, 2007, p. 1508):
#' \deqn{C_k = \frac{2\times\bar{f_k}-\sigma_{f_k}}{h^2}}
#' With \eqn{\bar{f_k}} being the average of the frequencies when using k bins, and \eqn{\sigma_{f_k}} the population variance. 
#' In formula notation:
#' \deqn{\bar{f_k}=\frac{\sum_{i=1}^k f_{i,k}}{k}}
#' \deqn{\sigma_{f_k}=\frac{\sum_{i=1}^k\left(f_{i,k}-\bar{f_k}\right)^2}{k}}
#' Where \eqn{f_{i,k}} is the frequency of the i-th bin when using k bins.
#' 
#' **Stone (stone)**
#' 
#' This method uses as a cost function (Stone, 1984, p. 3):
#' \deqn{C_k = \frac{1}{h}\times\left(\frac{2}{n-1}-\frac{n+1}{n-1}\times\sum_{i=1}^k\left(\frac{f_i}{n}\right)^2\right)}
#' 
#' **Knuth (knuth)**
#' 
#' This method uses as a profit function (Knuth, 2019, p. 8):
#' 
#' \deqn{P_k=n\times\ln\left(k\right) + \ln\Gamma\left(\frac{k}{2}\right) - k\times\ln\Gamma\left(\frac{1}{2}\right) - \ln\Gamma\left(n+\frac{k}{2}\right) + \sum_{i=1}^k\ln\Gamma\left(f_i+\frac{1}{2}\right)}
#' 
#' @references
#' Doane, D. P. (1976). Aesthetic frequency classifications. *The American Statistician, 30*(4), 181–183. https://doi.org/10.2307/2683757
#' 
#' Freedman, D., & Diaconis, P. (1981). On the histogram as a density estimator. *Zeitschrift Für Wahrscheinlichkeitstheorie Und Verwandte Gebiete, 57*(4), 453–476. https://doi.org/10.1007/BF01025868
#' 
#' Iman, R. L., & Conover, W. J. (1989). *Modern business statistics* (2nd ed.). Wiley.
#' 
#' Knuth, K. H. (2019). Optimal data-based binning for histograms and histogram-based probability density models. *Digital Signal Processing, 95*, 1–30. https://doi.org/10.1016/j.dsp.2019.102581
#' 
#' Lohaka, H. O. (2007). Making a grouped-data frequency table: Development and examination of the iteration algorithm \[Doctoral dissertation, Ohio University\]. https://etd.ohiolink.edu
#' 
#' Scott, D. W. (1979). On optimal and data-based histograms. *Biometrika, 66*(3), 605–610. https://doi.org/10.1093/biomet/66.3.605
#' 
#' Shimazaki, H., & Shinomoto, S. (2007). A method for selecting the bin size of a time histogram. *Neural Computation, 19*(6), 1503–1527. https://doi.org/10.1162/neco.2007.19.6.1503
#' 
#' Stone, C. J. (1984). An asymptotically optimal window selection rule for kernel density estimates. *The Annals of Statistics, 12*(4), 1285–1297.
#' 
#' Sturges, H. A. (1926). The choice of a class interval. *Journal of the American Statistical Association, 21*(153), 65–66. https://doi.org/10.1080/01621459.1926.10502161
#' 
#' Terrell, G. R., & Scott, D. W. (1985). Oversmoothed nonparametric density estimates. *Journal of the American Statistical Association, 80*(389), 209–214. https://doi.org/10.2307/2288074
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples 
#' #Example 1
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' ex1 = df1['age']
#' ex1 = replace(ex1, ex1=="89 OR OLDER", "90")
#' tab_nbins(ex1)
#' 
#' #Example 2
#' file2 = 'https://peterstatistics.com/Packages/ExampleData/StudentStatistics.csv'
#' df2 = read.csv(file2, sep=';', na.strings=c("", "NA"))
#' ex2 = df2['Gen_Age']
#' tab_nbins(ex2)
#' 
#' @export
tab_nbins <- function(data, method='src', maxBins=NULL, qmethod="cdf"){
    data = data.frame(data)
    data <- na.omit(data) # remove missing values
    
    data <- as.numeric(data[,1]) # convert to numeric
    # Number of scores
    n <- length(data)
    
    if (is.null(maxBins)){
        maxBins=n
    }

    #Square-root choice
    if(method=='src'){k <- sqrt(n)}

    #Sturges
    else if(method=='sturges'){k <- log2(n) + 1}

    else if(method=='qr'){k <- 2.5*n^(1/4)}

    #Rice
    else if(method=='rice'){k <- 2*(n^(1/3))}

    #Terrell and Scott
    else if(method=='ts'){k <- (2*n)^(1/3)} 

    #Exponential
    else if(method=='exp'){k <- log2(n)}  

    #Exponential
    else if(method=='velleman'){
        if(n<=100){k <- 2*n^0.5}
        else{k <- 10*log10(n)}
    }  

    #Doane
    else if(method=='doane'){
        sdPop <- sqrt((n-1)/n)*sd(data)
        varPop <- sdPop^2
        g1 <- sum((data-mean(data)/varPop)^3)/n
        sg1 <- sqrt(6*(n-2)/((n+1)*(n+3)))
        k <- 1 + log2(n) + log2(1+abs(g1)/sg1)
    }

    else {
        dataRange <- max(data)-min(data)

        #Scott
        if(method=='scott'){
            h <- 3.49*sd(data)/(n^(1/3))
            k <- dataRange/h
        }

        #Freedman-Diaconis
        else if(method=='fd'){
            iqr = me_quartile_range(data, method = qmethod)[1,3]
            h <- 2*iqr/(n^(1/3))
            k <- dataRange/h
        }

        else {
            k <- 2: maxBins
            C <- numeric(length(k))
            h <- C

            for (i in 1:length(k)) {
                h[i] <- diff(range(data))/k[i]
                edges = seq(min(data),max(data),length=k[i])
                hp <- hist(data, breaks = edges, plot=FALSE )
                ki <- hp$counts

                avg <- mean(ki)
                v <- sum((ki-avg)^2)/k[i]

                if(method=="shinshim"){C[i] <- (2*avg-v)/(h[i]^2)}
                else if(method=="stone"){C[i] <- 1/h[i] * (2/(n-1) - (n+1)/(n-1)*sum((ki/n)^2))}
                else if(method=="knuth"){
                    C[i] <- n*log(k[i]) + lgamma(k[i]/2) - k[i]*lgamma(1/2) - lgamma(n+k[i]/2) + sum(lgamma(ki+0.5)) #Profit Function
                    C[i] <- -C[i]}
                    }
                idx <- which.min(C)
                h <- h[idx]
                k <- ceiling(dataRange/h)
            }
        }
    return (ceiling(k))
}