#' Tetrachoric Correlation Coefficient
#' @importFrom fMultivar pnorm2d
#' 
#' @description 
#' In essence this attempts to mimic a correlation coefficient between two scale variables. It can be defined as "An estimate of the correlation between two random variables having a bivariate normal distribution, obtained from the information from a double dichotomy of their bivariate distribution" (Everitt, 2004, p. 372). 
#' 
#' This assumes the two binary variables have ‘hidden’ underlying normal distribution. If so, the combination of the two forms a bivariate normal distribution with a specific correlation between them. The quest is then to find the correlation, such that the cumulative density function of the z-values of the two marginal totals of the top-left cell (a) match that value.
#' 
#' This is quite tricky to do, so a few have proposed an approximation for this. These include Yule r, Pearson Q4 and Q5, Camp, Becker and Clogg, and Bonett and Price , all available and more with (\code{\link{es_bin_bin}}).
#' 
#' Besides closed form approximation formula's, various algorithms have been designed as well. The three most often mentioned are Brown (1977), Kirk (1973), and Divgi (1979), available in this function.
#' 
#' @param field1 : dataframe field with categories for the rows
#' @param field2 : dataframe field with categories for the columns
#' @param categories1 : optional list with selection and/or order for categories of field1
#' @param categories2 : optional list with selection and/or order for categories of field2
#' @param method method to use (see details). Either "divgi" (default), "search", "kirk", "brown"
#' 
#' @return Tetrachoric Correlation Coefficient
#' 
#' @details 
#' The "search" method does a binary search for rt using the bivariate normal distribution
#' from the *fMultivar* library.
#' 
#' "kirk" will use Kirk (1973) Fortran TET8 procedure, adapted by stikpet
#' 
#' "brown" will use Brown (1977) - Algorithm AS 116
#' 
#' "divgi" will use Divgi (1979) algorithm
#' 
#' Flow charts of these algorithms can be found at https://peterstatistics.com
#' 
#' @references 
#' Brown, M. B. (1977). Algorithm AS 116: The tetrachoric correlation and its asymptotic standard error. *Applied Statistics, 26*(3), 343. https://doi.org/10.2307/2346985
#' 
#' Divgi, D. R. (1979). Calculation of the tetrachoric correlation coefficient. *Psychometrika, 44*(2), 169–172. https://doi.org/10.1007/BF02293968
#' 
#' Kirk, D. B. (1973). On the numerical approximation of the bivariate normal (tetrachoric) correlation coefficient. *Psychometrika, 38*(2), 259–268. https://doi.org/10.1007/BF02291118
#' 
#' @author 
#' P. Stikker. [Companion Website](https://PeterStatistics.com), [YouTube Channel](https://www.youtube.com/stikpet), [Patreon donations](https://www.patreon.com/bePatron?u=19398076)
#' 
#' @examples
#' #Example: dataframe
#' dataFile = "https://peterstatistics.com/Packages/ExampleData/GSS2012a.csv"
#' df1 <- read.csv(dataFile, sep=",", na.strings=c("", "NA"))
#' r_tetrachoric(df1[['mar1']], df1[['sex']], categories1=c("WIDOWED", "DIVORCED"))
#' 
#' @export
r_tetrachoric <- function(field1, field2, categories1=NULL, categories2=NULL, method="divgi"){
  
  #Create a cross table first
  ct = tab_cross(field1, field2, order1=categories1, order2=categories2)
  
  #store the individual cells
  a = ct[1,1]
  b = ct[1,2]
  c = ct[2,1]
  d = ct[2,2]
  
  if (method=="search"){
    rt = r_tetrachoric_search(a, b, c, d)
  }
  
  else if (method=="kirk"){
    rt = r_tetrachoric_kirt(a, b, c, d)
  }
  
  else if (method=="brown"){
    rt = r_tetrachoric_brown(a, b, c, d)
  }
  
  else if (method=="divgi"){
    rt = r_tetrachoric_divgi(a, b, c, d)
  }
  
  return(rt)
  
}

r_tetrachoric_search <- function(a,b,c,d){
  
  #the row totals
  R1 <- a+b
  R2 <- c+d
  
  #the column totals
  C1 <- a+c
  C2 <- b+d
  
  n = a+b+c+d
  
  p1 <- R1/n
  p2 <- C2/n
  p <- a/n
  
  z1 <- qnorm(p1)
  z2 <- -qnorm(p2)
  
  #Iterate to find optimal value for rt
  
  nDecimals <- 10
  
  rt = -1
  for (nd in 1:nDecimals) {
    prt<-0
    i <- -1
    while (prt < p && i < 1) {
      rt <- rt + 1/(10**nd)
      prt <- fMultivar::pnorm2d(z1, z2, rt)[1]
      i <- i + 0.1
    }
    rt <- rt - 1/(10**nd)
  }
  
  return (rt)
}

r_tetrachoric_divgi <- function(a,b,c,d){
  R1 = a+b
  C1 = a+c
  n = a+b+c+d
  
  p1 <- R1/n
  p2 <- C1/n
  
  h <- qnorm(p1)
  k <- qnorm(p2)
  
  hAdj <- max(abs(h), abs(k))
  kAdj <- min(abs(h), abs(k))
  
  sg <- sign(h)*sign(k)
  
  #initial guess = Odds Ratio
  OR <- a*d/(b*c)
  
  dA <- 0.5/(1+(hAdj**2+kAdj**2)*(0.12454-0.27102*(1-hAdj/(sqrt(hAdj**2+kAdj**2)))))
  dB <- 0.5/(1+(hAdj**2+kAdj**2)*(0.82281-1.03514*kAdj/(sqrt(hAdj**2+kAdj**2)))) 
  dC <- 0.07557*hAdj+(hAdj-kAdj)**2*(0.51141/(hAdj+2.05793)-0.07557/hAdj)
  dD <- kAdj*(0.79289+4.28981/(1+3.30231*hAdj))
  
  alp <- dA + dB*(-1 + 1/(1 + dC*(log(OR)-dD)**2))
  
  r <- cos(pi/(1+OR**alp))
  
  for (i in 1:10) {
    L <- fMultivar::pnorm2d(h, k, r)[1]
    Ld <- exp(-(h**2-2*r*h*k+k**2)/(2*(1-r**2)))/(2*pi*sqrt(1-r**2))
    r <- r - (L - a/n)/Ld
  }
  r <- r*sg
  
  rt= r
  
  return (rt)
  
}

r_tetrachoric_kirt <- function(a, b, c, d){
    # P. Stikker adaptation of Fortran IV code from:
    # Kirk, D. B. (1973). On the numerical approximation of the bivariate normal (tetrachoric) correlation coefficient. Psychometrika, 38(2), 259–268. doi: 10.1007/BF02291118
    af = a
    bf = b
    cf = c
    df = d
    # Since:
    # P = the joint proportion for which both marginal values are < .5
    # We might have to swop the rows and/or columns
    
    
    if (af + bf > cf + df) {
      # swop rows
      oldTemp = af
      af = cf
      cf = oldTemp
      oldTemp = bf
      bf = df
      df = oldTemp
    }
    if (af + cf > bf + df) {
      # swop columns
      oldTemp = af
      af = bf
      bf = oldTemp
      oldTemp = cf
      cf = df
      df = oldTemp
    }
    # input parameters of original function
    n = af + bf + cf + df
    p = af / n
    fm1 = (af + bf) / n
    fm2 = (af + cf) / n
    # value settings
    c1 = 0.019855071751232
    c2 = 0.101666761293187
    c3 = 0.237233795041835
    c4 = 0.408282678752175
    c5 = 0.591717321247825
    c6 = 0.762766204958165
    c7 = 0.898333238706813
    c8 = 0.980144928248768
    c1sq = 0.000394223874247
    c2sq = 0.010336130351846
    c3sq = 0.056279873509951
    c4sq = 0.166694745769052
    c5sq = 0.350129388264702
    c6sq = 0.581812283426281
    c7sq = 0.807002607765472
    c8sq = 0.960684080371783
    w1 = 0.050614268145188
    w2 = 0.111190517226687
    w3 = 0.156853322938943
    w4 = 0.181341891689181
    rpi = 0.398942280401433
    rtpi = 2.506628274631
    r = 0
    
    # EVALUATION OF H AND K
    j = 1
    
    # CONVERGENCE CRITERION FOR H AND K
    eps = 1 * 10^-5
    x = fm1
    con = (0.5 - fm1) * rtpi
    
    #GO TO 200
    loop200 = TRUE;
    while (loop200) {
      loop200 = FALSE;
      if (x > 0.5) {
        # GO TO 600
        # 600
        r = 4;
        # GO TO 1000
      } 
      else {
        # HASTINGS' APPROXIMATION
        # 250
        e = sqrt(-2 * log(x));
        e1 = (0.010328 * e + 0.802853) * e + 2.515517;
        e2 = ((0.001308 * e + 0.189269) * e + 1.432788) * e + 1;
        est = e - e1 / e2;
        loop300 = TRUE;
        
        while (loop300) {
          loop300 = FALSE;
          # 300
          old = est;
          # ITERATION LOOP
          loop350 = TRUE;
          i = 1;
          
          while (loop350 && i <= 20) {
            loop350 = FALSE;
            # 350
            xx = old * old;
            if (j == 3 && abs(old) >= 1) {
              # GO TO 550
              # 550
              if (r == 2) {
                # GO TO 1000
              } 
              else {
                r = 2;
                old = 0.97 * sign(a);
                # GO TO 350
                i = i + 1;
                loop350 = TRUE;
              }
            } 
            else {
              if (j == 3) {
                #GO TO 400
                # 400
                # R INTEGRAND EVALUATION
                fnum = old * (w1 * (hEfn2(h2k2, hk2, x, xx, c1, c1sq) + hEfn2(h2k2, hk2, x, xx, c8, c8sq)) + w2 * (hEfn2(h2k2, hk2, x, xx, c2, c2sq) + hEfn2(h2k2, hk2, x, xx, c7, c7sq)) + w3 * (hEfn2(h2k2, hk2, x, xx, c3, c3sq) + hEfn2(h2k2, hk2, x, xx, c6, c6sq)) + w4 * (hEfn2(h2k2, hk2, x, xx, c4, c4sq) + hEfn2(h2k2, hk2, x, xx, c5, c5sq)));
                fnew = old - (fnum - con) / hEfn2(h2k2, hk2, x, xx, 1, 1);
              } 
              else {
                # H AND K INTEGRAND EVALUATION FOR GAUSSIAN QUADRATURE
                fnum = old * (w1 * (hEfn1(xx, c1sq) + hEfn1(xx, c8sq)) + w2 * (hEfn1(xx, c2sq) + hEfn1(xx, c7sq)) + w3 * (hEfn1(xx, c3sq) + hEfn1(xx, c6sq)) + w4 * (hEfn1(xx, c4sq) + hEfn1(xx, c5sq)));
                fnew = old - (fnum - con) / hEfn1(xx, 1);
                # GO TO 450
              }
              
              # TEST FOR CONVERGENCE
              #450
              if (abs(old - fnew) > eps) {
                # GO TO 500
                # 500
                old = fnew;
                i = i + 1;
                if (i > 20) {
                  r = 3;
                  # GO TO 1000
                } 
                else {
                  loop350 = TRUE;
                }
              } 
              else {
                if (j == 1) {
                  # GO TO 100
                  # 100
                  j = 2;
                  h = fnew;
                  h2 = h * h;
                  zh = rpi * exp(-1 * h2 / 2);
                  x = fm2;
                  con = (0.5 - fm2) * rtpi;
                  
                  # GO TO 200
                  loop200 = TRUE;
                } 
                else {
                  if (j == 2) {
                    # GO TO 150
                    # 150
                    fk = fnew;
                    fk2 = fk * fk;
                    h2k2 = h2 + fk2;
                    hk2 = h * fk * 2;
                    zk = rpi * exp(-1 * fk2 / 2);
                    j = 3;
                    
                    # STARTING ESTIMATE FOR R
                    a = p - fm1 * fm2;
                    con = a * 6.28318530717959;
                    zhk = zh * zk;
                    est = 2 * (sqrt(abs(hk2 * a / zhk + 1)) - 1) / hk2;
                    
                    if (abs(hk2) <= 1 * 10^-8) {
                      est = a / zhk;
                    }
                    if (abs(est) > 0.8) {
                      est = 0.8 * sign(a);
                    }
                    
                    # CONVERGENCE CRITERION FOR R
                    eps = 1 * 10^-4;
                    # GO TO 300
                    loop300 = TRUE;
                  } 
                  else {
                    # GO TO 650
                    r = fnew;
                    # GO TO 1000
                  }
                }
              }
            }
          }
        }
      }
    }
    # 200
    # 1000
    rt = r
    
    return (rt)
  
  
}


r_tetrachoric_brown <- function(a, b, c, d){
    x = c(0.9972638618, 0.9856115115, 0.9647622556, 0.9349060759, 
          0.8963211558, 0.8493676137, 0.794483796, 0.7321821187, 
          0.6630442669, 0.5877157572, 0.5068999089, 0.4213512761, 
          0.3318686023, 0.2392873623, 0.1444719616, 0.0483076657)
    w = c(0.00701861, 0.0162743947, 0.0253920653, 0.0342738629,
          0.042835898, 0.0509980593, 0.0586840935, 0.0658222228, 
          0.0723457941, 0.0781938958, 0.0833119242, 0.087652093, 
          0.0911738787, 0.0938443991, 0.0956387201, 0.096540085)
    sqt2pi = 2.50662827;
    rlimit = 0.9999;
    rcut = 0.95;
    uplim = 5;
    var_const = 1 * 10^(-60);
    chalf = 1 * 10^(-30);
    conv = 1 * 10^-8;
    citer = 1 * 10^-6;
    niter = 25;
    
    #INITIALIZATION
    r = 0;
    itype = 0;
    
    #CHECK IF ANY CELL FREQUENCY IS NEGATIVE
    if (a < 0 || b < 0 || c < 0 || d < 0) {
      print("error in input");
      r = 2;
    } 
    else {
      #CHECK IF ANY FREQUENCY IS ZERO AND SET KDELTA
      kdelta = 1;
      delta = 0;
      if (a == 0 || d == 0) {
        kdelta = 2;
      }
      if (b == 0 || c == 0) {
        kdelta = kdelta + 2;
      }
      
      #KDELTA=4 MEANS TABLE HAS ZERO ROW OR COLUMN, RUN IS TERMINATED
      #DELTA IS 0.0, 0.5 OR -0.5 ACCORDING TO WHICH CELL IS ZERO
      if (kdelta == 4) {
        #GO TO 92
        print("error row or column total is zero");
        r = 2;
      } 
      else {
        if (kdelta == 1) {
          #GOTO 4
        } 
        else {
          if (kdelta == 2) {
            #GOTO 1
            #1
            delta = 0.5;
            if (a == 0 && d == 0) {
              r = -1;
            }
            #GOTO 4
          } 
          else {
            if (kdelta == 3) {
              #GOTO 2
              #2
              delta = -0.5;
              if (b == 0 && c == 0) {
                r = 1;
              }
            }
          }
        }
        
        #4
        if (r != 0) {
          itype = 3;
        }
        
        #STORE FREQUENCIES IN AA, BB, CC AND DD
        aa = a + delta;
        bb = b - delta;
        cc = c - delta;
        dd = d + delta;
        tot = aa + bb + cc + dd;
        
        #CHECK IF CORREIATION IS NEGATIVE, ZERO, POSITIVE
        #COMPUTE PROBILITIES OF QUADRANT AND OF MARGINALS PRBMAA AND PROBAC CHOSEN SO THAT CORREIATION IS POSITIVE.  KSIGN INDICATES WHETHER QUADRANTS HAVE BEEN SWITCHED
        if (aa * dd - bb * cc < 0) {
          #GO TO 7
          #7
          probaa = bb / tot;
          probac = (bb + dd) / tot;
          ksign = 2;
        } 
        else {
          if (aa * dd - bb * cc == 0) {
            #GO TO 5
            #5
            itype = 4;
          }
          
          #6
          probaa = aa / tot;
          probac = (aa + cc) / tot;
          ksign = 1;
          
          #GO TO 8
        }
        
        #8
        probab = (aa + bb) / tot;
        
        #COMPUTE NORMALDEVIATES FOR THE MARGINAL FREQUENCIES
        #SINCE NO MARGINAL CAN BE ZERO, IE IS NOT CHECKED
        zac = qnorm(probac);
        zab = qnorm(probab);
        ss = exp(-0.5 * (zac^2 + zab^2)) / (2 * pi);
        
        #WHEN R IS 0.0, 1.0 OR -1.0, TRANSFER TO COMUTE SDZERO
        if (r != 0 || itype > 0) {
          #GOTO 85
          skip85 = FALSE;
        } 
        else {
          #WHEN MARGINALS ARE EQUAL, COSINE EVALUATION IS USED
          if (a == d && b == c) {
            #GOTO 60
            #60
            rr = -1 * cos(2 * pi * probaa);
            itype = 2;
          } 
          else {
            #INITIAL ESTIMATE OF CORRELATION IS YULES Y
            rr = (sqrt(aa * dd) - sqrt(bb * cc))^2 / abs(aa * dd - bb * cc);
            iter = 0;
            
            #IF RR EXCEEDS RCUT, GAUSSIANQUADRATURE IS USED
            loop10 = TRUE;
            while (loop10) {
              loop10 = FALSE;
              #10
              if (rr > rcut) {
                #GOTO 40
                skip40 = FALSE;
              } 
              else {
                
                #TETRACHORIC SERIES IS COMPUTED
                
                #INITIALIALIZATION
                va = 1;
                vb = zac;
                wa = 1;
                wb = zab;
                term = 1;
                iterm = 0;
                sum = probab * probac;
                deriv = 0;
                sr = ss;
                
                #15
                loop15 = TRUE;
                
                while (loop15) {
                  loop15 = FALSE;
                  if (abs(sr) > var_const) {
                    #GOTO 20
                  } 
                  else {
                    #RESCALE TERMS TO AVOID OVERFLOWS AND UNDERFLOWS
                    sr = sr / var_const;
                    va = va * chalf;
                    vb = vb * chalf;
                    wa = wa * chalf;
                    wb = wb * chalf;
                  }
                  
                  #FORM SUM AND DERIV OF SERIES
                  #20
                  dr = sr * va * wa;
                  sr = sr * rr / term;
                  cof = sr * va * wa;
                  
                  #ITERM COUNTS NO. OF CONSECUTIVE TERMS .LT. CONV
                  iterm = iterm + 1;
                  if (abs(cof) > conv) {
                    iterm = 0;
                  }
                  sum = sum + cof;
                  deriv = deriv + dr;
                  vaa = va;
                  waa = wa;
                  va = vb;
                  wa = wb;
                  vb = zac * va - term * vaa;
                  wb = zab * wa - term * waa;
                  term = term + 1;
                  if (iterm < 2 || term < 6) {
                    loop15 = TRUE;
                  }
                }
                #CHECK IF ITERATION CONVERGED
                if (abs(sum - probaa) > citer) {
                  #GOTO 25
                  #CALCULATE NEXT ESTIMATE OF CORRELATION
                  #25
                  iter = iter + 1;
                  
                  #IF TOO MANY ITERATIONS, RUN IS TERMINATED
                  if (iter >= niter) {
                    #GOTO 93
                    skip40 = TRUE;
                    skip70 = TRUE;
                    skip85 = FALSE;
                  } 
                  else {
                    delta = (sum - probaa) / deriv;
                    rrprev = rr;
                    rr = rr - delta;
                    if (iter == 1) {
                      rr = rr + 0.5 * delta;
                    }
                    if (rr > rlimit) {
                      rr = rlimit;
                    }
                    if (rr < 0) {
                      rr = 0;
                    }
                    
                    #GOTO 10
                    loop10 = TRUE;
                    skip40 = FALSE;
                  }
                } 
                else {
                  #ITERATION HAS CONVERGED,SET ITYPE
                  itype = term;
                  #GOTO 70
                  skip40 = TRUE;
                  skip70 = FALSE;
                }
              }
            }
            if (!skip40) {
              #GAUSSIAN QUADRATURE
              #40
              if (iter > 0) {
                #GOTO 41
              } 
              else {
                #INITIALIZATION IF THIS IS FIRST ITERATION
                sum = probab * probac;
                rrprev = 0;
              }
              
              #INITIALIZATION
              #41
              sumprv = probab - sum;
              prob = bb / tot;
              if (ksign == 2) {
                prob = aa / tot;
              }
              itype = 1;
              
              #LOOP TO FIND ESTIMATE OF CORRELATION
              #COMPUTATION OF INTEGRAL(SUM) BY QUADRATURE
              loop42 = TRUE;
              while (loop42) {
                loop42 = FALSE;
                
                #42
                rrsq = sqrt(1 - rr^2);
                amid = 0.5 * (uplim + zac);
                xlen = uplim - amid;
                sum = 0;
                for (iquad in range(1, 17)) {
                  xla = amid + x[iquad] * xlen;
                  xlb = amid - x[iquad] * xlen;
                  
                  #TO AVOID UNDERFLOWS, TEMPA AND TEMPB ARE USED
                  tempa = (zab - rr * xla) / rrsq;
                  if (tempa >= -6) {
                    sum = sum + w[iquad] * exp(-0.5 * xla^2) * pnorm(tempa);
                  }
                  tempb = (zab - rr * xlb) / rrsq;
                  if (tempb >= -6) {
                    sum = sum + w[iquad] * exp(-0.5 * xlb^2) * pnorm(tempb);
                  }
                }
                #44
                sum = sum * xlen / sqt2pi;
                
                #CHECX IF ITERATION HAS CONVERGED
                if (abs(prob - sum) <= citer) {
                  #GOTO 70
                  skip70 = FALSE;
                }
                else {
                  iter = iter + 1;
                }
                
                #IF TOO MANY ITERATIONS, RUN IS TERMINATED
                if (iter >= niter) {
                  #GOTO 93
                  skip70 = TRUE;
                  skip85 = FALSE;
                }
                
                #ESTIMATE CORRELATIONFOR NEXT ITERATION BY LINEAR INTERPOLATION
                rrest = ((prob - sum) * rrprev - (prob - sumprv) * rr) / (sumprv - sum);
                
                #IS ESTIMATE POSITIVE AND LESS THAN UPPER LIMIT
                if (rrest > rlimit) {
                  rrest = rlimit;
                }
                if (rrest < 0) {
                  rrest = 0;
                }
                rrprev = rr;
                rr = rrest;
                sumprv = sum;
                
                #IF ESTIMATE HAS SAME VALUE ON TWO ITERATIONS, STOP ITERATION
                if (rr == rrprev) {
                  #GOTO 70
                  skip70 = FALSE;
                }
                else {
                  #GOTO 42
                  loop42 = TRUE;
                }
              }
            }
          }
          if (!skip70) {
            #COMPUTE SDR
            #70
            r = rr;
            rrsq = sqrt(1 - r^2);
            if (kdelta > 1) {
              itype = -itype;
            }
            if (ksign == 1) {
              #GOTO 71
            } 
            else {
              r = -r;
              zac = -zac;
            }
            
            #71
            pdf = exp(-0.5 * (zac^2 - 2 * r * zac * zab + zab^2) / rrsq^2) / (2 * pi * rrsq);
            pac = pnorm((zac - r * zab) / rrsq) - 0.5;
            pab = pnorm((zab - r * zac) / rrsq) - 0.5;
            sdr = (aa + dd) * (bb + cc) / 4 + pab^2 * (aa + cc) * (bb + dd) + pac^2 * (aa + bb) * (cc + dd) + 2 * pab * pac * (aa * dd - bb * cc) - pab * (aa * bb - cc * dd) - pac * (aa * cc - bb * dd);
            if (sdr < 0) {
              sdr = 0;
            }
            sdr = sqrt(sdr) / (tot * pdf * sqrt(tot));
            skip85 = FALSE;
          }
        }
        if (!skip85) {
          #85
          sdzero = sqrt((aa + bb) * (aa + cc) * (bb + dd) * (cc + dd) / tot) / (tot^2 * ss);
          if (r == 0) {
            sdr = 0;
          }
        }
      }
    }
    rt = r
    
    return (rt)
}

hEfn1 <- function(xx, w) {
  f = exp(-1 * xx * w / 2);
  return (f)
}

hEfn2 <- function(h2k2, hk2, x, xx, v, w) {
  f = exp((-1 * h2k2 + hk2 * v * x) / (2 * (1 - w * xx))) / sqrt(1 - w * xx);
  return (f)
}

