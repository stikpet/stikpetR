me_variation <- function(data, 
                         levels=NULL, 
                         measure="std", 
                         ddof=1, 
                         center="mean", 
                         azs="square"){
  
  if (is.null(levels)){dataN = data}
  else{
    myFieldOrd = factor(na.omit(data), ordered = TRUE, levels = levels)
    dataN = as.numeric(myFieldOrd)
  }
  
  dataN = na.omit(dataN)
  n = length(dataN)
  if (measure=="std"){
    if (ddof==1){lbl = "standard deviation (sample)"}
    else if (ddof==0){lbl = "standard deviation (population)"}
    else {lbl=paste0("standard deviation corrected with ", ddof)}
    m = mean(dataN)
    res = (sum((dataN - m)**2)/(n-ddof))**0.5}
  else if (measure=="var"){
    if (ddof==1){lbl = "variance (sample)"}
    else if (ddof==0){lbl = "variance (population)"}
    else {lbl=paste0("variance corrected with ", ddof)}
    m = mean(dataN)
    res = sum((dataN - m)**2)/(n-ddof)}
  else if (measure=="mad"){
    lbl = "mean absolute deviation"
    m = mean(dataN)
    res = sum(abs(dataN - m))/n}
  else if (measure=="madmed"){
    lbl = "mean absolute deviation around median"
    m = median(dataN)
    res = sum(abs(dataN - m))/n}
  else if (measure=="medmad"){
    lbl = "median absolute deviation"
    m = median(dataN)
    res = median(abs(dataN - m))}
  else if (measure=="cv"){
    lbl = "coefficient of variation"
    mu = mean(dataN)
    s = me_variation(dataN, measure="std", ddof=ddof)
    res = s/mu}
  else if (measure=="stddm"){
    lbl = "standard deviation with decile mean"
    mu = me_mean(dataN, version="decile")
    res = (sum((dataN - mu)**2)/(n-ddof))**0.5}
  else if (measure=="cd"){
    lbl = "coefficient of deviation"
    mu = me_mean(dataN, version="decile")
    s = (sum((dataN - mu)**2)/(n-ddof))**0.5
    res = s/mu}
  else{
    if (center=="mean"){
      lbl = "mean"
      mu = mean(dataN)}
    else if (center=="median"){
      lbl = "median"
      mu = median(dataN)}
    else if (center=="mode"){
      lbl = "mode"
      mu = me_mode(dataN)[1,1]}
    else {
      lbl = center
      mu = center}
    
    if (azs=="square"){
      lbl = paste("sum squared deviation around ", lbl)
      res = sum((dataN - mu)**2)}
    else if (azs=="abs"){
      lbl = paste("sum absolute deviation around ", lbl)
      res = sum(abs(dataN - mu))
    }
  }
  results = data.frame(res, lbl)
  colnames(results)<-c("value", "measure")    
  return (res)
}