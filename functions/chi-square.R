# calculate chi.sq test on object 'chi.square.data'
# Gene names in column 1, cell numbers in columns 2 and 3.
# clean function removes column 1 and any extra columns.  It returns data.frame then used for other analyses.
get.chisq.clean <- function(data){
  data <- data[, 2:3]
  new.data <- data
  counter <- 1
  mincells <- 1
  for (cell in rownames(data)){
    if (data[cell,1] < mincells & data[cell,2] < mincells){
      new.data <- new.data[-counter,]
    } else {
      counter <- counter + 1
    }
  }
  return(new.data)
}
get.chisq.pvalue <- function(data) {
  chi.square.analysis <- chisq.test(data)
  p.value <- chi.square.analysis$p.value
  message("p value is ",format(p.value,digits = 3))
  return(p.value)
}

get.chisq.df <- function(data){
  chi.square.analysis <- chisq.test(data)
  df <- chi.square.analysis$parameter
  message("df is ",df)
  return(df)
}

get.chisq.statistic <- function (data) {
  chi.square.analysis <- chisq.test(data)
  statistic <- chi.square.analysis$statistic
  message("Chi-square statistics is ", round(statistic,1))
  return(statistic)
}

get.chisq.expected <- function(data) {
  chi.square.analysis <- chisq.test(data)
  expected <- chi.square.analysis$expected
  return(expected)
}