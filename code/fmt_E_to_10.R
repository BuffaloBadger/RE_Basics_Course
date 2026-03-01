# format numeric columns in a kable before rendering
fmt_E_to_10 <- function(value, sigfig, maxdigits, maxzeros) {
  value <- signif(value,sigfig)
  curVal <- as.numeric(value)
  if (is.na(curVal)) {
    newVal <- ""
  } else if (abs(curVal) == 0.0) {
    newVal <- "0.0"
  } else if (abs(curVal) < 10^-maxzeros) {
    newVal <- format(curVal, scientific = TRUE)
  } else if (abs(curVal) < 10^maxdigits) {
    newVal <- format(curVal, scientific = FALSE, big.mark = ",")
  } else {
    newVal <- format(curVal, scientific = TRUE)
  }
  # add code to replace e or E with x 10^
  if (grepl('e',newVal,ignore.case = TRUE)) {
    newVal <-  sub('E\\+0',' x 10^',newVal,ignore.case = TRUE)
    newVal <-  sub('E-0',' x 10^-',newVal,ignore.case = TRUE)
    newVal <-  sub('E\\+',' x 10^',newVal,ignore.case = TRUE)
    newVal <-  sub('E-',' x 10^-',newVal,ignore.case = TRUE)
    newVal <- paste(newVal,"^",sep='')
  }
  newVal
}
