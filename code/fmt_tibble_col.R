# format numeric columns in a kable before rendering
fmt_tibble_col <- function(df, col, sigfig, maxdigits, maxzeros) {
  nRows <- as.integer(count(df))
  for (iCol in col) {
    df[iCol] <- signif(df[iCol],sigfig)
    newCol <- c()
    for (iRow in 1:nRows) {
      curVal <- as.numeric(df[iRow,iCol])
      if (is.na(curVal)) {
        newRowVal <- ""
      } else if (abs(curVal) < 10^-maxzeros) {
        newRowVal <- format(curVal, scientific = TRUE)
      } else if (abs(curVal) < 10^maxdigits) {
        newRowVal <- format(curVal, scientific = FALSE, big.mark = ",")
      } else {
        newRowVal <- format(curVal, scientific = TRUE)
      }
      # add code to replace e or E with x 10^
      if (grepl('e',newRowVal,ignore.case = TRUE)) {
        newRowVal <-  sub('E\\+0',' x 10^',newRowVal,ignore.case = TRUE)
        newRowVal <-  sub('E-0',' x 10^-',newRowVal,ignore.case = TRUE)
        newRowVal <-  sub('E\\+',' x 10^',newRowVal,ignore.case = TRUE)
        newRowVal <-  sub('E-',' x 10^-',newRowVal,ignore.case = TRUE)
        newRowVal <- paste(newRowVal,"^",sep='')
      }
      newCol <- c(newCol,newRowVal)
      }
    df[iCol] <- newCol
  }
  df
}
