#' Copy A table in Excel Format
#'

copyxl <- function(df, sep="\t", dec=".", max.size=(200*1000)){
  Title <- deparse(substitute(df))
  df2 <- list(Title,df)
  write.table(df2, paste0("clipboard-", formatC(max.size, format="f", digits=3)), sep=sep, col.names=TRUE, row.names=TRUE, dec=dec)
}
