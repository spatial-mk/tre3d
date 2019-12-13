#' Compute absolute difference in percent between two values
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description  Compute absolute difference in percent between two values
#' @param x1 A value
#' @param x2 Another value
#' @return Difference between x1 and x2 in absolute (always positive) percent
#' @export
#' @examples
#' diff_percent(10,12) # 18.18%
diff_percent <- function(x1,x2){
  abs_diff <- abs(x1-x2)
  mean <- mean(c(x1,x2), na.rm = T)
  diff = (abs_diff/mean) * 100
  return(diff)
}
