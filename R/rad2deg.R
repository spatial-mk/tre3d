#' Radian to degree conversion.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @describeIn Converts an angle in radian (0-1) to degree (0-360)
#' @param angle A positive angle in radian (0-1)
#' @return The angle in degree (0-360)
#' @examples
#' rad2deg(0.7853982) # 45.0
rad2deg <- function(angle=0){

  angle = angle * 180 / pi
  return(angle)

} # END-OF-FUNCTION
