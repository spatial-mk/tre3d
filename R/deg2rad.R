#' Degree to radian conversion.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Converts an angle in degree (0-360) to radian (0-1)
#' @param angle A positive angle in degree (0-360)
#' @return The angle in radian
#' @examples
#' deg2rad(45.0) # 0.7853982
deg2rad <- function(angle=0){
  angle = angle * pi / 180
  return(angle)
}
