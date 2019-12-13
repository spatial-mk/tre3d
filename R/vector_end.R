#' Compute end point of given point in 3d with direction and length
#' @author Matthias Kunz, last updated: 29.04.2019
#' @description Compute end point of given point in 3d with direction and length
#' @param p_start A vector with the start point of the vector, e.g. a cylinder, in 3D (x,y,z)
#' @param orient Directional vector of the vector, e.g. a cylinder, axis (x1,y1,z1)
#' @param length Length of the vector, e.g. cylinder length (in meter)
#' @return The end point (p_end) of the vector in 3D (x,y,z)
#' @examples
#' vector_end(p_start = c(0,0,0), orient = c(1,1,1), length = 1)
vector_end <- function(p_start=c(0,0,0), orient=c(1,1,1), length=1){
  p_end = c(0,0,0)
  Len = sqrt( (orient[1]^2) + (orient[2]^2) + (orient[3]^2) )
  dx = orient[1] / Len # Normalize direction vector if not already 1
  dy = orient[2] / Len # Normalize direction vector if not already 1
  dz = orient[3] / Len # Normalize direction vector if not already 1
  p_end[1] = p_start[1] + (length * dx)
  p_end[2] = p_start[2] + (length * dy)
  p_end[3] = p_start[3] + (length * dz)
  return(p_end)
}
