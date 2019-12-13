# /*
# Software License Agreement (BSD License)
#
# Copyright (c) 2017, Matthias Kunz.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#  *
#  * * Redistributions of source code must retain the above copyright
#* notice, this list of conditions and the following disclaimer.
#* * Redistributions in binary form must reproduce the above
#* copyright notice, this list of conditions and the following
#* disclaimer in the documentation and/or other materials provided
#* with the distribution.
#* * Neither the name of Willow Garage, Inc. nor the names of its
#* contributors may be used to endorse or promote products derived
#* from this software without specific prior written permission.
#*
#  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
#* COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#                                                            * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#                                                            * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#* POSSIBILITY OF SUCH DAMAGE.
#*
#  */

## Function that plots a cylinder in three-dimensional space (x,y,z) using:
## the start point (x,y,z) and directional vector of the cylinder, it's radius and length
## BETA VERSION, 19.04.2017

#' Visualise a 3D cylinder
#' @author Matthias Kunz, last updated: 29.04.2019
#' @description Visualise a 3D cylinder
#' @param p_start A vector with the start point of the cylinder in 3D (x,y,z)
#' @param orient Directional vector of the cylinder axis (x1,y1,z1)
#' @param radius Radius of the cylinder (in meter)
#' @param length Length of the cylinder (in meter)
#' @param order Branch order level of the cylinder. This is optional for coloring. Default 0
#' @param n_seg Number of segments that draw the cylinder. Default 12
#' @param axis Display the cylinder axis. Default to FALSE.
#' @param cylinders Display the cylinder. Defaut to TRUE.
#' @param circles Display circles at the end of the cylinder. Default to FALSE.
#' @param transparent Makes the cylinder transparent. Default to TRUE with transparency of 0.3
#' @return Does not return anything at the moment.
#' @export
#' @examples
#' cylinder3d(p_start = c(0,0,0), orient = c(1,1,1), radius = 0.25, length = 1)
cylinder3d <- function(p_start = c(0,0,0), # Start point of the cylinder
                       orient = c(1,1,1),  # Directional vector of the cylinder axis
                       radius = 0.25,      # Radius of the cylinder
                       length = 1,         # Length of the cylinder
                       order = 0,          # Branch order level of the cylinder
                       n_seg = 12,         # Number of segments that draw the cylinder
                       axis = F,           # Display the cylinder axis
                       cylinders = T,      # Display the cylinder
                       circles = F,        # Display circles at the end of the cylinder
                       transparent = T)    # Makes the cylinder transparent
  {

  ## --------------------------------------------------------------------------
  ## Calculate end point of cylinder (p_end)
  p_end = c(0,0,0)
  Len = sqrt( (orient[1]^2) + (orient[2]^2) + (orient[3]^2) )
  dx = orient[1] / Len # Normalize direction vector if not already 1
  dy = orient[2] / Len # Normalize direction vector if not already 1
  dz = orient[3] / Len # Normalize direction vector if not already 1
  p_end[1] = p_start[1] + (length * dx)
  p_end[2] = p_start[2] + (length * dy)
  p_end[3] = p_start[3] + (length * dz)

  ## Check if distance between points and length match
  dist = sqrt((p_end[1]-p_start[1])^2 + (p_end[2]-p_start[2])^2 + (p_end[3]-p_start[3])^2 )

  ## Compute unit vector (must have length of 1)
  ## Note: p[1]->x, p[2]->y, p[3]->z
  orient_u = orient / sqrt(sum(orient^2))

  ## Get inclination/zenith (theta)
  theta_rad = acos(orient_u[3] / sqrt(sum(orient_u^2)))
  theta_deg = theta_rad * (180 / pi)

  ## Get azimuth (phi)
  phi_rad = atan2(orient_u[2], orient_u[1])
  phi_deg = phi_rad * (180 / pi)

  ## Set parameters / angle of circle (t)
  degvec <- seq(0, 2*pi, length=n_seg)

  ## Compute perpendicular and unit vectors
  u = c(-sin(phi_rad),
        cos(phi_rad),
        0)
  n = c(cos(phi_rad) * sin(theta_rad),
        sin(theta_rad) * sin(phi_rad),
        cos(theta_rad))
  nxu = c(cos(theta_rad) * cos(phi_rad),
          cos(theta_rad) * sin(phi_rad),
          -sin(theta_rad))

  ## Set color based on branch order
  col_code = data.frame(order=seq(from=0, to=9),
                        color=c("chocolate4","red","blue","darkgreen","darkorchid",
                                "deeppink1","black","black","black","black"))
  color = as.character(col_code[col_code == order,]$color)

  ## --------------------------------------------------------------------------
  ## Draw cylinder axis with start and end points
  if (axis == T) {
    rgl::segments3d(c(p_start[1], p_end[1]),
               c(p_start[2], p_end[2]),
               c(p_start[3], p_end[3]),
               color=color, lwd=2.5,
               add = T)
  }

  ## --------------------------------------------------------------------------
  ## Compute circle points at start and end of cylinder
  nn = n_seg+1
  cyl = data.frame(x_start=rep(NA, nn), y_start=rep(NA, nn), z_start=rep(NA, nn),
                   x_end=rep(NA, nn), y_end=rep(NA, nn), z_end=rep(NA, nn))

  ii = 1
  for (i in degvec){
    P_start = (radius * cos(i) * u) + (radius * sin(i) * nxu) + p_start
    P_end = (radius * cos(i) * u) + (radius * sin(i) * nxu) + p_end

    ## Store first point to add at the end twice for closing the cylinder
    if (ii == 1){
      P_first_start = P_start
      P_first_end = P_end
    }

    ## Store quadrilaterals information
    cyl$x_start[ii] = P_start[1]
    cyl$y_start[ii] = P_start[2]
    cyl$z_start[ii] = P_start[3]
    cyl$x_end[ii] = P_end[1]
    cyl$y_end[ii] = P_end[2]
    cyl$z_end[ii] = P_end[3]

    ii = ii + 1 ## Counter
  }

  ## Add first point at the end again for closing the cylinder
  cyl$x_start[nn] = P_first_start[1]
  cyl$y_start[nn] = P_first_start[2]
  cyl$z_start[nn] = P_first_start[3]
  cyl$x_end[nn] = P_first_end[1]
  cyl$y_end[nn] = P_first_end[2]
  cyl$z_end[nn] = P_first_end[3]

  ## --------------------------------------------------------------------------
  ## Draw cylinder
  # quads3d draws polygons with 4 corners, you have to specify the four corners
  # quads3d(c(1,1,0,0),c(0,0,0,0),c(-1,1,1,-1),col="green")
  if (cylinders == T) {
    for (i in 1:n_seg){
      if (transparent == F){alpha_value=1} else {alpha_value=0.3}
      rgl::quads3d(c(cyl$x_start[i], cyl$x_end[i], cyl$x_end[i+1], cyl$x_start[i+1]), # X values of quadrilaterals
              c(cyl$y_start[i], cyl$y_end[i], cyl$y_end[i+1], cyl$y_start[i+1]), # Y values of quadrilaterals
              c(cyl$z_start[i], cyl$z_end[i], cyl$z_end[i+1], cyl$z_start[i+1]), # Z values of quadrilaterals
              col=color, alpha = alpha_value, add=T)
    }
  }

  if (circles == T) {
    for (i in 1:n_seg-1){
      rgl::points3d(cyl$x_start[i], cyl$y_start[i], cyl$z_start[i], col=color, add=T)
    }
  }

  ## --------------------------------------------------------------------------
  ## Check calculation
  x = 1 * sin(theta_rad) * cos(phi_rad)
  y = 1 * sin(theta_rad) * sin(phi_rad)
  z = 1 * cos(theta_rad)
  #check <- c(x,y,z)
  #if(all.equal(check, orient_u)){cat("OK\n")}

  ## Print results
  #  cat("Cylinder base:", p_start, "\n")
  #  cat("Cylinder end:", p_end, "\n")
  #  cat("Orientation vector [x,y,z]:", orient,"\n")
  #  cat("Unit vector [x,y,z]:", orient_u,"\n")
  #  cat("Zenith [?]:",theta_deg)
  #  cat("\nAzimuth [?]:",phi_deg,"\n")
  #  cat(x,y,z)

  ## Show axes and labels in plot
  #axes3d(c('x', 'y', 'z'))
  #title3d('','','x','y','z')
}
