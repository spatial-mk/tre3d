## ****************************************************************************
## Function that visualises a voxel at a given local cartesian coordinate
## Matthias Kunz, 18.02.2020

#' Visualise 3d point as voxel
#' @author Matthias Kunz, last updated: 18.02.2020
#' @description Visualise 3d point as voxel. Usage with rgl: open3d(); cube3d(1,1,1,...)
#' @param x x coordinate
#' @param y y coordinate
#' @param z z coordinate
#' @param scale Size of the cube, e.g. 0.1 when 10 cm voxels are used. Default 1 (meter).
#' @param bordered Logical, Draws a border
#' @param filled Logical, Fills the cube with a color
#' @param lwd Line width of the border. Default 2
#' @param bordercol Color of the border
#' @param alpha Transparancy of the cube
#' @param fillcol ## Color of the cube. Default gray.
#' @return A cube object that can be plotted with rgl
#' @export
#' @examples
#' xyz2cube3d(data.frame(x=runif(20),y=runif(20), z=runif(20)), alpha=1.0)
## ****************************************************************************

## Additional libraries
#library(rgl) ## Plotting in 3D
#library(magrittr) ## for use of %>%
#library(scales) ## for transparent colors in R

## Function
xyz2cube3d <- function(x=0,y=0,z=0, ## Coordinates of the cube
                      bordered=TRUE, ## Logical, Draws a border
                      filled = TRUE, ## Logical, Fills the cube with a color
                      lwd=2, ## Line width of the border
                      scale=1, ## Size of the cube, e.g. 0.1 when 10 cm voxels are used
                      fillcol = "gray", ## Color of the cube
                      bordercol ='black', ## Color of the border
                      alpha=0.5, ## Transparancy of the cube
                      ...) {
  mycube <- rgl::cube3d(color=fillcol, alpha=alpha)

  # Reduce size to unit
  mycube$vb[4,] <- mycube$vb[4,]/scale*2

  for (i in 1:length(x)){
    # Add cube border
    if (bordered) {
      bcube <- mycube
      bcube$material$lwd <- lwd
      bcube$material$front <- 'line'
      bcube$material$back <- 'line'
      bcube %>% translate3d(x[i], y[i], z[i]) %>% shade3d
    }

    # Add cube fill
    if (filled) {
      fcube <- mycube
      fcube$vb[4,] <- fcube$vb[4,]*1.01
      fcube$material$col <- alpha(fillcol, alpha)
      fcube %>% translate3d(x[i], y[i], z[i]) %>% shade3d
    }
  }
}
