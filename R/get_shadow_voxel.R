#' Simple identification of shadow voxel inside the stem segment of the tree voxelmodel.
#' @author Carsten Hess, last updated: 14.02.2020
#' @description If the voxel size is much smaller than DBH a simple voxelization will not reconsider the space inside the stem.
#' @param input_vox Voxelmodel of the trees stem segment (X,Y,Z)
#' @return Voxelmodel with added and marked (TRUE, FALSE) shadow voxel.
#' @export
#' @examples
#'

get_shadow_voxel <- function(input_vox=NA, voxel_size=0.2) {

  # Set colnames to X Y Z
  colnames(input_vox) <- c("X","Y","Z")

  input_vox$shadow <- FALSE

  # Check for occluded shadow voxel
  top <- max(input_vox$Z)
  bottom <- min(input_vox$Z)

  # table for potential shadow voxels
  shadow_voxel <- data.frame()

  # loop over voxel layers along stem height
  for(height in seq(bottom, top, by=voxel_size)) {

    # ??? height as numeric leads to empty results for some heights ???
    plane <- subset(input_vox, Z == as.character(height), select=c(X,Y))

    # height layer with pcd voxels
    if(nrow(plane) > 0) {

      plane$pcd <- TRUE

      # estimate bounding box dimension
      rx <- range(plane$X)
      xlen <- length(seq(rx[1], rx[2], voxel_size))

      ry <- range(plane$Y)
      ylen <- length(seq(ry[1], ry[2], voxel_size))

      # create 2D voxelraster inside box
      x <- rep(seq(rx[1], rx[2], voxel_size), ylen)
      y <- rep(seq(ry[1], ry[2], voxel_size), each=xlen)

      bbox <- data.frame(X=x, Y=y)

      # add pcd voxels
      area <- merge(bbox, plane, by=c("X","Y"), all.x=TRUE)

      # filter: remove pcd voxel
      candidates <- subset(area, is.na(area$pcd)==TRUE, select=c(X,Y))

      # check neighbors if there are empty voxel in bbox
      if(nrow(candidates)>0) {

        # (re)set num as counter
        candidates$num <- 0

        # check each candidate for neighbors
        for(v in 1:nrow(candidates)) {

          pt <- candidates[v,]
          # > X
          candidates[v, "num"] <- ifelse( nrow(subset(plane, Y==pt$Y & X>pt$X)) > 0, candidates[v, "num"]+1, candidates[v, "num"])
          # < X
          candidates[v, "num"] <- ifelse( nrow(subset(plane, Y==pt$Y & X<pt$X)) > 0, candidates[v, "num"]+1, candidates[v, "num"])
          # > Y
          candidates[v, "num"] <- ifelse( nrow(subset(plane, X==pt$X & Y>pt$Y)) > 0, candidates[v, "num"]+1, candidates[v, "num"])
          # < Y
          candidates[v, "num"] <- ifelse( nrow(subset(plane, X==pt$X & Y<pt$Y)) > 0, candidates[v, "num"]+1, candidates[v, "num"])
        }

        # filter: only voxel with 4 neighbors
        cn4 <- subset(candidates, num==4, select=c(X,Y))

        # add shadow voxels to table
        if(nrow(cn4)> 0) {

          shadows <- data.frame(cn4, Z=height)
          shadows$shadow <- TRUE
          shadow_voxel <- rbind(shadow_voxel, shadows)
        }

      } # end check candidates

    } # end plane

  } # end loop along height

  output_vox <- rbind(input_vox, shadow_voxel)

  return(output_vox)


} # END-OF-FUNCTION


