#' Voxelization of an 3D alphashape object
#' @author Carsten Hess, last updated: 10.02.2020
#' @description Voxelization of an 3D alphashape crown model into a voxelgrid model
#' @param input_ashape3d Tree crowns ashape3d object based on alphashape3d R package (see get_3d_alpha_shape)
#' @param voxel_size Size of the voxels (in meter) of the underlying voxel grid. Default 0.2
#' @return Simple voxelgrid point cloud (X, Y, Z) of all voxels inside the given 3D alphashape.
#' @export
#' @examples
#'

ashape2vox <- function(input_ashape3d="", voxel_size=0.2) {

  # first simple check that an ashape3d object was delivered
  if(length(input_ashape3d)>0) {
    as3d <- input_ashape3d
  } else {
    stop('Missing ashape3d object! \n\n')
  }

  # get underlying point cloud of the ahsape3d
  pcd <- as.data.frame(as3d$x)

  # Set colnames to x y z
  colnames(pcd) <- c("x","y","z")

  # create closed bounding box voxelraster including the 3d alpha shape
  # estimate range dimensions of this bounding box
  xmin <- plyr::round_any(min(pcd$x), voxel_size, f=floor) - (2*voxel_size)
  xmax <- plyr::round_any(max(pcd$x), voxel_size, f=ceiling) + (2*voxel_size)
  xlen <- length(seq(xmin, xmax, voxel_size))

  ymin <- plyr::round_any(min(pcd$y), voxel_size, f=floor) - (2*voxel_size)
  ymax <- plyr::round_any(max(pcd$y), voxel_size, f=ceiling) + (2*voxel_size)
  ylen <- length(seq(ymin, ymax, voxel_size))

  zmin <- plyr::round_any(min(pcd$z), voxel_size, f=floor) - (2*voxel_size)
  zmax <- plyr::round_any(max(pcd$z), voxel_size, f=ceiling) + (2*voxel_size)
  zlen <- length(seq(zmin, zmax, voxel_size))

  # create data frame bounding box voxelraster as data frame
  x <- rep(seq(xmin, xmax, voxel_size), ylen*zlen)
  y <- rep(rep(seq(ymin, ymax, voxel_size), each=xlen), zlen)
  z <- rep(seq(zmin, zmax, voxel_size),each=xlen*ylen)

  # voxelraster of bounding box
  bbox <- data.frame(X=x, Y=y, Z=z)

  bbox.mat <- as.matrix(bbox)

  tryCatch({
    in3d <- alphashape3d::inashape3d(as3d, indexAlpha=1, points = bbox.mat)
  }, warning = function(war) {
    in3d <- 0
  }, error = function(err) {
    in3d <- -1
  }, finally = {
    in3d <- alphashape3d::inashape3d(as3d, indexAlpha=1, points = bbox.mat)
  })

  # Check IN/OUT TRUE/FALSE blue/red
  inoutcolors <- ifelse(in3d, "blue", "red")

  # Add inoutcolors to points
  bbox$cols <- inoutcolors

  # Subset of IN (blue) points
  bbox_in <- subset(bbox, cols=="blue", select=c(X,Y,Z))

  # fix possible rounding issues => to many digits
  bbox_in$X <- plyr::round_any(bbox_in$X, voxel_size)
  bbox_in$Y <- plyr::round_any(bbox_in$Y, voxel_size)
  bbox_in$Z <- plyr::round_any(bbox_in$Z, voxel_size)

  ## Voxelraster of a-shaped crown
  as3d_voxel <- unique(bbox_in)

  return(as3d_voxel)

} # END-OF-FUNCTION

