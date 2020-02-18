#' Create a voxelgrid tree model including the voxel attributs tree height, segment (crown, stem) and type of occupation (pcd, ashape3d)
#' @author Carsten Hess, last updated: 10.02.2020
#' @description Merge voxelmodels of tree point cloud and crown 3D alphashape.
#' @param pcd_vox Voxelmodel of the tree point cloud (Size of the voxels (see pcd2vox)
#' @param as3d_vox Voxelmodel of the tree crown 3D alphashape (see get_3d_alpha_shape)
#' @param alpha_value Documented modelparameter of the pre-processed tree crown 3D alphashape (ashaped3d).
#' @param voxel_size Documented modelparameter of the pre-prcessed voxelmodels (pcd, ashape).
#' @return List object of tree with underlying modelparameter and Voxelgrid tree model (X,Y,Z,type)
#' @export
#' @examples
#'

voxelgrid_tree_model <- function(tree_id="Dummy", pcd_vox=NA, as3d_vox=NA, alpha_value=0.5, voxel_size=0.2) {

  ## point cloud voxelmidel
  if(is.data.frame(pcd_vox) == TRUE) {

    ## precaution check > fix roundings in case they were missing
    pcd_vox$X <- plyr::round_any(pcd_vox$x, voxel_size)
    pcd_vox$Y <- plyr::round_any(pcd_vox$y, voxel_size)
    pcd_vox$Z <- plyr::round_any(pcd_vox$z, voxel_size)
    ## precaution check > only unique voxel
    pcd <- unique(subset(pcd_vox, select=c(X,Y,Z)))
    ## generate unique voxel ID based on coordinates
    pcd$id <- paste(as.character(pcd$X*100), as.character(pcd$Y*100), as.character(pcd$Z*100), sep="_")

  } else {
    stop("Missing general point cloud voxelmodel")
  }

  ## ashape3d voxelmodel
  if(is.data.frame(as3d_vox) == TRUE) {

    ## precaution check > fix roundings in case they were missing
    as3d_vox$X <- plyr::round_any(as3d_vox$X, voxel_size)
    as3d_vox$Y <- plyr::round_any(as3d_vox$Y, voxel_size)
    as3d_vox$Z <- plyr::round_any(as3d_vox$Z, voxel_size)

    ## precaution check > only unique voxel
    as3d <- unique(as3d, select=c(X,Y,Z))

    ### generate unique voxel ID
    as3d$id <- paste(as.character(as3d$X*100), as.character(as3d$Y*100), as.character(as3d$Z*100), sep="_")

  }

  ## remove voxel in ashape3d that are already derived by pcd voxelmodel
  if(nrow(as3d) > 0) {
    as3d_only <- subset(as3d, !(id %in% pcd$id))
    as3d <- subset(as3d_only, select=c(X,Y,Z))
  } else {
    as3d <- data.frame(X=c(),Y=c(),Z=c())
  }

  ## modify and combine both parts
  pcd_vox <- subset(pcd, select=c(X,Y,Z))
  pcd_vox$type <- "pcd"

  as3d_vox <- subset(as3d, select=c(X,Y,Z))
  as3d_vox$type <- "ashape"

  tree_vox <- rbind(pcd_vox,as3d_vox) # voxelgrid
  ## Alternative output as list object
  #treedata <- list(id=tree_id, alpha=alpha_value, voxel=voxel_size, voxelgrid=tree_vox)
  tree_vox$id <- tree_id
  tree_vox$alpha <- alpha_value
  tree_vox$voxel <- voxel_size
  ## output optimised for direct generation of voxelgrid_plot_datatable
  treedata <- subset(tree_vox, select=c(id,X,Y,Z,type,alpha,voxel))

  return(treedata)

} # END-OF-FUNCTION
