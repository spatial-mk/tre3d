#' Voxel grid reduction of a point cloud.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Converts point cloud (x y z) to a voxel grid with defined voxel size.
#' @param input_cloud A file or data.frame containing a tree point cloud with x y z. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param voxel_size Size of the voxels (in meter) in the reduced (gridded) point cloud. Default 0.03
#' @param set_z_to_zero Set lowest tree point (z axis) to zero. Default to FALSE.
#' @return The reduced point cloud voxel grid as data.frame with column names x y z.
#' @export
#' @examples
#' cloud2vox("tree1.xyz", 0.03) # 0.7853982
cloud2vox <- function(input_cloud="*.xyz", voxel_size=0.03, set_z_to_zero=FALSE){

  ## Check what type of input is provided for tree_i: data.frame or file
  if(is.data.frame(input_cloud)){
    tree = input_cloud
  } else{
    ## Check if file exists. If not, give warning
    if (!file.exists(input_cloud)){
      stop(paste0("Point cloud file: ", input_cloud, " does not exist."))
    } else {
      ## Read the file, assuming column order: x y z
      options(readr.num_columns = 0) ## suppress the output of readr
      tree <- readr::read_delim(file = input_cloud, delim = " ", col_names = c("x","y","z"))
    }
  }

  ## Set lowest tree point (z axis) to zero
  if (set_z_to_zero==TRUE){
    tree[,3] <- tree[,3] - min(tree[,3])
  }

  ## Convert to voxelgrid with specified edge length
  tree_vox <- VoxR::vox(tree, res = voxel_size)

  ## Rename and use only first three columns
  tree_vox <- tree_vox[,1:3]
  colnames(tree_vox) <- c("x","y","z")

  ## Return the data.frame
  return(tree_vox)

} # END-OF-FUNCTION
