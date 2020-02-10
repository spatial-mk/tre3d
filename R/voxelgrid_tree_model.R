#' Create a voxelgrid tree model including the voxel attributs tree height, segment (crown, stem) and type of occupation (pcd, ashape3d)
#' @author Carsten Hess, last updated: 10.02.2020
#' @description Merge voxelmodels of tree point cloud and crown 3D alphashape.
#' @param input_ashape3d_voxelmodel Voxelmodel of the tree crown 3D alphashape (see get_3d_alpha_shape)
#' @param input_pcd_voxelmodel Voxelmodel of the tree point cloud (Size of the voxels (see pcd2vox)
#' @param alpha_value Required value of the pre-processed tree crown 3D alphashape (ashaped3d).
#' @param voxel_size Required value of the pre-prcessed voxelmodels (pcd, ashape).
#' @param cbh Required value for input_pcd_voxel segmentation into crown and stem.
#' @return Voxelgrid tree model
#' @export
#' @examples
#'
