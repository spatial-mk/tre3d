#' Create a voxelgrid datatable on plot level to further analyse of the spatial canopy occupation.
#' @author Carsten Hess, last updated: 20.02.2020
#' @description Join the voxelgrid models of the DGM and all trees as one big voxelgrid datatable.
#' @param plot_id Name/ID of the plot/forest stand
#' @param path_voxelgrid_folder Plot based path to the folder which contains all relevant tree and DGM voxelgrid models.
#' @param clip_area Clip the underlying dgm to the spatial extend of the projected tree voxelgrid models
#' @return Saved datagrid table which contains all plot based voxelgrid models and data.table object with total voxel volumes for each type of occupation.
#' @export
#' @examples
#'

voxelgrid_plot_model <- function(plot_id="Dummy", path_voxelgrid_folder=NA, voxel_size=0.2, clip_area=TRUE) {

  ## voxelgrid tree models
  path_tree_files <- list.files(path=file.path(path_voxelgrid_models), pattern=".xyz$", full.names=T, ignore.case=T)
  ## voxelgrid DGM
  path_dgm_file <- list.files(path=file.path(path_voxelgrid_models), pattern="DTM.txt$", full.names=T, ignore.case=T)
  ## create destination folder inside sourde folder
  dir.create(file.path(path_voxelgrid_folder, "datagrid"), showWarnings = FALSE, recursive = T)

  ## read and rbind single trees
  trees <- data.frame()

  for(treefile in path_tree_files) {
    ## columns voxelgrid_tree_model: id,X,Y,Z,type,alpha,voxel
    treemodel <- read.table(treefile, header=T, dec =".", sep=";", stringsAsFactors=FALSE)
    trees <- rbind(trees, treemodel)
  }

  ## define bounding box in XY for a simple DGM clipping
  area_range_x <- c(min(trees$X)-2*voxel_size, max(trees$X)+2*voxel_size)
  area_range_y <- c(min(trees$Y)-2*voxel_size, max(trees$Y)+2*voxel_size)

  ## read DGM
  dgm <- read.table(path_dgm_file[1], header=T, dec =".", sep=";", stringsAsFactors=FALSE)
  dgm_area <- subset(dgm, dgm$X %between% area_range_x & dgm$Y %between% area_range_y)

  ## add and map Z_grd to each tree voxel
  ## dgm defines first layer of type "ground"
  ground <- subset(dgm_area, select=c(X, Y, Z))
  setnames(ground, "Z", "Z_grd")

  ## manipulate data.frame dummies
  dy1_trees <- merge(trees, ground, by=c("X","Y"), all.x=TRUE)
  dy1_trees$height_grd <- plyr::round_any( (dy1_trees$Z - dy1_trees$Z_grd), voxel_size)

  ## filter out tree voxel: only voxel above ground
  dy2_trees <- subset(dy1_trees, height_grd > 0, select=c(X, Y, Z, id, type, height_grd))

  ## +++ create datatable ++++++++++++++++++++++++++++++++++++++++++++++
  ## add missing fields to dgm_area
  dgm_area$id <- NA
  dgm_area$type <- "ground"
  dgm_area$height_grd <- 0

  datagrid <- rbind(dgm_area, dy2_trees)

  rdata_file <- file.path(path_voxelgrid_folder, "datagrid", paste(plot_id, "datagrid.RData", sep="_"))
  csv_file <-  file.path(path_voxelgrid_folder, "datagrid", paste(plot_id, "datagrid.csv", sep="_"))

  save(datagrid, file=file.path(datafile))
  write.table(datagrid, file=file.path(csvfile), quote = FALSE, sep=";", dec = ".", row.names = F, col.names = T)

  ## simple aggregation of voxel positions to estimate total volumes for each type of occupation (single,twofold,threefold ...)
  occupations <- data.frame(num=1:5, occ=c("single","twofold","threefold","fourfold","fivefold"))
  dy1 <- subset(datagrid, type!="ground")
  dy1$num <- 1
  ## sum(num) for each Voxel position
  dy2 <- aggregate(num ~ X + Y + Z, dy1, sum)
  ## add "factor" occupation for each voxel position
  dy3 <- merge(dy2, occupations, by=c("num"), all.x=TRUE)
  dy3$num <- 1

  agg <- aggregate(num ~ occupations, dy3, sum)
  agg$volume <- agg$num * voxel_size

  return(agg)

} # END-OF-FUNCTION
