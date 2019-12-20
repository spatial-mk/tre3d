## R Skript that visualise QSM data from Pasi Raumonen
## Matthias Kunz, 03.07.2017

#' Visualize 3D QSM model
#' @author Matthias Kunz, last updated: 29.04.2019
#' @description Visualize 3D QSM model from TreeQSM 2.30
#' @references \url{https://github.com/InverseTampere/TreeQSM}
#' @param qsm qsm A file containing a QSM, e.g. cyl_data_*.txt
#' @param n_seg Number of segments that draw the cylinder. Default 12
#' @param axis Display the cylinder axis. Default to FALSE.
#' @param cylinders Display the cylinder. Defaut to TRUE.
#' @param circles Display circles at the end of the cylinder. Default to FALSE.
#' @param transparent Makes the cylinder transparent. Default to TRUE with transparency of 0.3
#' @return Does not return anything at the moment.
#' @export
#' @examples
#' display_qsm("cyl_data*.txt")

display_qsm <- function(qsm="cyl_data*.txt", n_seg=12, axis=FALSE, cylinders=TRUE, circles=FALSE, transparent=TRUE){

  ## ----------------------------------
  ## Read in QSM and check if file exists. If not, give warning
  if (!file.exists(qsm)){
    stop(paste0("QSM file: ", qsm, " does not exist."))
  } else {
    ## Read the file, assuming column order: x y z
    options(readr.num_columns = 0) ## suppress the output of readr
    df <- read.table(file = qsm, sep= ",", header = T)
  }

  ## Rename columns for better understanding
  df <- dplyr::rename(df, r=radius, l=length,
                      x=start_1, y=start_2, z=start_3,
                      xdir=axis_1, ydir=axis_2, zdir=axis_3)

  ## ----------------------------------
  ## Print the QSM
  while (rgl::rgl.cur() > 0) { rgl::rgl.close() }
  rgl::rgl.open()# Open a new RGL device
  rgl::rgl.bg(color = "white") # Setup the background color
  rgl::par3d(windowRect=c(50,50,1200,800), mouseMode = "trackball") # Change size of rgl window
  for (i in 1:nrow(df)){
    tre3d::cylinder3d(p_start = c(df[i,c("x")], df[i,c("y")], df[i,c("z")]), # Start point of the cylinder
               orient = c(df[i,c("xdir")], df[i,c("ydir")], df[i,c("zdir")]), # Directional vector of the cylinder axis
               radius = df[i,c("r")],#/1000,      # Radius of the cylinder
               length = df[i,c("l")],#/1000,         # Length of the cylinder
               order = df[i,c("BranchOrder")],          # Branch order level of the cylinder
               n_seg = n_seg,         # Number of segments that draw the cylinder
               axis = axis,           # Display the cylinder axis
               cylinders = cylinders,      # Display the cylinder
               circles = circles,        # Display circles at the end of the cylinder
               transparent = transparent)
  }


} ## END-OF-FUNCTIONS
## ****************************************************************************
