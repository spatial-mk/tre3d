#' Display CPA.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Print/Display tree crown projection area.
#' @param input_cloud A 3D point cloud file, e.g. "tree1.xyz".
#'     Assumes first three columns to be x y z without header. Futher columns are ignored.
#' @return Return a PNG file
#' @export
#' @examples
#' print_CPA("tree1.xyz")

print_CPA <- function(input_cloud="*.xyz"){
  i = input_cloud
  wd <- dirname(i) #getwd()
  tree <- tools::file_path_sans_ext(i) #substr(filenames[i],1,6)
  tree_cloud <- paste(wd, "/",tree,".xyz",sep="")
  analysis <- paste(wd, "/analysis/",tree,"_analysis.txt",sep="")
  gd_cloud <- paste(wd, "/analysis/",tree,"_gd_cloud.xyz",sep="")
  dbh_cloud <- paste(wd, "/analysis/",tree,"_dbh_cloud.xyz",sep="")
  cpa <- paste(wd, "/analysis/",tree,"_cpa_as.txt",sep="")
  cpa_cen <- paste(wd, "/analysis/",tree,"_analysis.txt",sep="")
  output <- paste(wd, "/analysis/",tree,"_CHECK.png",sep="")
  if(file.exists(output)){
    file.remove(output)
    cat("\nFile exists (delete). Plot new file: ", basename(output))
  }
  else{
    cat("\nPlot file: ", basename(output))
  }

  ## Make a temporary copy of the point cloud for faster processing
  temp_copy <- paste("C:/Temp/",tree,".xyz",sep="")
  file.copy(tree_cloud, temp_copy)
  tree_cloud = temp_copy

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Read in data and set variables

  ## Read in polygon points from alpha shape in extra file
  cpa_as_poly <- read.table(cpa ,sep = " ", skip=1, col.names = c("x","y","id_tree","tudtag"))
  centroid_as <- read.table(cpa_cen, sep = " ", skip=7, nrows=1, col.names = c("name","CPA_cen_X","CPA_cen_Y","id_tree","tudtag"))
  tree_pos <- read.table(analysis , sep = " ", nrows=1, col.names = c("name","X_TLS","Y_TLS","ALTITUDE_TLS","id_plot","id_tree","tudtag"))
  dbh_pos <- read.table(analysis , sep = " ", skip=2, nrows=1, col.names = c("name","DBH","X_TLS","Y_TLS","tudtag"))
  gd <- read.table(analysis, sep = " ", nrows=1, skip=1, col.names = c("name","GD_TLS","tudtag"))
  gd_cloud_xyz <- read.table(gd_cloud, sep=" ", col.names = c("x","y","z"))
  dbh_cloud_xyz <- read.table(dbh_cloud, sep=" ", col.names = c("x","y","z"))
  tree_cloud_xyz <- read.table(tree_cloud, sep=" ", col.names = c("x","y","z"))

  circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }

  if(dbh_pos$X_TLS==-9999){
    dbh_pos$X_TLS = tree_pos$X_TLS
    dbh_pos$Y_TLS = tree_pos$Y_TLS
    dbh_pos$DBH = 0
  }

  stem <- circleFun(c(tree_pos$X_TLS,tree_pos$Y_TLS),gd$GD_TLS,npoints = 100)
  stem_dbh <- circleFun(c(dbh_pos$X_TLS,dbh_pos$Y_TLS),dbh_pos$DBH,npoints = 100)

  ## Remove temp copy of point cloud
  file.remove(temp_copy)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Create plot

  ## Reduce tree cloud to 2.5 cm for faster plotting
  tree_cloud_xyz <- VoxR::vox(tree_cloud_xyz, 0.01)[,1:3]
  tree_cloud_xyz <- dplyr::rename(tree_cloud_xyz, x=data...1., y=data...2., z=data...3.)

  ## Adjust min max to fit to plot
  tree_xmin <- round(min(tree_cloud_xyz$x),1)
  tree_xmax <- round(max(tree_cloud_xyz$x),1)
  tree_ymin <- round(min(tree_cloud_xyz$y),1)
  tree_ymax <- round(max(tree_cloud_xyz$y),1)
  tree_zmin <- round(min(tree_cloud_xyz$z),1)
  tree_zmax <- round(max(tree_cloud_xyz$z),1)
  gd_xmin <- round(min(gd_cloud_xyz$x),1)
  gd_xmax <- round(max(gd_cloud_xyz$x),1)
  gd_ymin <- round(min(gd_cloud_xyz$y),1)
  gd_ymax <- round(max(gd_cloud_xyz$y),1)

  p<-ggplot2::ggplot() +
    ## Plot CPA polygons and centroid
    ggplot2::geom_polygon(data=cpa_as_poly, ggplot2::aes(x=x,y=y), alpha=0.2,size=.25, color="black", fill="darkgreen") +
    ggplot2::coord_equal() + ggplot2::xlab("x (local coordinates)") + ggplot2::ylab("y (local coordinates)") +
    ## Plot point cloud of tree
    ggplot2::geom_point(data=tree_cloud_xyz, ggplot2::aes(x=x,y=y), shape=20, color="darkgrey", size=0.4) +
    ## Plot point cloud of ground stem
    ggplot2::geom_point(data=gd_cloud_xyz, ggplot2::aes(x=x,y=y), shape=20, color="red", size=0.4) +
    ## Plot tree positions
    ggplot2::geom_point(data=tree_pos, ggplot2::aes(x=X_TLS,y=Y_TLS), shape=16, size=3.5, color="darkgreen") +
    ## Plot displacement vector
    ggplot2::geom_segment(data=tree_pos,ggplot2::aes(x=X_TLS,y=Y_TLS,xend=centroid_as$CPA_cen_X,yend=centroid_as$CPA_cen_Y),
                          arrow=ggplot2::arrow(length=ggplot2::unit(0.2,"cm")), size=1) +
    ## Plot the fitted ground diameter
    ggplot2::geom_path(data=stem, ggplot2::aes(x=x,y=y), color="green", size=1) +
    ## Color bars, legend, title, axis, layout
    ggplot2::ggtitle(tree) + ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(lineheight=.8, face="bold")) +
    ggplot2::scale_x_continuous(breaks=seq(tree_xmin,tree_xmax,abs(tree_xmax-tree_xmin)/10)) +
    ggplot2::scale_y_continuous(breaks=seq(tree_ymin,tree_ymax,abs(tree_ymax-tree_ymin)/10))

  p2<-ggplot2::ggplot() +
    ## Plot point cloud of ground stem
    ggplot2::geom_point(data=gd_cloud_xyz, ggplot2::aes(x=x,y=y), shape=20, color="red", size=0.5) +
    ## Plot the fitted ground diameter
    ggplot2::geom_path(data=stem, ggplot2::aes(x=x,y=y), color="lightgreen", size=1.25) +
    ## Plot tree positions
    ggplot2::geom_point(data=tree_pos, ggplot2::aes(x=X_TLS,y=Y_TLS), shape=3, size=3.5, color="darkgreen") +
    ggplot2::coord_equal() +
    ## Color bars, legend, title, axis, layout
    ggplot2::ggtitle("GD and position check") + ggplot2::theme_bw()

  p3<-ggplot2::ggplot() +
    ## Plot point cloud of tree
    ggplot2::geom_point(data=tree_cloud_xyz, ggplot2::aes(x=x,y=z), shape=20, color="darkgrey", size=0.4) +
    ## Plot point cloud of ground stem
    ggplot2::geom_point(data=gd_cloud_xyz, ggplot2::aes(x=x,y=z), shape=20, color="red", size=0.4) +
    ## Plot point cloud of ground stem
    ggplot2::geom_point(data=dbh_cloud_xyz, ggplot2::aes(x=x,y=z), shape=20, color="cadetblue3", size=0.4) +
    ## Plot tree positions
    ggplot2::geom_point(data=tree_pos, ggplot2::aes(x=X_TLS,y=tree_zmin), shape=16, size=3.5, color="darkgreen") +
    ## Plot DBH positions
    ggplot2::geom_point(data=dbh_pos, ggplot2::aes(x=X_TLS,y=tree_zmin+1.3), shape=16, size=3.5, color="blue2") +
    ## Color bars, legend, title, axis, layout
    ggplot2::ggtitle("Side view - X") + ggplot2::theme_bw() + ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(breaks=seq(tree_xmin,tree_xmax,abs(tree_xmax-tree_xmin)/5)) +
    ggplot2::scale_y_continuous(breaks=seq(tree_zmin,tree_zmax,abs(tree_zmax-tree_zmin)/10))

  p4<-ggplot2::ggplot() +
    ## Plot point cloud of tree
    ggplot2::geom_point(data=tree_cloud_xyz, ggplot2::aes(x=y,y=z), shape=20, color="darkgrey", size=0.4) +
    ## Plot point cloud of ground stem
    ggplot2::geom_point(data=gd_cloud_xyz, ggplot2::aes(x=y,y=z), shape=20, color="red", size=0.4) +
    ## Plot point cloud of dbh stem
    ggplot2::geom_point(data=dbh_cloud_xyz, ggplot2::aes(x=y,y=z), shape=20, color="cadetblue3", size=0.4) +
    ## Plot tree positions
    ggplot2::geom_point(data=tree_pos, ggplot2::aes(x=Y_TLS,y=tree_zmin), shape=16, size=3.5, color="darkgreen") +
    ## Plot DBH positions
    ggplot2::geom_point(data=dbh_pos, ggplot2::aes(x=Y_TLS,y=tree_zmin+1.3), shape=16, size=3.5, color="blue2") +
    ## Color bars, legend, title, axis, layout
    ggplot2::ggtitle("Side view - Y") + ggplot2::theme_bw() + ggplot2::coord_equal() +
    ggplot2::scale_x_continuous(breaks=seq(tree_ymin,tree_ymax,abs(tree_ymax-tree_ymin)/5)) +
    ggplot2::scale_y_continuous(breaks=seq(tree_zmin,tree_zmax,abs(tree_zmax-tree_zmin)/10))

  p5<-ggplot2::ggplot() +
    ## Plot point cloud of DBH stem
    ggplot2::geom_point(data=dbh_cloud_xyz, ggplot2::aes(x=x,y=y), shape=20, color="cadetblue3", size=0.5) +
    ## Plot the fitted BDH diameter
    ggplot2::geom_path(data=stem_dbh, ggplot2::aes(x=x,y=y), color="blue2", size=1.25) +
    ## Plot tree positions
    ggplot2::geom_point(data=dbh_pos, ggplot2::aes(x=X_TLS,y=Y_TLS), shape=3, size=3.5, color="blue2") +
    ## Color bars, legend, title, axis, layout
    ggplot2::ggtitle("DBH and position check") + ggplot2::coord_equal() + ggplot2::theme_bw()

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Plot to file with specific resolution
  png(output, width=4000, height=1500, res=200, type="cairo")
  #grid.arrange(p2, p5, p3, p4, p, ncol=2, nrow=3)
  gridExtra::grid.arrange(p, p3, p4, p2, p5, ncol=5)
  dev.off()
  #grid.arrange(p2, p5, p3, p4, p, ncol=2, nrow=3)
  gridExtra::grid.arrange(p, p3, p4, p2, p5, ncol=5)
}
