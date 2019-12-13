## Script that computes the alpha shape, area, and centroid of a 2D point set
## Matthias Kunz, 12.11.2018

#' Compute area, perimeter and centroid of an 2d alpha-shape for a set of points.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Compute area, perimeter and centroid of an 2d alpha-shape.\cr
#'     Alpha-shape is checked for connectivity and circular graphs.\cr
#'     If problems occur alpha-value is increased by 0.2. To avoid collinearity a small random jitter is added to the points before alpha-shape computation.
#' @param input_cloud A file or data.frame containing x and y points. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param alpha Alpha value for alpha-shape computation. Default 0.5
#' @param plot Should the alpha shape be plotted. Default to FALSE.
#' @return A data.frame containing: area, perimeter, cpa_cen_x, cpa_cen_y.
#' @export
#' @examples
#' get_2d_alpha_shape(data.frame(x=runif(20),y=runif(20)), alpha=1.0)

## Actual function
get_2d_alpha_shape <- function(input_cloud=data.frame(x=runif(20),y=runif(20)), alpha=1.0, plot=FALSE){

  ## Check what type of input is provided for tree_i: data.frame or file
  if(is.data.frame(input_cloud)){
    df = input_cloud
    #df <- df[,1:2] # reduce to x and y coordinates only
  } else{
    ## Check if file exists. If not, give warning
    if (!file.exists(input_cloud)){
      stop(paste0("Point cloud file: ", input_cloud, " does not exist."))
    } else {
      ## Read the file, assuming column order: x y z
      options(readr.num_columns = 0) ## suppress the output of readr
      df <- readr::read_delim(file = input_cloud, delim = " ", col_names = c("x","y","z"))
      #df <- df[,1:2] # reduce to x and y coordinates only
    }
  }

  ##---
  ## Use LAStools to reduce the number of points, i.e. only use highest point in a grid
  ## Store the data to a file
  write.table(df, "df_temp.txt", col.names = F, sep = " ", row.names = F)
  df_temp<-paste0(getwd(),"/df_temp.txt")
  grid_temp<-paste0(getwd(),"/grid_temp.txt")
  ## Change thin step based on number of points
  if (nrow(df)>1000000){step_size=0.02} # >1 Mio. --> 2 cm
  else if (nrow(df)>5000000){step_size=0.03} # >5 Mio. --> 3 cm
  else if (nrow(df)>10000000){step_size=0.04} # >10 Mio. --> 4 cm
  else if (nrow(df)>25000000){step_size=0.05} # >25 Mio. --> 5 cm
  else {step_size=0.01} # Default to 1 cm
  system("cmd.exe", input=paste0("C:\\LAStools\\bin\\lasthin.exe -i ", df_temp,
                                 " -iparse xy -oparse xy -step ", step_size," -central -o ", grid_temp))
  df<-readr::read_delim(file = grid_temp, delim = " ", col_names = c("x","y"))
  unlink(c(df_temp,grid_temp))
  ##---

  ## Compute alpha-shape
  ## add random jitter to points
  df[,1] = df[,1] + rnorm(nrow(df),mean = 0.001, sd=0.0005)
  df[,2] = df[,2] + rnorm(nrow(df),mean = 0.001, sd=0.0005)

  ## Compute alpha-shape
  as <- alphahull::ashape(df, alpha = alpha)

  ##https://rpubs.com/geospacedman/alphasimple
  asg = igraph::graph.edgelist(cbind(as.character(as$edges[, "ind1"]),
                             as.character(as$edges[, "ind2"])), directed = FALSE)

  ## Check if graph is correct, if not recompute alpha-shape and give warning
  if (!igraph::is.connected(asg)) {
    as <- alphahull::ashape(df, alpha = alpha+0.2)
    asg = igraph::graph.edgelist(cbind(as.character(as$edges[, "ind1"]),
                               as.character(as$edges[, "ind2"])), directed = FALSE)
    warning("Graph not connected")
  }
  if (any(igraph::degree(asg) != 2)) {
    as <- alphahull::ashape(df, alpha = alpha+0.2)
    asg = igraph::graph.edgelist(cbind(as.character(as$edges[, "ind1"]),
                               as.character(as$edges[, "ind2"])), directed = FALSE)
    warning("Graph not circular")
  }
  if (igraph::clusters(asg)$no > 1) {
    as <- alphahull::ashape(df, alpha = alpha+0.2)
    asg = igraph::graph.edgelist(cbind(as.character(as$edges[, "ind1"]),
                               as.character(as$edges[, "ind2"])), directed = FALSE)
    warning("Graph composed of more than one circle")
  }

  cutg = asg - igraph::E(asg)[1]
  # find chain end points
  ends = names(which(igraph::degree(cutg) == 1))
  path = igraph::get.shortest.paths(cutg, ends[1], ends[2])[[1]]
  # this is an index into the points
  pathX = as.numeric(names(unlist(path)))
  # join the ends
  pathX = c(pathX, pathX[1])

  ## Convert to polygon
  p = sp::Polygon(df[pathX,c("x","y")])
  ps = sp::Polygons(list(p),1)
  sps = sp::SpatialPolygons(list(ps))
  #plot(sps) # Plot the polygon

  ## Get polygon coordinates
  as_poly <- df[pathX,c("x","y")]

  ## Get area
  area <- sps@polygons[[1]]@Polygons[[1]]@area

  ## Get perimeter of polygon
  perimeter <- spatialEco::polyPerimeter(sps)

  ## Get centroid
  centroid <- sps@polygons[[1]]@Polygons[[1]]@labpt

  ## Plotting
  if (plot==T){
    plot(df, pch=20, cex=0.3, asp=1)
    polygon(as_poly, col = rgb(0, 1, 0, 0.75))
    lines(as_poly, lwd=2)
    points(centroid[1], centroid[2], pch=7, cex=1)
    title(paste0("Area [m^2]: ", round(area,2), "\nPerimeter [m]: ", round(perimeter,2),"\nAlpha: ",alpha))
  }

  ## Convert polygon to simple data.frame
  hull_poly <- ggplot2::fortify(sps)
  hull_poly <- dplyr::rename(hull_poly, x=long, y=lat)
  hull_poly <- dplyr::select(hull_poly, x, y)

  ## Return result
  result <- data.frame(area=area, perimeter=perimeter, cpa_cen_x=centroid[1], cpa_cen_y=centroid[2])
  result <- list("area"=result$area, "perimeter"=result$perimeter,
                "cpa_cen_x"=result$cpa_cen_x,  "cpa_cen_y"=result$cpa_cen_y,
                "polygon"=hull_poly)
  return(result)
} ## END-OF-FUNCTION
## ****************************************************************************
