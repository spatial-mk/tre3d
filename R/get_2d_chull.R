## Function that computes the convex hull of a 2d set of points
## Return the area, perimeter and centroid coordinates
## Matthias Kunz
## 26.03.2019

#' Compute area, perimeter and centroid of an 2d convex hull for a set of points.
#' @author Matthias Kunz, last updated: 27.04.2019
#' @description Compute area, perimeter and centroid of an 2d convex hull.
#' @param input_cloud A file or data.frame containing x and y points. If file, it assumes first three columns to be x y z without header. Futher columns are ignored.
#' @param plot Should the convex hull be plotted. Default to FALSE.
#' @return A data.frame containing: area, perimeter, cpa_cen_x, cpa_cen_y.
#' @examples
#' get_2d_chull(data.frame(x=runif(20),y=runif(20)))

## Actual function
get_2d_chull<-function(input_cloud=data.frame(x=runif(20),y=runif(20)), plot=FALSE){

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

  ## Compute convex hull
  ## add random jitter to points
  df <- as.data.frame(df)
  df[,1] = df[,1] + rnorm(nrow(df),mean = 0.001, sd=0.0005)
  df[,2] = df[,2] + rnorm(nrow(df),mean = 0.001, sd=0.0005)

  hull<-chull(df[,1],df[,2])

  # join the ends
  hull = c(hull, hull[1])

  ## Convert to polygon
  p = sp::Polygon(df[hull,c("x","y")])
  ps = sp::Polygons(list(p),1)
  sps = sp::SpatialPolygons(list(ps))
  #plot(sps) # Plot the polygon

  ## Get polygon coordinates
  as_poly <- df[hull,c("x","y")]

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
    points(centroid[1], centroid[2], pch=7, cex=2)
    title(paste0("Area [m^2]: ", round(area,2), "\nPerimeter [m]: ", round(perimeter,2)))
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
}
