## A function that computes the mean fractal dimension (MeanFrac)
## See: Ehbrecht et al., 2017, Agricultural and Forest Meteorology, DOI:10.1016/j.agrformet.2017.04.012
## Matthias Kunz, 15.10.2019
## Note: Angular resolution depends on the scanner settings during the scanning

## ****************************************************************************
#' Computation of mean fractal dimensions for a point cloud with definded cross section based on angular resolution
#' @author Matthias Kunz, last updated: 19.12.2019
#' @description Computation of mean fractal dimensions for a point cloud with definded cross section based on angular resolution
#' @param input_df Data frame that holds the point cloud
#' @param angular_resolution Angular resolution of the scan, i.e. 360/number of points. Default 360/2560.
#' @param plotting Logical. Should the MeanFrac polygons be plotted. Note: This may takes a long time. Default FALSE.
#' @return List with MeanFrac and fractal dimension for each cross section
#' @export
#' @importFrom foreach "%dopar%"
#' @examples
#' MeanFrac(input_cloud=df, angular_resolution=360/2560, plotting=F)
## ****************************************************************************

## Required R packages
#library(sp) # Polygon functions
#library(rgl) # 3D plotting
#library(foreach) ## R-libraries needed for paralell processing
#library(doParallel) ## R-libraries needed for paralell processing
#library(data.table) ## rbindlist function

MeanFrac <- function(input_cloud=df,
                 angular_resolution=360/2560,
                 plotting=F){

  ## Rename first three variables to x,y,z
  df <- input_cloud
  colnames(df) <- c("x","y","z")

  ## ----------------------------------------------------------------------------
  ## Convert coordinates from cartesian to polar system
  cat("\n>>> Convert from cartesian to polar coordinates")
  df$r = sqrt(df$x^2 + df$y^2 + df$z^2)
  df$inclination = acos(df$z / df$r)
  df$inclination_deg = tre3d::rad2deg(df$inclination)
  df$azimuth = atan2(df$y,df$x)
  df$azimuth_deg = ifelse(rad2deg(df$azimuth)<0, tre3d::rad2deg(df$azimuth)+360, tre3d::rad2deg(df$azimuth))

  ## Check reverse conversion
  df$x1 <- df$r * sin(df$inclination) * cos(df$azimuth)
  df$y1 <- df$r * sin(df$inclination) * sin(df$azimuth)
  df$z1 <- df$r * cos(df$inclination)
  if (identical(round(df$x,4), round(df$x1,4)) &
      identical(round(df$y,4), round(df$y1,4)) &
      identical(round(df$z,4), round(df$z1,4))){
    cat("\nConversion successfull (x,y,z to r,inclination,azimuth)\n")
  }
  df$x1 <- NULL
  df$y1 <- NULL
  df$z1 <- NULL

  ## ----------------------------------------------------------------------------
  ## Extract cross-sections and compute MeanFrac
  cat(crayon::bold(crayon::blue("\n>>> Extract crosssections and FRAC\n")))
  n_section = 360 / angular_resolution

  #n_section <- 2560 # number of cross-sections, i.e. 360/2560=0.140625 deg
  sections <- data.frame(sections = seq(from = 0, to = 180-(360/n_section), length.out = n_section/2) ,
                         col = colorRampPalette(c("black", "blue3", "blue","cyan", "darkseagreen" , "green", "yellow", "red", "darkorchid"))(n_section/2),
                         FRAC=NA, area=NA, perimeter=NA, section_id=seq(1:n_section/2))

  ## Create key on a data.table for faster query
  dt <- data.table::as.data.table(df)
  data.table::setkey(dt, azimuth_deg)

  ## Register a parallel backend and the number of cores to use
  cat(format(Sys.time(),usetz = TRUE),"\n")
  nc = parallel::detectCores() # Number of cores
  cl = parallel::makeCluster(nc) # Number of clusters
  doParallel::registerDoParallel(cl) # Register clusters

  ## PARELLEL PROCESSING of sections
  result<-foreach::foreach(i = sections$sections) %dopar% {
    ## Get crosssection
    #cs <- dplyr::filter(df, azimuth_deg >= i & azimuth_deg < i+(360/n_section) | azimuth_deg >= i+180 & azimuth_deg < i+(360/n_section)+180)
    cs <- dt[dt$azimuth_deg >= i & dt$azimuth_deg < i+(360/n_section) | dt$azimuth_deg >= i+180 & dt$azimuth_deg < i+(360/n_section)+180, ]
    cs <- as.data.frame(cs)

    cs$side <- ifelse(cs$azimuth_deg > i & cs$azimuth_deg < i+(360/n_section), 1, 2)

    ## Make angles on one side negative, as inclination has 0 degree in zenith
    cs$inclination_deg[cs$side==1] <- cs$inclination_deg[cs$side==1] * -1

    ## Rotate inclination degrees by 90 degrees to make 0 deg the horizontal (so zenith should be 90 deg now)
    cs$inclination_deg = cs$inclination_deg + 90

    ## Convert to cartesian coordinates to get area and perimeter of crosssection polygon
    #deg2rad <- function(angle=0){angle = angle * pi / 180; return(angle)}
    #rad2deg <- function(angle=0){angle = angle * 180 / pi; return(angle)}
    cs$xx <- cs$r * cos(tre3d::deg2rad(cs$inclination_deg))
    cs$yy <- cs$r * sin(tre3d::deg2rad(cs$inclination_deg))

    ## Generate polygon from points ordered by inclination angle
    cs <- cs[order(cs$inclination_deg),]

    ## Get the sorted coordinates, start and end point is the scanner center at (0,0)
    cs_poly <- data.frame(x=c(0, cs$xx, 0), y=c(0, cs$z, 0))

    ## Convert to polygon
    p = sp::Polygon(cs_poly)
    ps = sp::Polygons(list(p),1)
    sps = sp::SpatialPolygons(list(ps))

    ## Get area and perimeter of the crosssection polygon
    area <- sps@polygons[[1]]@Polygons[[1]]@area
    perimeter <- spatialEco::polyPerimeter(sps)

    ## Get centroid
    centroid <- sps@polygons[[1]]@Polygons[[1]]@labpt

    ## Compute fractal dimension index, see p.3 in Ehbrecht et al. 2017 AFM
    FRAC = (2*log(0.25*perimeter)) / log(area)

    if (plotting == T){
      ## Plot the crosssection polygon in 2D
      png(paste0("section_",sections[sections$sections==i,]$section_id,".png"), res = 100, width=1000, height=1000, type="cairo")
      plot(cs_poly, pch='.', axes=T, ylab="", xlab="", asp=1)
      polygon(cs_poly, col=rgb(1, 0, 0, 0.5)) # Plot the polygon
      title(paste0("Area [m^2]: ", round(area,2), "\nPerimeter [m]: ", round(perimeter,2), "\nFRAC:",round(FRAC,2)))
      dev.off()

      ## Plot in 3D
      #plot3d(cs[,1:3], add=T, col=sections[sections$sections==i,]$col,
      #       xlim=c(-20,20), ylim=c(-20,20), zlim=c(-20,20),
      #       aspect = "iso", axes=T)
    }

    result_section <- data.frame(area=area,perimeter=perimeter, FRAC=FRAC)
    #sections[sections$sections==i,]$area = area
    #sections[sections$sections==i,]$perimeter = perimeter

  } # END-OF-FOREACH-LOOP

  section_result <- as.data.frame(data.table::rbindlist(result))
  sections$area <- section_result$area
  sections$perimeter <- section_result$perimeter
  sections$FRAC <- section_result$FRAC

  parallel::stopCluster(cl) ## Stop paralell processing
  foreach::registerDoSEQ() ## Make processing sequential again
  cat(format(Sys.time(),usetz = TRUE))

  rm(df)
  rm(dt)
  gc() # Garbage collection

  ## Compute MeanFrac
  MeanFrac = mean(sections$FRAC, na.rm = T)
  result_MeanFrac <- data.frame(MeanFrac=MeanFrac)
  return(list(result_MeanFrac=MeanFrac, Crosssection=sections))
}
