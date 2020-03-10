#' tre3d version information
#' @author Matthias Kunz, last updated: 09.03.2020
#' @description Gives information of the tre3d version and changelog
#' @export
#' @examples
#' tre3d_version()
tre3d_version <- function(){
  cat(crayon::bold(crayon::blue("\ntre3d version 0.9 'Pretty Pinguin', last updated 09.03.2020\n")))
}
