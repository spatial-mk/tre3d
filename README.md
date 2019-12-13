# tre3d

Beta version!

13.12.2019

Matthias Kunz, matthias.kunz[at]tu-dresden.de

Requirements:
+ Install Rtools (https://cran.r-project.org/bin/windows/Rtools/)
+ Install LAStools under C: (https://rapidlasso.com/lastools/; C:/LAStools)
+ Install CloudCompare under C:/Program Files/ (https://www.danielgm.net/cc/; C:/Program Files/CloudCompare)

Installation:

> install.packages("devtools")

> devtools::install_github("spatial-mk/tre3d")

Test:

#> df = tre3d::tree1

#> tre3d::get_2d_alpha_shape(input_cloud = df, plot = T)
