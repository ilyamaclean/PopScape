#' A 50 m resolution SpatRaster of habitat suitability for a hypothetical woodland species
#'
#' A spatial dataset of habitat suitability for the Lizard Peninsula: an area bounded by 159170,
#' 182970, 10740, 30940 (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference
#' system (CRS: 27700). Different land cover types are presented by numeric values as given in
#' the inbuilt dastaset `lclass`.
#'
#' @format A PackedSpatRaster object with 404 rows and 476 columns
#' @source derived from `landcover` and `veghgt`
"habsuit"

#' A 50 m resolution SpatRaster of landcover
#'
#' A spatial dataset of land cover types for the Lizard Peninsula: an area bounded by 159170,
#' 182970, 10740, 30940 (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference
#' system (CRS: 27700). Different land cover types are presented by numeric values as given in
#' the inbuilt dastaset `lclass`.
#'
#' @format A PackedSpatRaster object with 404 rows and 476 columns
#' @source derived from LCM2021 \url{https://www.ceh.ac.uk/data/ukceh-land-cover-maps}
"landcover"

#' Data.frame giving legend for landcover
#'
#' A data,frame with the following columns
#' \describe{
#'  \item{Habitat}{Description of habitat type}
#'  \item{Class}{Numeric Value corresponding to values in inbuilt dataset `landcover`}
#'  }
"lclass"

#' Simulated species occurrences in habitat patches across the Lizard Peninsula
#'
#' A data,frame with the following columns
#' \describe{
#'   \item(x){OSGB grid reference - Easting}
#'   \item(y){OSGB grid reference - Northing}
#'   \item(pa){presence (1) or absence (0)}
#'   \item{Area}{Patch size (ha)}
#' }
"metapopdata"

#' A 50 m resolution SpatRaster of vegetation height
#'
#' A spatial dataset of vegetation height (m) for the Lizard Peninsula: an area bounded by 159170,
#' 182970, 10740, 30940 (xmin, xmax, ymin, ymax) using the Ordance Survey GB Grid Reference
#' system (CRS: 27700).
#'
#' @format A PackedSpatRaster object with 404 rows and 476 columns
#' @source derived from LiDAR data \url{http://www.tellusgb.ac.uk/}
"veghgt"

