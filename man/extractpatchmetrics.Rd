% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metapopcode.R
\name{extractpatchmetrics}
\alias{extractpatchmetrics}
\title{Extracts patch metrics to data.frame of occurrence data}
\usage{
extractpatchmetrics(
  occdata,
  landcover,
  conpredictors,
  habsuit = NA,
  Oj = NA,
  thresh = 0.5,
  alpha = 10,
  aos = 0.1,
  minden = 0.01,
  maxden = 20
)
}
\arguments{
\item{occdata}{a data.frame of occurrence data with x & y coordinates as returned
by [subsample()].}

\item{landcover}{a SpatRaster landcover classes.}

\item{conpredictors}{a single or multi-layer SpatRaster of other continuous predictors}

\item{habsuit}{optionally - a SpatRaster of predicted occupancy probability. Set to
NA on first iteration.}

\item{Oj}{a vector of 0s and 1s indicating which patches are occupied - set to NA
on first iteration.}

\item{thresh}{threshold of suitability below which patch density is assumed to equal minden}

\item{alpha}{inverse of average dispersal distance (km) - used to calculate connectivity}

\item{aos}{parameter scaling logit `habsuit` to population density}

\item{minden}{minimum potential density in suitable habitat (individuals per hectare)}

\item{maxden}{maximum density (individuals per hectare)}
}
\value{
a data.frame matching `occdata` but with the following columns added:
\describe{
 \item{habclass}{land cover type extracted from landcover}
 \item{alpha}{inverse of average dispersal distance (same units as x and y)}
 \item{Area}{Patch area in ha}
 \item{Connectivity}{Connectivity to neighboring occupied patches}
 \item{names(conpredictors)}{Additional data extracted from conpredictors}
}
}
\description{
Extracts habitat patch size, connectivity, landcover type and other relevent
environmental data from provided SpatRasters. Called repeatedly by [fitmetapopsdm()]
}
