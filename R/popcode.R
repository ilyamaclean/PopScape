#' @title Simulates spatially explicit population model
#' @description Simulates population dynamics over a landscape, including metapopulation
#' dynamics
#' @param popden multi-layer SpaRaster of maximum potential population densities in each timestep
#' @param birthrate multi-layer SpaRaster of number of births as a function of environmental
#' suitability within each timestep
#' @param survival multi-layer SpaRaster of number of survival rates as a function of environmental
#' suitability within each timestep
#' @param avdispdist average dispersal distance capability (km)
#' @param fracdisp fraction of population emigrating from each population annually
#' @param em_surv fraction of emmigrants that typically survive
#' @param Nt optionally initial population size. Set at carrying capacity if unknown
#' stochastically (TRUE) or determenistically (FALSE)
#' @return a SpatRaster of predicting population density within each time-step
#' @rdname PopDynamicsSim
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib PopScape, .registration = TRUE
#' @export
#' @description. See vignette for further details
PopDynamicsSim <- function(popden, birthrate, survival,  avdispdist, fracdisp, emmsurv, Nt = NA) {
  # Delineate habitat patches
  hab <- popden[[1]] * 0 + 1
  patch<- patches(hab, directions = 8)
  patchm <- as.matrix(patch, wide = TRUE)
  patchm <- renumberPatchIds(patchm)
  patchm[patchm == 0]<-NA
  # Calculate patch area
  Ai <- (calculatePatchSizes(patchm) * res(patch)[1]^2) / (100 * 100) # convert to Ha
  # Calculate matrix of birth rates
  n <- dim(popden)[3]
  birthm <- arraytomat(as.array(birthrate), patchm)
  survm <-  arraytomat(as.array(survival), patchm)
  Km <- arraytomat(as.array(popden), patchm)
  for (i in 1:n) Km[,i] <- Km[,i] * Ai
  Km <- round(Km, 0)
  Km[Km == 0] <- 1
  # Calculate Matrix of patch distances weighted by average suitability
  apopden <- apply3D(as.array(popden))
  dij <- patchdistcpp(patchm, apopden) * (res(popden)[1] / 1000) # convert to Km
  # Set initial population size
  if (class(Nt) == "logical") {
    Nt <- Km[,1]
  } else Nt <- calculateAveragePatchQuality(patchm, as.matrix(Nt, wide = TRUE))
  # Run population model in c++
  N <- PopSim(Nt, Km, birthm, survm, dij, avdispdist, fracdisp, emmsurv, Ai)
  # Produce raster output
  Nr <- .rast(mattoarray(N, patchm), patch)
  Kr <- .rast(mattoarray(Km, patchm), patch)
  Kf <- Nr / Kr
  popdeno <- Kf * popden
  return(popdeno)
}
