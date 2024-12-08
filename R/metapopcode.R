#' @title Parameterises a metapopulation model
#' @description Derives parameters for standard IFM metapopulation model searching
#' for best alpha
#' @param metapopdata metapopdata data,frame with columns x & y (coordinates of
#' patches), A (area) & pa (presence absence as 1 & 0)
#' @param alphamin minimum possible value of alpha
#' @param alphamax maximum possible value of alpha
#' @param n number of alphas to try between alphamin & alphamax  (ignored if alphamin = alphamax)
#' @return a list of the following:
#' \describe{
#'  \item{mu}{extinction probability of a patch of unit size (same units as A)}
#'  \item{x}{scaling of extinction risk with patch area}
#'  \item{ygamma}{parameter controlling how fast the colonisation probability approaches unity with increasing }
#'  \item{alpha}{inverse of average dispersal distance (same units as x and y)}
#' }
#' @rdname fitmetapop
#' @export
#' @description. See vignette for further details of fitting a metapopulation model
#' and the meaning of parameters. If alpha is known, set `alphamin` to be the same
#' value as `alphamax`
fitmetapop <- function(metapopdata, alphamin = 0.1, alphamax = 10, n = 100) {
  d <- dist(cbind(metapopdata$x, metapopdata$y))/1000  # Calculates matrix of distances
  if (alphamin == alphamax) {
    edis <- as.matrix(exp(-alphamin * d))
    diag(edis) <- 0
    edis <- sweep(edis, 2, metapopdata$Area, "*")
    S <- rowSums(edis[, metapopdata$pa > 0])
    mod <- suppressWarnings(glm(pa ~ offset(2 * log(S)) + log(Area), family = binomial, data = metapopdata))
    alphas <- alphamin
    i <- 1
  } else {
    aa<-seq(log(alphamin),log(alphamax),length.out = n)
    alphas<-exp(aa)
    aic<-0
    i<-1
    for (alpha in alphas) {
      edis <- as.matrix(exp(-alpha * d))
      diag(edis) <- 0
      edis <- sweep(edis, 2, metapopdata$Area, "*")
      S <- rowSums(edis[, metapopdata$pa > 0])
      mod <- suppressWarnings(summary(glm(pa ~ offset(2 * log(S)) + log(Area), family = binomial, data = metapopdata)))
      aic[i]<-mod$aic
      i<-i+1
    }
    i<-which.min(aic)[1]
    edis <- as.matrix(exp(-alphas[i] * d))
    diag(edis) <- 0
    edis <- sweep(edis, 2, metapopdata$Area, "*")
    S <- rowSums(edis[, metapopdata$pa > 0])
    mod <- suppressWarnings(glm(pa ~ offset(2 * log(S)) + log(Area), family = binomial,data = metapopdata))
  }
  be <- coef(mod)
  xhat <- be[2]
  A0 <- min(metapopdata$A[metapopdata$p > 0])
  muga <- exp(-be[1])
  mu <- A0^xhat
  ga <- muga/mu
  # save parameters
  params<-list(mu=as.numeric(mu),x=as.numeric(xhat),ygamma=as.numeric(ga),alpha=alphas[i])
  return(params)
}
#' @title  Calculates patch area
#' @description Calculates patch area in hectares
#' @param habsuit a SpatRaster of habitat suitability with unsuitable areas set to NA
#' @param asraster optional logical indicating whether to return area as a SpatRaster (TRUE)
#' or as a vector of areas.
#' @return A SpatRaster or vector of areas (Ha) of each habitat patch
#' @details the x and y units of `habsuit` must be in metres.
#' @rdname calcarea
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib PopScape, .registration = TRUE
#' @export
calcarea <- function(habsuit, thresh = 0, asraster = TRUE) {
  hab <- habsuit * 0 + 1
  hab[habsuit < thresh] <- 0
  patch<- patches(hab, directions = 8)
  patchm <- as.matrix(patch, wide = TRUE)
  patchm <- renumberPatchIds(patchm)
  patchm[patchm==0]<-NA
  A <- (calculatePatchSizes(patchm) * res(patch)[1]^2) / (100*100) # convert to Ha
  if (asraster) A <- .rast(assignPatchValues(patchm, A), patch)
  return(A)
}
#' @title  Calculates suitability-weighted inter-patch distances
#' @description Calculates a suitability-weighted measure of connectivity
#' of each habitat patch to other occupied habitat patches.
#' @param habsuit a SpatRaster of habitat suitability with unsuitable areas set to NA
#' @param thresh theshold suitability value below which patches are assumed unoccupied
#' @param alpha inverse of average dispersal distance (km)
#' @param aos parameter scaling `habsuit` logistically transformed to population density (individuals per ha)
#' @param minden minimum possible density value (individuals per Ha)
#' @param maxden maximum possible density value (individuals per Ha)
#' @param Oj optional vector indicating which patches are occupied (1) or unoccopied (0)
#' asraster optional logical indicating whether to return area as a SpatRaster (TRUE)
#' @param asraster optional logical indicating whether to return connectivity as a SpatRaster (TRUE)
#' or as a vector of areas.
#' @return A SpatRaster or vector of connectvity of each habitat patch to other occupied patches
#' @rdname calcconectivity
#' @details To calculate connectivity, first the suitability-weighted distance of each
#' pixel in patch i to each pixel in patch j is calculated and averaged. A matrix of
#  suitability-weighted averaged distances from each patch to every other is calculated.
#' The parameter alpha is then used to define a dispersal kernel representing the likelihood of
#' individuals from patch j reaching each focal patch i. For a given focal patch, this value
#' is then summed and weighted by patch carrying capacity to give a value that scales
#' with the expected total number of immigrants o any given focal patch
#' @importFrom Rcpp sourceCpp
#' @useDynLib PopScape, .registration = TRUE
#' @import terra
#' @export
calcconectivity <- function(habsuit, thresh = 0.5, alpha = 10, aos = 0.1, minden = 0.01,
                            maxden = 20, Oj = NA, asraster = TRUE) {
  # Calculate patches
  hab <- habsuit * 0 + 1
  hab[habsuit < thresh] <- 0
  patch<- patches(hab, directions = 8)
  patchm <- as.matrix(patch, wide = TRUE)
  patchm <- renumberPatchIds(patchm)
  patchm[patchm==0]<-NA
  habsuitm<-as.matrix(habsuit, wide = TRUE)
  # Calculate density with assumed abundance-occupancy relationship
  mx <- max(as.vector(habsuit), na.rm = TRUE)
  mn <- min(as.vector(habsuit), na.rm = TRUE)
  rge <- mx - mn
  if (rge > 1) habsuit <- (habsuit - mn) / rge
  lhabsuit <- suppressWarnings(log(habsuitm / (1- habsuitm)) + (thresh - 0.5))
  densty <- lhabsuit^aos
  densty[is.na(densty)] <- minden
  densty[densty < minden] <- minden
  densty[densty > maxden] <- maxden
  # Create matrix of density-weighted distances in km
  dij<-patchdistcpp(patchm, densty)*res(habsuit)[1] / 1000 # convert to km
  densty[habsuitm < thresh] <- 0
  # Calculate area of each patch  in Ha
  Ai <- (calculatePatchSizes(patchm) * res(habsuit)[1]^2) / (100*100) # convert to Ha
  if (class(Oj) == "logical") Oj <- rep(1, length(Ai))
  # Calculate connectvity of occupied patches
  S <- calculateConnectivity(Ai, Oj, dij, alpha=alpha)
  # return as SpatRaster
  if (asraster) S <- .rast(assignPatchValues(patchm,S), habsuit)
  return(S)
}
#' @title  Simulates metapopulation dynamics in a real landscape
#' @description Runs a stochastic patch occupancy model given a set of metapopulation
#' input parameters
#' @param habsuit a SpatRaster of habitat suitability with unsuitable areas set to NA
#' @param mu extinction probability in patch of unit carrying capacity.
#' @param x parameter scaling patch carrying capacity to extinction
#' @param alpha inverse of average dispersal distance (km)
#' @param ygamma parameter scaling colonisation probability to connectivity
#' @param aos parameter scaling `habsuit` logistically transformed to population density (individuals per ha)
#' @param timestep number of tiem-steps over which to run the model.
#' @param minden minimum possible density value (individuals per Ha)
#' @param maxden maximum possible density value (individuals per Ha)
#' @param asprob optional logical indicating whether to return expected probability of
#' occurrence (TRUE) or presence absence (FALSE)
#' @return A SpatRaster of presences / absences or occurrence probabilities.
#' @rdname MetaPopSim
#' @import terra
#' @importFrom Rcpp sourceCpp
#' @useDynLib PopScape, .registration = TRUE
#' @export
MetaPopSim <- function(habsuit, mu, x, alpha, ygamma, aos, timesteps = 100, minden = 0.01, maxden = 20, asprob = TRUE) {
  # Calculate habitat patches
  habsuit[habsuit == 0]<-NA
  hab <- habsuit * 0 + 1
  patch<- patches(hab, directions = 8)
  patchm <- as.matrix(patch, wide = TRUE)
  patchm <- renumberPatchIds(patchm)
  patchm[patchm==0]<-NA
  # Create density layer
  rden <- log(habsuit / (1- habsuit))^aos
  rden[is.na(rden)] <- minden
  rden[rden < minden] <- minden
  rden[rden > maxden] <- maxden
  rdenm <-as.matrix(rden, wide = TRUE)
  dij<-patchdistcpp(patchm, rdenm) * res(patch)[1] / 1000 # convert to km
  Oij <- runmetapopmodel(patchm, rdenm, res(patch)[1], dij, mu, x, alpha, ygamma, timesteps)
  if (asprob) {
    # Calculate area
    Ai <- (calculatePatchSizes(patchm) * res(patch)[1]^2) / (100*100) # convert to Ha
    A <- .rast(assignPatchValues(patchm, Ai), patch)
    # Calculate extinction probability
    Ei <- mu / (A * rden)^x
    Ei[Ei > 1] <-1
    # Calculate colonisation probability
    Ci <- ColProb(Ai, Oij, dij, alpha, ygamma)
    Ci <- .rast(assignPatchValues(patchm, Ci), patch)
    occo <- Ci / (Ci + Ei - Ci * Ei)
  } else {
    occo <- .rast(assignPatchValues(patchm, Oij), patch)
  }
  return(occo)
}
#' @title Extracts patch metrics to data.frame of occurrence data
#' @description Extracts habitat patch size, connectivity, landcover type and other relevent
#' environmental data from provided SpatRasters. Called repeatedly by [fitmetapopsdm()]
#' @param occdata a data.frame of occurrence data with x & y coordinates as returned
#' by [subsample()].
#' @param landcover a SpatRaster landcover classes.
#' @param conpredictors a single or multi-layer SpatRaster of other continuous predictors
#' @param habsuit optionally - a SpatRaster of predicted occupancy probability. Set to
#' NA on first iteration.
#' @param Oj a vector of 0s and 1s indicating which patches are occupied - set to NA
#' on first iteration.
#' @param thresh threshold of suitability below which patch density is assumed to equal minden
#' @param alpha inverse of average dispersal distance (km) - used to calculate connectivity
#' @param aos parameter scaling logit `habsuit` to population density
#' @param minden minimum potential density in suitable habitat (individuals per hectare)
#' @param maxden maximum density (individuals per hectare)
#' @return a data.frame matching `occdata` but with the following columns added:
#' \describe{
#'  \item{habclass}{land cover type extracted from landcover}
#'  \item{alpha}{inverse of average dispersal distance (same units as x and y)}
#'  \item{Area}{Patch area in ha}
#'  \item{Connectivity}{Connectivity to neighboring occupied patches}
#'  \item{names(conpredictors)}{Additional data extracted from conpredictors}
#' }
#' @importFrom Rcpp sourceCpp
#' @useDynLib PopScape, .registration = TRUE
#' @import terra
#' @export
extractpatchmetrics <- function(occdata, landcover, conpredictors, habsuit = NA, Oj = NA,
                                thresh = 0.5, alpha = 10, aos = 0.1, minden = 0.01,
                                maxden = 20) {
  # first iteration
  if (class(habsuit) == "logical") {
    # extract land cover type
    occdata$habclass <- extract(landcover, data.frame(x = occdata$x, y = occdata$y))[,2]
    # extract continious environmental predictors
    nlayers <- dim(conpredictors)[3]
    nms <- names(conpredictors)
    n <- dim(occdata)[1]
    mout <- matrix(0, nrow = n, ncol = nlayers)
    for (i in 1:nlayers) {
      r <- conpredictors[[i]]
      mout[,i] <- extract(r, data.frame(x = occdata$x, y = occdata$y))[,2]
    }
    mout <- as.data.frame(mout)
    names(mout)<-nms
    occdata<-cbind(occdata, mout)
    # work out which habitat types have presences or absences
    suithabs <- unique(occdata$habclass)
    suithabs <- suithabs[is.na(suithabs) == FALSE]
    hab <- landcover * NA
    for (i in 1:length(suithabs)) hab[landcover == suithabs[i]] <- 1
    # work out which areas are wihin range of continious predictors
    for (i in 1:nlayers) {
      v <- mout[,i]
      mn <- min(v, na.rm = TRUE)
      mx <- max(v, na.rm = TRUE)
      msk <- conpredictors[[i]]
      msk[msk > mx] <- NA
      msk[msk < mn] <- NA
      hab <- mask(hab, msk)
    }
    habsuit <- hab * 0 + (1 + thresh) / 2
    Ai <- calcarea(habsuit, asraster = FALSE)
    Oj <- Ai * 0 + 1
  }
  # extract patch area
  Ar <- calcarea(habsuit)
  occdata$Area <- extract(Ar, data.frame(x = occdata$x, y = occdata$y))[,2]
  # extract patch density
  rden <- log(habsuit / (1- habsuit))^aos
  rden[is.na(rden)] <- minden
  rden[rden < minden] <- minden
  rden[rden > maxden] <- maxden
  occdata$Density <- extract(rden, data.frame(x = occdata$x, y = occdata$y))[,2]
  # extract patch connectivity
  Sr <- calcconectivity(habsuit, thresh, alpha, aos, minden, maxden, Oj)
  occdata$Connectivity<-extract(Sr, data.frame(x = occdata$x, y = occdata$y))[,2]
  return(occdata)
}
#' @title  Fit hybrid metapopulation-species distribution model
#' @description Iteratively fits a hybrid incidentfunction-species distribution model
#' @param occdata a data.frame of occurrence data with x & y coordinates as returned
#' by [subsample()].
#' @param landcover a SpatRaster landcover classes.
#' @param conpredictors a single or multi-layer SpatRaster of other continuous predictors
#' @param alpha inverse of average dispersal distance (km) - used to calculate connectivity
#' @param aos parameter scaling `habsuit` logistically transformed to population density (individuals per ha)
#' @param bwgt back-weighting used during iteration process in range 0 to <1. Higher = slower, but more likely to converge
#' @param tol tolerence limited for convergence (Maximum difference in predicted probability of occurance between sucessive iterations)
#' @param maxiter maximum number of iterations - will return results
#' @param minden minimum potential density in suitable habitat (individuals per hectare)
#' @param maxden maximum density (individuals per hectare)
#' @param Bayesian optional logical indicating whether to use a normal glm (FALSE - much faster)
#' or constrain the x parameter by setting Bayesian priors (TRUE - much slower)
#' @return a list of two SpatRasters: (1) habsuit- the habitat suitability and (2) pocco -
#' estimated probability of occurence.
#' @importFrom Rcpp sourceCpp
#' @useDynLib PopScape, .registration = TRUE
#' @import terra
#' @export
fitmetapopsdm <- function(occdata, landcover, conpredictors, alpha = 10, aos = 0.1,
                          bwgt = 0, tol = 0.01, maxiter = 100, minden = 0.01, maxden = 20, Bayesian = FALSE) {
  # first iteration
  thresh <- 0.5  # can be adjusted if needed
  names(conpredictors) <- "lveghgt"
  occdata <- extractpatchmetrics(occdata, landcover, conpredictors, habsuit = NA, Oj = NA,
                                 thresh, alpha, aos, minden, maxden)
  occdata$Area[occdata$Area == 0] <- 0.000001
  habsuit <- landcover * NA
  lc <- unique(occdata$habclass )
  habsuit[landcover == lc[1]] <- 0.75
  habsuit[landcover == lc[2]] <- 0.75
  mn <- min(occdata$lveghgt)
  mx <- max(occdata$lveghgt)
  habsuit[conpredictors > mx] <- NA
  habsuit[conpredictors < mn] <- NA
  Ai <- calcarea(habsuit, asraster = FALSE)
  Oj <- Ai * 0 + 1
  Ar <- calcarea(habsuit)
  Sr <- calcconectivity(habsuit, thresh, alpha, aos, minden, maxden, Oj)
  preds <- .fitone(occdata, landcover, conpredictors, habsuit, Ar, Sr, aos, thresh, minden, maxden)
  habsuit <- preds$habsuit
  pocco <- preds$pocco
  # recalculate patch metrics
  hab <- habsuit * 0 + 1
  hab[habsuit < 0.5] <- 0
  patch<- patches(hab, directions = 8)
  patchm <- as.matrix(patch, wide = TRUE)
  patchm <- renumberPatchIds(patchm)
  patchm[patchm==0]<-NA
  poccopatch <- calculateAveragePatchQuality(patchm, as.matrix(pocco, wide = TRUE))
  Oj <- ifelse(poccopatch < thresh, 0, 1)
  occdata <- extractpatchmetrics(occdata, landcover, conpredictors, habsuit, Oj,
                                 thresh, alpha, aos, minden, maxden)
  occdata$Area[occdata$Area == 0] <- 0.000001
  # set error etc
  errr <- tol * 2
  itr = 0
  while (errr < tol) {
    preds <- .fitone(occdata, landcover, conpredictors, habsuit, Ar, Sr, aos, thresh, minden, maxden)
    habsuit <- bwgt * habsuit + (1 - bwgt) * preds$habsuit
    nocco <- bwgt * pocco + (1 - bwgt) * preds$pocco
    dif <- abs(nocco - pocco)
    plot(dif)
    pocco <- nocco
    errr <- max(as.vector(dif), na.rm = TRUE)
    itr <- itr + 1
    if (itr == 100) errr <- 0
    # recalculate patch metrics
    hab <- habsuit * 0 + 1
    hab[habsuit < thresh] <- 0
    patch<- patches(hab, directions = 8)
    patchm <- as.matrix(patch, wide = TRUE)
    patchm <- renumberPatchIds(patchm)
    patchm[patchm==0]<-NA
    poccopatch <- calculateAveragePatchQuality(patchm, as.matrix(pocco, wide = TRUE))
    Oj <- ifelse(poccopatch < thresh, 0, 1)
    occdata <- extractpatchmetrics(occdata, landcover, conpredictors, habsuit, Oj,
                                   thresh, alpha, aos, minden, maxden)
    occdata$Area[occdata$Area == 0] <- 0.000001
  }
  # derive metapopulation parameters
  occdata$ccapacity <- occdata$Area * occdata$Density
  mod <- suppressWarnings(glm(pa ~ offset(2 *log(Connectivity)) + log(ccapacity) +
                                as.factor(habclass) + lveghgt, family = binomial, data = occdata))
  be <- coef(mod)
  x <- as.numeric(be[2])
  A0 <- min(occdata$Area[occdata$pa > 0])
  muga <- as.numeric(exp(-be[1]))
  mu <- A0^x
  ga <- muga/mu
  params <- list(mu = mu, x = x, alpha = alpha, ygamma = ga, aos = aos)
  return(list(params = params, pocco = pocco))
}
