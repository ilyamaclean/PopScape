#' converts matrix to SpatRaster using template
#' @import terra
#' @export
.rast <- function(m, rtemplate) {
  r <- rast(m)
  ext(r) <- ext(rtemplate)
  crs(r) <- crs(rtemplate)
  return(r)
}
#' @title  Extracts sub-sample of points from simulated metapopulation model output
#' @description Extracts a sub-sample of points from simulated metapopulation model output
#' for subsequent testing
#' @param occo SpatRaster of predicted probability of occupancy as returned by [MetaPopSim()]
#' @param landcover SpatRaster of land cover type
#' @param n number of points ot extract
#' @return a data.frame with the following columns
#' \describe {
#'  \item{x}{randomly generated x-coordinate}
#'  \item{y}{randomly generated y-coordinate}
#'  \item{pa}{presence or absence at location given by x and y}
#' }
#' @rdname subsample
#' @import terra
#' @export
subsample<-function(occo, landcover, n) {
  # calculate which land cover types have predicted occupancy greater than 0
  occo[occo == 0] <- NA
  landcovers <- mask(landcover, occo)
  suithabs <- unique(as.vector(landcovers))
  suithabs <- suithabs[is.na(suithabs) == FALSE]
  hab <- landcover * NA
  for (i in 1:length(suithabs)) hab[landcover == suithabs[i]] <- 1
  landcover <- mask(landcover, hab)
  # calculate number of datapoints ot extract to ensure enough if dataset has NAs
  v <- as.vector(landcover)
  mu <- length(v) / length(v[is.na(v)==FALSE])
  n2 <- trunc(n * mu *2) # twice as many as needed to ensure definately enough
  e <- ext(occo)
  # generate random xs and ys
  x <- runif(n2, e$xmin, e$xmax)
  y <- runif(n2, e$ymin, e$ymax)
  pa <- extract(occo, data.frame(x,y))[,2]
  pa <- ifelse(pa < 0.5, 0, 1)
  s <- which(is.na(pa) == FALSE)
  s2 <- round(runif(n * 2,1,length(s)),0)[1:n]
  s <- s[s2]
  dfout <- data.frame(x = x[s],y = y[s], pa = pa[s][1:n])
  return(dfout)
}
#' Fits one iteration of hybrid metapopulation-sdm model
.fitone <- function(occdata, landcover, conpredictors, habsuit, Ar, Sr, aos = 0.1, thresh = 0.5, minden = 0.01, maxden = 20, Bayesian = FALSE) {
  # Calculate carrying capacity
  occdata$ccapacity <- occdata$Area * occdata$Density
  # fit model
  if (Bayesian) {
    mod <- brm(
      pa ~ offset(2 * log(Connectivity)) + log(ccapacity) + as.factor(habclass) + lveghgt,
      family = bernoulli(link = "logit"),
      data = occdata,
      prior = c(
        prior(normal(0.5, 1), class = "b", coef = "logccapacity"),
        prior(uniform(-10, 10), class = "b")  # Default priors for other coefficients
      )
    )
    be<-as.numeric(fixef(mod)[,1])
  }  else {
    mod <- suppressWarnings(glm(pa ~ offset(2 *log(Connectivity)) + log(ccapacity) +
                                  as.factor(habclass) + lveghgt, family = binomial, data = occdata))
    be<-coef(mod)
  }
  # identify the two landcover types
  lcs <- unique(occdata$habclass)
  # derive habitat suitability
  lhabsuit <- be[1] + 2 * mean(log(occdata$Connectivity)) + be[2] * mean(log(occdata$ccapacity)) +
    be[4] * conpredictors
  lhabsuit[landcover == lcs[2]] <- lhabsuit[landcover == lcs[2]] + be[3]
  lhabsuit[landcover != lcs[1] & landcover != lcs[2]] <- NA
  # derive occurence prediction
  densty <- lhabsuit^aos
  densty[is.na(densty)] <- minden
  densty[densty < minden] <- minden
  densty[densty > maxden] <- maxden
  lpocc <- be[1] + 2 * log(Sr) + be[2] * mean(log(densty * Ar)) +
    be[4] * conpredictors
  lpocc[landcover == lcs[1]] <- lpocc[landcover == lcs[1]] + be[3]
  # perform prevailance adjustment
  prev <- sum(occdata$pa) / length(occdata$pa)
  lprev <- log(prev / (1 - prev))
  toadd <- lprev - mean(as.vector(lpocc), na.rm = TRUE)
  lpocc <- lpocc + toadd
  lhabsuit <- lhabsuit + toadd
  # Inverse logit transofrm
  habsuit <- 1 / (1+exp(-lhabsuit))
  pocco <- 1 / (1+exp(-lpocc))
  return(list(habsuit = habsuit, pocco= pocco))
}
