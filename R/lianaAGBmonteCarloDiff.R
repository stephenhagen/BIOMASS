library("VGAM")
#' Propagating liana above-ground biomass (AGB) or carbon (AGC) errors to the stand level
#'
#' Propagation of the errors throughout the steps needed to compute AGB or AGC.
#'
#' @param Dsizebins Vector of liana counts by size bin 
#' @param sizeClassLow Vector of low end of liana size bins (in cm)
#' @param sizeClassHigh Vector of high end of liana size bins (in cm)
#' @param n Number of iterations. Cannot be smaller than 50 or larger than 1000. By default `n = 1000`
#' @param Carbon (logical) Whether or not the propagation should be done up to the carbon value (FALSE by default).
#' @param Dlim (optional) Minimum diameter (in cm) for which above-ground biomass should be calculated (all diameter below
#' `Dlim` will have a 0 value in the output).
#' @param plot (optional) Plot ID, must be either one value, or a vector of the same length as D. This argument is used to build
#' stand-specific HD models.
#'
#' @details See Rejou-Mechain et al. (2017) for all details on the error propagation procedure.
#'
#' @return Returns a list  with (if Carbon is FALSE):
#'   - `meanAGB`: Mean stand AGB value following the error propagation
#'   - `medAGB`: Median stand AGB value following the error propagation
#'   - `sdAGB`: Standard deviation of the stand AGB value following the error propagation
#'   - `credibilityAGB`: Credibility interval at 95\% of the stand AGB value following the error propagation
#'   - `AGB_simu`: Matrix with the AGB of the trees (rows) times the n iterations (columns)
#'
#' @references Chave, J. et al. (2004). _Error propagation and scaling for tropical forest biomass estimates_.
#' Philosophical Transactions of the Royal Society B: Biological Sciences, 359(1443), 409-420.
#' @references Rejou-Mechain et al. (2017).
#' _BIOMASS: An R Package for estimating above-ground biomass and its uncertainty in tropical forests_.
#' Methods in Ecology and Evolution, 8 (9), 1163-1167.
#'
#' @author Maxime REJOU-MECHAIN, Bruno HERAULT, Camille PIPONIOT, Ariane TANGUY, Arthur PERE
#'
#' @examples
#' # Load a database
#' data(NouraguesHD)
#' data(KarnatakaForest)
#'
#' # Modelling height-diameter relationship
#' HDmodel <- modelHD(D = NouraguesHD$D, H = NouraguesHD$H, method = "log2")
#'
#' # Retrieving wood density values
#' \donttest{
#' KarnatakaWD <- getWoodDensity(KarnatakaForest$genus, KarnatakaForest$species,
#'   stand = KarnatakaForest$plotId
#' )
#' }
#'
#' # Propagating errors with a standard error in wood density in one plot
#' filt <- KarnatakaForest$plotId == "BSP20"
#' set.seed(10)
#' \donttest{
#' resultMC <- AGBmonteCarlo(
#'   D = KarnatakaForest$D[filt], WD = KarnatakaWD$meanWD[filt],
#'   errWD = KarnatakaWD$sdWD[filt], HDmodel = HDmodel
#' )
#' str(resultMC)
#' }
#'
#' # If only the coordinates are available
#' lat <- KarnatakaForest$lat[filt]
#' long <- KarnatakaForest$long[filt]
#' coord <- cbind(long, lat)
#' \donttest{
#' resultMC <- AGBmonteCarlo(
#'   D = KarnatakaForest$D[filt], WD = KarnatakaWD$meanWD[filt],
#'   errWD = KarnatakaWD$sdWD[filt], coord = coord
#' )
#' str(resultMC)
#' }
#'
#' # Propagating errors with a standard error in wood density in all plots at once
#' \donttest{
#' KarnatakaForest$meanWD <- KarnatakaWD$meanWD
#' KarnatakaForest$sdWD <- KarnatakaWD$sdWD
#' resultMC <- by(
#'   KarnatakaForest, KarnatakaForest$plotId,
#'   function(x) AGBmonteCarlo(
#'       D = x$D, WD = x$meanWD, errWD = x$sdWD,
#'       HDmodel = HDmodel, Dpropag = "chave2004"
#'     )
#' )
#' meanAGBperplot <- unlist(sapply(resultMC, "[", 1))
#' credperplot <- sapply(resultMC, "[", 4)
#' }
#'
#' @keywords monte carlo
#' @importFrom stats pnorm qnorm runif
#' @export

lianaAGBmonteCarloDiff <- function(Dsizebins, sizeClassLow, sizeClassHigh, Dsizebins2 = NULL, 
                                   n = 1000, Carbon = FALSE, Dlim = NULL, plot = NULL,
                                   useMeasurementError = TRUE, useModelError = TRUE) {

  # parameters verification -------------------------------------------------

  if (n > 1000 | n < 50) {
    stop("n cannot be smaller than 50 or larger than 1000")
  }

  # function truncated random gausien law -----------------------------------
  myrtruncnorm <- function(n, lower = -1, upper = 1, mean = 0, sd = 1) {
    qnorm(runif(n, pnorm(lower, mean = mean, sd = sd), pnorm(upper, mean = mean, sd = sd)), mean = mean, sd = sd)
  }

  # adjust the counts in each size class due to (a) miscounting and (b) mis-sizing
  # uses a double exponential random variable; adjust this by comparing the liana counts of the same plot
  #   conducted by two separate measurement teams
  adjustBins <- function(Dsizebins) {
    scaleFactor <- 2.0 # make this dependent on plot size
    adjustArray <- round(rlaplace(length(Dsizebins), location = 0, scale = 1)/scaleFactor)
    DsizebinsAdj <- Dsizebins + adjustArray
    DsizebinsAdj[DsizebinsAdj < 0] <- 0
    return(DsizebinsAdj)
  }

  D_mat = array(NA, dim=c(round(sum(Dsizebins)*1.2),n))
  maxLen = 0
  if (!is.null(Dsizebins2)) {
    D_mat2 = array(NA, dim=c(round(sum(Dsizebins2)*1.2),n))
    maxLen2 <- 0
  }
  for (j in 1:n) {
    if (useMeasurementError == FALSE) {
      DsizebinsAdj <- Dsizebins
      D <- c()
      for (i in 1:length(sizeClassLow)) {
        tempLianas <- rep(mean(sizeClassLow[i], sizeClassHigh[i]), DsizebinsAdj[i])
        D <- append(D, tempLianas)
      }
    } else {
      DsizebinsAdj <- adjustBins(Dsizebins)
      # change from count in size class bins to a list of all liana stems, and use a random uniform
      #   draw from the size bin range
      D <- c()
      for (i in 1:length(sizeClassLow)) {
        tempLianas <- runif(DsizebinsAdj[i], min=sizeClassLow[i], max=sizeClassHigh[i])
        D <- append(D, tempLianas)
      }
    }
    D_mat[1:length(D),j] = D
    if (length(D) > maxLen) {
      maxLen <- length(D)
    }
    if (!is.null(Dsizebins2)) {
      if (useMeasurementError == FALSE) {
        DsizebinsAdj2 <- Dsizebins2
        D2 <- c()
        for (i in 1:length(sizeClassLow)) {
          tempLianas <- rep(mean(sizeClassLow[i], sizeClassHigh[i]), DsizebinsAdj2[i])
          D2 <- append(D2, tempLianas)
        }
      } else {
        DsizebinsAdj2 <- adjustBins(Dsizebins2)
        # change from count in size class bins to a list of all liana stems, and use a random uniform
        #   draw from the size bin range
        D2 <- c()
        for (i in 1:length(sizeClassLow)) {
          tempLianas <- runif(DsizebinsAdj2[i], min=sizeClassLow[i], max=sizeClassHigh[i])
          D2 <- append(D2, tempLianas)
        }
      }
      D_mat2[1:length(D2),j] = D2
      if (length(D2) > maxLen2) {
        maxLen2 <- length(D2)
      }
    }
  }
  D_simu <- D_mat[1:maxLen,] 
  len <- maxLen
  if (!is.null(Dsizebins2)) {
    D_simu2 <- D_mat2[1:maxLen2,] 
    len2 <- maxLen2
  }
  ### Propagate error with Markov Chain Monte Carlo approach

  # --------------------- D ---------------------

#  D_simu <- suppressWarnings(replicate(n, myrtruncnorm(len, mean = D, sd = Dpropag, lower = 0.1, upper = 500)))


  # --------------------- AGB ---------------------

  if (useModelError == TRUE) {
    param_4 <- BIOMASS::param_4liana
    selec <- sample(1:nrow(param_4), n)
    RSE <- param_4[selec, "sd"]

    # Posterior model parameters
    Ealpha <- param_4[selec, "intercept"]
    Ebeta <- param_4[selec, "logagbt"]

    # Construct a matrix where each column contains random errors taken from N(0,RSEi) with i varying between 1 and n
    matRSE <- mapply(function(y) {
      rnorm(sd = y, n = len)
    }, y = RSE)
  } else {
    matRSE <- 0
    Ebeta <- 2.556
    Ealpha <- -1.431
  }

  # Propagation of the error using simulated parameters
  Comp <- t(log(D_simu)) * Ebeta + Ealpha
  Comp <- t(Comp) + matRSE

  # Backtransformation
  AGB_simu <- exp(Comp) / 1000

  if (!is.null(Dlim)) AGB_simu[D < Dlim, ] <- 0
  AGB_simu[ which(is.infinite(AGB_simu)) ] <- NA


  if (!is.null(Dsizebins2)) {
# THIS ERROR NEEDS TO BE GENERATED SEPERATELY FOR TIME 2 - THIS IS DIFFERENT THAN THE TREE UNCERTAINTY; INCONSISTENT; FIX
    if (useModelError == TRUE) {
      # Construct a matrix where each column contains random errors taken from N(0,RSEi) with i varying between 1 and n
      matRSE2 <- mapply(function(y) {
        rnorm(sd = y, n = len2)
      }, y = RSE)
    } else {
      matRSE2 <- 0
      Ebeta <- 2.556
      Ealpha <- -1.431
    }

    # Propagation of the error using simulated parameters
    Comp2 <- t(log(D_simu2)) * Ebeta + Ealpha
    Comp2 <- t(Comp2) + matRSE2

    # Backtransformation
    AGB_simu2 <- exp(Comp2) / 1000

    if (!is.null(Dlim)) AGB_simu2[D2 < Dlim, ] <- 0
    AGB_simu2[ which(is.infinite(AGB_simu2)) ] <- NA
  }

  if (Carbon == FALSE) {
    sum_AGB_simu <- colSums(AGB_simu, na.rm = TRUE)
    res1 <- list(
      meanAGB = mean(sum_AGB_simu),
      medAGB = median(sum_AGB_simu),
      sdAGB = sd(sum_AGB_simu),
      credibilityAGB = quantile(sum_AGB_simu, probs = c(0.025, 0.975)),
      sum_AGB_simu = sum_AGB_simu,
      AGB_simu = AGB_simu
    )
    if (!is.null(Dsizebins2)) {
      sum_AGB_simu2 <- colSums(AGB_simu2, na.rm = TRUE)
      sum_Diff_AGB_simu <- sum_AGB_simu2 - sum_AGB_simu
      res2 <- list(
        meanAGB2 = mean(sum_AGB_simu2),
        medAGB2 = median(sum_AGB_simu2),
        sdAGB2 = sd(sum_AGB_simu2),
        credibilityAGB2 = quantile(sum_AGB_simu2, probs = c(0.025, 0.975)),
        sum_AGB_simu2 = sum_AGB_simu2,
        AGB_simu2 = AGB_simu2,
        
        meanDiff_AGB = mean(sum_Diff_AGB_simu),
        medDiff_AGB = median(sum_Diff_AGB_simu),
        sdDiff_AGB = sd(sum_Diff_AGB_simu),
        credibilityDiff_AGB = quantile(sum_Diff_AGB_simu, probs = c(0.025, 0.975)),
        sum_Diff_AGB_simu = sum_Diff_AGB_simu
      )
    }
  } else {
    # Biomass to carbon ratio calculated from Thomas and Martin 2012 forests data stored in DRYAD (tropical
    # angiosperm stems carbon content)
    AGC_simu <- AGB_simu * rnorm(mean = 47.13, sd = 2.06, n = n * len) / 100
    sum_AGC_simu <- colSums(AGC_simu, na.rm = TRUE)
    res1 <- list(
      meanAGC = mean(sum_AGC_simu),
      medAGC = median(sum_AGC_simu),
      sdAGC = sd(sum_AGC_simu),
      credibilityAGC = quantile(sum_AGC_simu, probs = c(0.025, 0.975)),
      sum_AGC_simu = sum_AGC_simu,
      AGC_simu = AGC_simu
    )
    if (!is.null(Dsizebins2)) {
      AGC_simu2 <- AGB_simu2 * rnorm(mean = 47.13, sd = 2.06, n = n * len2) / 100
      sum_AGC_simu2 <- colSums(AGC_simu2, na.rm = TRUE)
      sum_Diff_AGC_simu <- sum_AGC_simu2 - sum_AGC_simu
      res2 <- list(
        meanAGC2 = mean(sum_AGC_simu2),
        medAGC2 = median(sum_AGC_simu2),
        sdAGC2 = sd(sum_AGC_simu2),
        credibilityAGC2 = quantile(sum_AGC_simu2, probs = c(0.025, 0.975)),
        sum_AGC_simu2 = sum_AGC_simu2,
        AGC_simu2 = AGC_simu2,

        meanDiff_AGC = mean(sum_Diff_AGC_simu),
        medDiff_AGC = median(sum_Diff_AGC_simu),
        sdDiff_AGC = sd(sum_Diff_AGC_simu),
        credibilityDiff_AGC = quantile(sum_Diff_AGC_simu, probs = c(0.025, 0.975)),
        sum_Diff_AGC_simu = sum_Diff_AGC_simu
      )
    }
  }
  res <- c(res1, res2)
  return(res)
}
