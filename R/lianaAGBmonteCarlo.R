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

lianaAGBmonteCarlo <- function(Dsizebins, sizeClassLow, sizeClassHigh, n = 1000, Carbon = FALSE, Dlim = NULL, plot = NULL) {

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

  D_mat = array(NA, dim=c(round(sum(Dsizebins)*1.1),n))
  maxLen = 0
  for (j in 1:n) {
    DsizebinsAdj <- adjustBins(Dsizebins)
    # change from count in size class bins to a list of all liana stems, and use a random uniform
    #   draw from the size bin range
    D <- c()
    for (i in 1:length(sizeClassLow)) {
      tempLianas <- runif(DsizebinsAdj[i], min=sizeClassLow[i], max=sizeClassHigh[i])
      D <- append(D, tempLianas)
    }
    D_mat[1:length(D),j] = D
    if (length(D) > maxLen) {
      maxLen <- length(D)
    }
  }
  D_simu <- D_mat[1:maxLen,] 
  len <- maxLen
  ### Propagate error with Markov Chain Monte Carlo approach

  # --------------------- D ---------------------

#  D_simu <- suppressWarnings(replicate(n, myrtruncnorm(len, mean = D, sd = Dpropag, lower = 0.1, upper = 500)))


  # --------------------- AGB ---------------------

  param_4 <- BIOMASS::param_4liana
  selec <- sample(1:nrow(param_4), n)
  RSE <- param_4[selec, "sd"]

  # Construct a matrix where each column contains random errors taken from N(0,RSEi) with i varying between 1 and n
  matRSE <- mapply(function(y) {
    rnorm(sd = y, n = len)
  }, y = RSE)

  # Posterior model parameters
  Ealpha <- param_4[selec, "intercept"]
  Ebeta <- param_4[selec, "logagbt"]

  # Propagation of the error using simulated parameters
  Comp <- t(log(D_simu)) * Ebeta + Ealpha
  Comp <- t(Comp) + matRSE

  # Backtransformation
  AGB_simu <- exp(Comp) / 1000

  if (!is.null(Dlim)) AGB_simu[D < Dlim, ] <- 0
  AGB_simu[ which(is.infinite(AGB_simu)) ] <- NA

  if (Carbon == FALSE) {
    sum_AGB_simu <- colSums(AGB_simu, na.rm = TRUE)
    res <- list(
      meanAGB = mean(sum_AGB_simu),
      medAGB = median(sum_AGB_simu),
      sdAGB = sd(sum_AGB_simu),
      credibilityAGB = quantile(sum_AGB_simu, probs = c(0.025, 0.975)),
      AGB_simu = AGB_simu
    )
  } else {
    # Biomass to carbon ratio calculated from Thomas and Martin 2012 forests data stored in DRYAD (tropical
    # angiosperm stems carbon content)
    AGC_simu <- AGB_simu * rnorm(mean = 47.13, sd = 2.06, n = n * len) / 100
    sum_AGC_simu <- colSums(AGC_simu, na.rm = TRUE)
    res <- list(
      meanAGC = mean(sum_AGC_simu),
      medAGC = median(sum_AGC_simu),
      sdAGC = sd(sum_AGC_simu),
      credibilityAGC = quantile(sum_AGC_simu, probs = c(0.025, 0.975)),
      AGC_simu = AGC_simu
    )
  }
  return(res)
}
