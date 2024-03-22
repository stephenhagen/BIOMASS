#' Propagating above-ground biomass (AGB) or carbon (AGC) errors to the stand level
#'
#' Propagation of the errors throughout the steps needed to compute AGB or AGC.
#'
#' @param D Vector of tree diameters (in cm)
#' @param WD Vector of wood density estimates (in g/cm3)
#' @param errWD Vector of error associated to the wood density estimates (should be of the same size as `WD`)
#' @param H (option 1) Vector of tree heights (in m). If set, `errH` must be set too.
#' @param errH (if `H`) Residual standard error (RSE) of a model or vector of errors (sd values) associated to tree height
#' values (in the latter case the vector should be of the same length as `H`).
#' @param HDmodel (option 2) Model used to estimate tree height from tree diameter (output from [modelHD()], see example).
#' @param coord (option 3) Coordinates of the site(s), either a vector giving a single site (e.g. c(longitude, latitude))
#' or a matrix/dataframe with two columns (e.g. cbind(longitude, latitude)). The coordinates are used to predict
#' height-diameter allometry with bioclimatic variables.
#' @param Dpropag This variable can take three kind of values, indicating how to propagate the errors on diameter measurements:
#' a single numerical value or a vector of the same size as `D`, both representing the standard deviation associated
#' with the diameter measurements or `"chave2004"` (an important error on 5 percent of the measures, a smaller error on
#' 95 percent of the trees).
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


# WHEN USING REAL DATA AT TWO TIMES, D AND D2 WILL HAVE TO BE OF THE SAME LENGTH, EVEN WITH
#    MORTALITY AND RECRUITMENT. DBH WILL 0 IN TIME 2 FOR MORTALITY. DBH WILL BE 0 IN TIME 1 
#    FOR RECRUITMENT.
AGBmonteCarloDiff <- function(D, D2 = NULL, WD = NULL, errWD = NULL, H = NULL, errH = NULL,
                          HDmodel = NULL, coord = NULL, Dpropag = NULL, n = 1000,
                          Carbon = FALSE, Dlim = NULL, plot = NULL, Dthresh = NULL, 
                          paramAGBread = NULL, paramAGBwrite = NULL) {
  len <- length(D)

  # parameters verification -------------------------------------------------

  if (n > 1000 | n < 50) {
    stop("n cannot be smaller than 50 or larger than 1000")
  }

  if (!is.null(Dpropag)) {
    if ((is.numeric(Dpropag) && !(length(Dpropag) %in% c(1, len)) || (!is.numeric(Dpropag) && tolower(Dpropag) != "chave2004"))) {
      stop('Dpropag should be set to one of these options:
             - "chave2004"
             - a single sd value that will be applied to all trees
             - a vector of sd values of the same length as D')
    }
  }

  if (is.null(WD) || is.null(errWD)) {
    stop("The WD and errWD arguments must be not NULL")
  }

  if (len != length(WD) || len != length(errWD)) {
    stop("One of vector WD or errWD does not have the same length as D")
  }

  if (is.null(HDmodel) & is.null(coord) & is.null(H)) {
    stop("Input missing, you need to provide one of the following arguments:
             - H
             - HDmodel
             - coord")
  }

  if ((!is.null(HDmodel) && !is.null(coord)) || (!is.null(HDmodel) && !is.null(H)) || (!is.null(coord) && !is.null(H))) {
    stop("Too many input, choose one input among those arguments:
              - H and Herr
              - HDmodel
              - coord")
  }

  if (!is.null(H)) {
    if (is.null(errH)) {
      stop("Cannot propagate height errors without information on associated errors (errH is null),
         if you do not want to propagate H errors please set errH to 0")
    }
    if (length(H) != len || !(length(errH) %in% c(1, len))) {
      stop("H must be the same length as D and errH must be either one value or the same length as D")
    }
  }

  if (!is.null(coord) && ((is.vector(coord) && length(coord) != 2) || (is.matrix(coord) && nrow(coord) != len))) {
    stop("coord should be either
             - a vector (e.g. c(longitude, latitude))
             - a matrix with two columns (longitude and latitude)
             having the same number of rows as the number of trees (length(D))")
  }

  # the length of the plot is tested in predictHeight
  # the names of the plot and the names of the model is tested in predictHeight
  if (!is.null(plot) && is.null(HDmodel)) {
    stop("The 'plot' vector must be with 'model' argument")
  }

  if (!is.null(Dthresh) && (Dthresh < 10.0 || Dthresh > 100.0)) {
    stop("The Dthresh dendrometer threshold must be between 10 and 100 cm")
  }


  # function truncated random gausien law -----------------------------------
  myrtruncnorm <- function(n, lower = -1, upper = 1, mean = 0, sd = 1) {
    qnorm(runif(n, pnorm(lower, mean = mean, sd = sd), pnorm(upper, mean = mean, sd = sd)), mean = mean, sd = sd)
  }





  ### Propagate error with Markov Chain Monte Carlo approach

  # --------------------- D ---------------------

  if (!is.null(Dpropag)) {
    if (length(Dpropag) == 1 && tolower(Dpropag) == "chave2004") {
      # Propagation of the measurement error on D: based on Chave et al. 2004 (p.412) Phil. Trans. R. Soc. Lond. B.
      fivePercent <- round(len * 5 / 100)
      chaveError <- function(x, len, Dthresh) {
        ## Use dendrometers on trees with DBH over 30 cm
        #dendro_threshold <- 20.0
        
        ## Assigning large errors on 5% of the trees
        largeErrSample <- sample(len, fivePercent)

        D_sd <- 0.0062 * x + 0.0904 # Assigning small errors on the remaining 95% trees
        # Assume that dendrometers reduce the error - see Butt et al 2013 - https://pubmed.ncbi.nlm.nih.gov/23865241/
#        if (!is.null(Dthresh)) {
#          if (Dthresh >= 10.0 && Dthresh <= 100.0) {
#            D_sd[x > dendro_threshold] <- 0.0062 * x[x > dendro_threshold] + 0.0704 # Assigning smaller errors on dendrobanded trees
#          }
#        }
        D_sd[largeErrSample] <- 4.64

        D_adj <- myrtruncnorm(n = len, mean = x, sd = D_sd, lower = 0.1, upper = 500)
        # Assume the measurement error associated with dendrometers is 0.03 cm lower than with tape 
        #   dendrometers reduce the error - see Butt et al 2013 - https://pubmed.ncbi.nlm.nih.gov/23865241/
        if (!is.null(Dthresh)) {
          D_err <- D_adj - x
          signArray <- D_err
          signArray[signArray > 0] <- 1
          signArray[signArray <= 0] <- -1
          D_err_Den <- abs(D_err) 
          D_err_Den[x >= Dthresh] <- D_err_Den[x >= Dthresh] - 0.031 
          D_err_Den[D_err_Den < 0] <- 0
          D_err_Den[signArray <= 0] <- -1.0 * D_err_Den[signArray <= 0] 
          D_adj <- x + D_err_Den
        }
        return(D_adj)
      }
      D_simu <- suppressWarnings(replicate(n, chaveError(D, len, Dthresh)))
      if (!is.null(D2)) {
        D2_simu <- suppressWarnings(replicate(n, chaveError(D2, len, Dthresh)))
      }
    }
    else {
      D_simu <- suppressWarnings(replicate(n, myrtruncnorm(len, mean = D, sd = Dpropag, lower = 0.1, upper = 500)))
      if (!is.null(D2)) {
        # DELETE THE NEXT LINE IF DBH AT TIME 2 IS AN ACTUAL MEASUREMENT 
#        D2_simu <- suppressWarnings(replicate(n, myrtruncnorm(len, mean = D2, sd = sqrt(2 * Dpropag^2), lower = 0.1, upper = 500)))
        D2_simu <- suppressWarnings(replicate(n, myrtruncnorm(len, mean = D2, sd = Dpropag, lower = 0.1, upper = 500)))
      }
    }
  } else {
    D_simu <- replicate(n, D)
      if (!is.null(D2)) {
        D2_simu <- replicate(n, D2)
      }
  }





  # --------------------- WD ---------------------


  #### Below 0.08 and 1.39 are the minimum and the Maximum WD value from the global wood density database respectively
  WD_simu <- suppressWarnings(replicate(n, myrtruncnorm(n = len, mean = WD, sd = errWD, lower = 0.08, upper = 1.39)))






  # --------------------- H ---------------------

  # if there is data for H
  if (!is.null(HDmodel) | !is.null(H)) {
    if (!is.null(HDmodel)) {
      # Propagation of the error thanks to the local model of H
      H_simu <- apply(D_simu, 2, function(x) predictHeight(x, model = HDmodel, err = TRUE, plot = plot))
      if (!is.null(D2)) {
        H2_simu <- apply(D2_simu, 2, function(x) predictHeight(x, model = HDmodel, err = TRUE, plot = plot))
      }
    } else {
      # Propagation of the error using the errH value(s)
      upper <- max(H, na.rm = TRUE) + 15
      H_simu <- suppressWarnings(replicate(n, myrtruncnorm(len, mean = H, sd = errH, lower = 1.3, upper = upper)))
#      if (!is.null(D2)) {
#        # PROBLEM HERE: NEED TO COME UP WITH A SOLUTION. ANOTHER SET OF HEIGHT MEASUREMENTS?
#      }
    }

    # --------------------- AGB ---------------------

    param_4 <- BIOMASS::param_4
    if (!is.null(paramAGBread)) {
      selec <- read.csv(paramAGBread)
    } else {
      selec <- sample(1:nrow(param_4), n)
    }
    if (!is.null(paramAGBwrite)) {
      write.csv(selec, paramAGBwrite, row.names=FALSE)
    }
    RSE <- param_4[selec, "sd"]

    # Construct a matrix where each column contains random errors taken from N(0,RSEi) with i varying between 1 and n
    matRSE <- mapply(function(y) {
      rnorm(sd = y, n = len)
    }, y = RSE)

    # Posterior model parameters
    Ealpha <- param_4[selec, "intercept"]
    Ebeta <- param_4[selec, "logagbt"]

    # Propagation of the error using simulated parameters
    Comp <- t(log(WD_simu * H_simu * D_simu^2)) * Ebeta + Ealpha
    Comp <- t(Comp) + matRSE

    # Backtransformation
    AGB_simu <- exp(Comp) / 1000

    if (!is.null(D2)) {
      # DO THESE RANDOM ERRORS NEED TO BE ESTIMATED SEPARATELY AT T1 AND T2? ASSUMING NO FOR NOW 
#      matRSE2 <- mapply(function(y) {
#        rnorm(sd = y, n = len)
#      }, y = RSE)
      Comp2 <- t(log(WD_simu * H_simu * D2_simu^2)) * Ebeta + Ealpha
#      Comp2 <- t(Comp2) + matRSE2
      Comp2 <- t(Comp2) + matRSE
      AGB2_simu <- exp(Comp2) / 1000
    }

  }

  # --------------------- Coordinates ---------------------

  # If there is no data for H, but site coordinates
  if (!is.null(coord)) {
    if (is.null(dim(coord))) {
      coord <- as.matrix(t(coord))
    }

    bioclimParams <- getBioclimParam(coord) # get bioclim variables corresponding to the coordinates

    if (nrow(bioclimParams) == 1) {
      bioclimParams <- bioclimParams[rep(1, len), ]
    }

    # Equ 7
    # Log(agb) = -1.803 - 0.976 (0.178TS - 0.938CWD - 6.61PS) + 0.976log(WD) + 2.673log(D) -0.0299log(D2)
    param_7 <- BIOMASS::param_7
#    selec <- sample(1:nrow(param_7), n)
    if (!is.null(paramAGBread)) {
      selec <- read.csv(paramAGBread)
      selec <- selec$x
    } else {
      selec <- sample(1:nrow(param_7), n)
    }
    if (!is.null(paramAGBwrite)) {
      write.csv(selec, paramAGBwrite, row.names=FALSE)
    }

    # Posterior model parameters
    RSE <- param_7[selec, "sd"] # vector of simulated RSE values

    # Recalculating n E values based on posterior parameters associated with the bioclimatic variables
    Esim <- tcrossprod(as.matrix(param_7[selec, c("temp", "prec", "cwd")]), as.matrix(bioclimParams))

      # DO THESE RANDOM ERRORS NEED TO BE ESTIMATED SEPARATELY AT T1 AND T2? ASSUMING NO FOR NOW
#    selec2 <- sample(1:nrow(param_7), n)
#    RSE2 <- param_7[selec2, "sd"] # vector of simulated RSE values
#    Esim2 <- tcrossprod(as.matrix(param_7[selec2, c("temp", "prec", "cwd")]), as.matrix(bioclimParams))

    # Applying AGB formula over simulated matrices and vectors
    AGB_simu <- t(t(log(WD_simu)) * param_7[selec, "logwsg"] +
      t(log(D_simu)) * param_7[selec, "logdbh"] +
      t(log(D_simu)^2) * param_7[selec, "logdbh2"] +
      Esim * -param_7[selec, "E"] +
      param_7[selec, "intercept"])

    if (!is.null(D2)) {
      # DO THESE RANDOM ERRORS NEED TO BE ESTIMATED SEPARATELY AT T1 AND T2? ASSUMING NO FOR NOW
#      AGB2_simu <- t(t(log(WD_simu)) * param_7[selec2, "logwsg"] +
#        t(log(D2_simu)) * param_7[selec2, "logdbh"] +
#        t(log(D2_simu)^2) * param_7[selec2, "logdbh2"] +
#        Esim2 * -param_7[selec2, "E"] +
#        param_7[selec2, "intercept"])
      AGB2_simu <- t(t(log(WD_simu)) * param_7[selec, "logwsg"] +
        t(log(D2_simu)) * param_7[selec, "logdbh"] +
        t(log(D2_simu)^2) * param_7[selec, "logdbh2"] +
        Esim * -param_7[selec, "E"] +
        param_7[selec, "intercept"])
    }

    # Construct a matrix where each column contains random errors taken from N(0,RSEi) with i varying between 1 and n
    matRSE <- mapply(function(y) {
      rnorm(sd = y, n = len)
    }, y = RSE)
    AGB_simu <- AGB_simu + matRSE
    AGB_simu <- exp(AGB_simu) / 1000
    if (!is.null(D2)) {
      # DO THESE RANDOM ERRORS NEED TO BE ESTIMATED SEPARATELY AT T1 AND T2? ASSUMING NO FOR NOW
#      matRSE2 <- mapply(function(y) {
#        rnorm(sd = y, n = len)
#      }, y = RSE2)
#      AGB2_simu <- AGB2_simu + matRSE2
      AGB2_simu <- AGB2_simu + matRSE
      AGB2_simu <- exp(AGB2_simu) / 1000
    }
  }

  if (!is.null(Dlim)) {
    AGB_simu[D < Dlim, ] <- 0
    if (!is.null(D2)) {
      AGB2_simu[D2 < Dlim, ] <- 0
    }
  }
  AGB_simu[ which(is.infinite(AGB_simu)) ] <- NA
  if (!is.null(D2)) {
    AGB2_simu[ which(is.infinite(AGB2_simu)) ] <- NA
  }

  if (Carbon == FALSE) {
    sum_AGB_simu <- colSums(AGB_simu, na.rm = TRUE)
    res <- list(
      meanAGB = mean(sum_AGB_simu),
      medAGB = median(sum_AGB_simu),
      sdAGB = sd(sum_AGB_simu),
      credibilityAGB = quantile(sum_AGB_simu, probs = c(0.025, 0.975)),
      AGB_simu = AGB_simu
    )
    if (!is.null(D2)) {
      sum_AGB2_simu <- colSums(AGB2_simu, na.rm = TRUE)
      diff_AGB_simu <- AGB2_simu - AGB_simu
      sum_Diff_simu <- colSums(diff_AGB_simu, na.rm = TRUE)
      res2 <- list(
        meanAGB2 = mean(sum_AGB2_simu),
        medAGB2 = median(sum_AGB2_simu),
        sdAGB2 = sd(sum_AGB2_simu),
        credibilityAGB2 = quantile(sum_AGB2_simu, probs = c(0.025, 0.975)),
        AGB2_simu = AGB2_simu,
        meanDiff = mean(sum_Diff_simu),
        medDiff = median(sum_Diff_simu),
        sdDiff = sd(sum_Diff_simu),
        credibilityDiff = quantile(sum_Diff_simu, probs = c(0.025, 0.975)),
        diff_AGB_simu = diff_AGB_simu
      )
      res <- c(res, res2)
    }
  } else {
    # Biomass to carbon ratio calculated from Thomas and Martin 2012 forests data stored in DRYAD (tropical
    # angiosperm stems carbon content)
    AGC_simu <- AGB_simu * rnorm(mean = 47.13, sd = 2.06, n = n * len) / 100
    sum_AGC_simu <- colSums(AGC_simu, na.rm = TRUE)
    if (!is.null(D2)) {
      AGC2_simu <- AGB2_simu * rnorm(mean = 47.13, sd = 2.06, n = n * len) / 100
      sum_AGC2_simu <- colSums(AGC2_simu, na.rm = TRUE)
      diff_AGC_simu <- AGC2_simu - AGC_simu
      sum_Diff_simu <- colSums(diff_AGC_simu, na.rm = TRUE)
    }
    res <- list(
      meanAGC = mean(sum_AGC_simu),
      medAGC = median(sum_AGC_simu),
      sdAGC = sd(sum_AGC_simu),
      credibilityAGC = quantile(sum_AGC_simu, probs = c(0.025, 0.975)),
      AGC_simu = AGC_simu
    )
    if (!is.null(D2)) {
      res2 <- list(
        meanAGC2 = mean(sum_AGC2_simu),
        medAGC2 = median(sum_AGC2_simu),
        sdAGC2 = sd(sum_AGC2_simu),
        credibilityAGC2 = quantile(sum_AGC2_simu, probs = c(0.025, 0.975)),
        AGC2_simu = AGC2_simu,
        meanDiff = mean(sum_Diff_simu),
        medDiff = median(sum_Diff_simu),
        sdDiff = sd(sum_Diff_simu),
        credibilityDiff = quantile(sum_Diff_simu, probs = c(0.025, 0.975)),
        diff_AGC_simu = diff_AGC_simu,
      )
      res <- c(res, res2)
    }
  }
  return(res)
}
