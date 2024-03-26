library(dplyr)
library(stringr)


##### READ IN 1 ha PLOTS
#dbhTree <- read.csv("/home/sch/data/devine/external/field-data/dendrometer-test/tree-dbh.csv")

#dbhTree1 <- dbhTree[dbhTree$sample == 1, 'dbh']
#dbhTree2 <- dbhTree[dbhTree$sample == 2, 'dbh']

##### GENERATE 1 ha PLOTS
peruTrees <- read.csv("/home/sch/data/devine/external/field-data/tree-data-BIOMASS.csv")
#peruWD <- getWoodDensity(peruTrees$Genus, peruTrees$Species, stand = peruTrees$SubplotID_unique)
latPeru <- peruTrees$lat
lonPeru <- peruTrees$lon
coordPeru <- cbind(lonPeru, latPeru)

# so we have measured 2.08 ha in total in 2022-2023
# biomass of this area appears to be 511 to 610 Mg Biomass

subplotSize_ha <- 20*20/10000
subplots <- unique(peruTrees$SubplotID_unique)
nsubplots <- length(subplots)
allSubplotsSize <- subplotSize_ha * nsubplots
targetPlotSize_ha <- 1.0
nSubplots <- targetPlotSize_ha / subplotSize_ha

dendroCost <- 25
plots <- 250
years <- 10
carbonPrice <- 10


# create two "plots" from the existing subplots so we get exactly 1.0 ha measured for each
peruTrees %>% count(SubplotID_unique)
fullPlotList <- c("CON3-1-1", "CON3-1-2", "CON3-1-3", "CON3-2-1", "CON3-2-2", "CON3-2-3", 
                  "CON3-3-1", "CON3-3-2", "CON3-3-3", "CON4-1-1", "CON4-1-2", "CON4-1-3", 
                  "CON4-2-1", "CON4-2-2", "CON4-2-3", "CON4-3-1", "CON4-3-2", "CON4-3-3", 
                  "MA1-1-1", "MA1-1-2", "MA1-1-3", "MA1-2-1", "MA1-2-2", "MA1-2-3", 
                  "MA1-3-1", "MA1-3-2", "MA1-3-3", "MA2-1-1", "MA2-1-2", "MA2-1-3", 
                  "MA2-2-1", "MA2-2-2", "MA2-2-3",
                  "MA2-3-1", "MA2-3-2", "MA2-3-3", "CA1-2-1", "CA1-2-2", "CA1-2-3",
#                  "CA1-3-1", "CA1-3-2", "CA1-3-3", "CA1-1-3", "FM1-1-1", "FM1-1-2",
                  "CA1-3-1", "CA1-3-2", "CA1-3-3", "FM1-1-1", "FM1-1-2",
                  "FM1-1-3", "FM1-2-1", "FM1-2-2", "FM1-2-3", "FM1-3-1", "FM1-3-2",
                  "FM1-3-3")

generatePlots <- function(fullPlotList, n) {
    plotList1 <- sample(fullPlotList, n, replace=TRUE)
    peruTrees1 <- subset(peruTrees, SubplotID_unique %in% plotList1)
    plotList2 <- sample(fullPlotList, n, replace=TRUE)
    peruTrees2 <- subset(peruTrees, SubplotID_unique %in% plotList2)
    return(list(peruTrees1, peruTrees2, plotList1, plotList2))
}


####################################
##### FUNCTIONS ######
myrtruncnorm <- function(n, lower = -1, upper = 1, mean = 0, sd = 1) {
  qnorm(runif(n, pnorm(lower, mean = mean, sd = sd), pnorm(upper, mean = mean, sd = sd)), mean = mean, sd = sd)
}

chaveError <- function(x, len, Dthresh) {
  fivePercent <- round(len * 5 / 100)
  largeErrSample <- sample(len, fivePercent)
  D_sd <- 0.0062 * x + 0.0904 # Assigning small errors on the remaining 95% trees
  D_sd[largeErrSample] <- 4.64
  D_adj <- myrtruncnorm(n = len, mean = x, sd = D_sd, lower = 0.1, upper = 500)

  D_err <- D_adj - x

  signArray <- D_err
  signArray[signArray > 0] <- 1
  signArray[signArray <= 0] <- -1
  D_err_Den <- abs(D_err) 
  D_err_Den[x >= Dthresh] <- D_err_Den[x >= Dthresh] - 0.031 
  D_err_Den[D_err_Den < 0] <- 0
  D_err_Den[signArray <= 0] <- -1.0 * D_err_Den[signArray <= 0] 
  D_adj_Den <- x + D_err_Den

  ret <- list(D_adj, D_adj_Den)
  return(ret)
}

#kg
biomass <- function(E, row, dbh) {
  dbh2 <- dbh^2 
  b <- exp(-1.803 - 0.976 * E + 0.976 * log(row) + 2.673 * log(dbh) - 0.0299 * log(dbh2))
  return(b)
}

E <- -0.103 # average of the six values from the sites
row <- 0.569 # wood specific gravity in g / cm3 from Chave data in Peru
metaN <- 1000


#D <- dbhTree2
################# ADD ANOTHER LOOP OR TWO HERE ################################
# ONE LOOP OVER Dthresh AND ANOTHER OVER RANDOM GENERATION OF PLOTS
# THE OUTER LOOP SHOULD BE RANDOM PLOTS AND THE INNER SHOULD BE THRESH
dThreshList <- c(10.0, 20.0, 30.0, 50.0, 80.0, 100.0)

simulationResults <- data.frame(matrix(ncol = 11, nrow = 0))
nReps <- 10
colnames(simulationResults) <- c('rep', 'dThresh', 'sdStandardGrowth', 'sdDendroGrowth',  
                                 'standardAnnualGrowth10th', 'dendroAnnualGrowth10th',
                                 'dendroBenefit', 'countDBH', 'dendroCost_dbhThresh',
                                 'dendroBenefit_perHa', 'applicationAreaRequired')

for (dThresh in dThreshList) {
  for (iRep in 1:nReps) {
    dbhTreesRand <- generatePlots(fullPlotList, nSubplots)
    D <- dbhTreesRand[[1]]$dbh_T1
    len <- length(D)

    metaListDiff <- c()
    metaListStandard <- c()
    metaListDendrometer <- c()

    for (metaI in 1:metaN) {
      n <- 1000
      D_simu_both <- replicate(n, chaveError(D, len, dThresh))
      D_simu <- D_simu_both[1,]
      D_simu_Den <- D_simu_both[2,]

      carbonList <- c()
      carbonDenList <- c()
      for (i in 1:n) {
        biom <- biomass(E, row, D_simu[[i]])
        carbonList <- c(carbonList, sum(biom)/1000 * 0.47) 
        biom_Den <- biomass(E, row, D_simu_Den[[i]])
        carbonDenList <- c(carbonDenList, sum(biom_Den)/1000 * 0.47)
      } 

      biomTruth <- biomass(E, row, D)
  
      sdStandard <- sd(carbonList)
      sdDen <- sd(carbonDenList)
      sdDiff <- sdStandard - sdDen
      metaListDiff <- c(metaListDiff, sdDiff)
      metaListStandard <- c(metaListStandard, sdStandard)
      metaListDendrometer <- c(metaListDendrometer, sdDen)
    }

    sdStandardMeta <- mean(metaListStandard)
    sdDendroMeta <- mean(metaListDendrometer)
  
    # Now consider when we have two measurements in time and
    #   we want to estimate the difference

    meanAnnualGrowth <- 5.5 #tC per ha per year
    sdStandardGrowth <- sqrt(sdStandardMeta^2 + sdStandardMeta^2) 
    sdDendroGrowth <- sqrt(sdDendroMeta^2 + sdDendroMeta^2) 

    tVal_10thpct <- qt(1.0 - 0.90, 1000)

    standardAnnualGrowth10th <- meanAnnualGrowth + tVal_10thpct * sdStandardGrowth
    dendroAnnualGrowth10th <- meanAnnualGrowth + tVal_10thpct * sdDendroGrowth
    dendroBenefit <- dendroAnnualGrowth10th - standardAnnualGrowth10th

    countDBH <- length(D[D > dThresh])

    dendroCost_dbhThresh <- dendroCost * countDBH * plots
    dendroBenefit_perHa <- dendroBenefit * years * carbonPrice # 10 years and 10 dollars
    applicationAreaRequired <- dendroCost_dbhThresh / dendroBenefit_perHa

    simulationResults[nrow(simulationResults) + 1,] = c(iRep, dThresh, sdStandardGrowth,
                                 sdDendroGrowth, standardAnnualGrowth10th, dendroAnnualGrowth10th,
                                 dendroBenefit, countDBH, dendroCost_dbhThresh,
                                 dendroBenefit_perHa, applicationAreaRequired)

    print(c(iRep, dThresh, sdStandardGrowth, sdDendroGrowth,  
            standardAnnualGrowth10th, dendroAnnualGrowth10th,
            dendroBenefit, countDBH, dendroCost_dbhThresh,
            dendroBenefit_perHa, applicationAreaRequired))

  }
}

write.csv(simulationResults, "/home/sch/data/devine/external/field-data/dendrometerCostBenefit.csv", row.names=FALSE)

plotColors <- c('black', 'blue', 'red', 'green', 'cyan')
plotNum <- 1
benefitList <- c()
treeList <- c()
costList <- c()
areaList <- c()
#for (iStudy in uniqueStudyList) {
for (iDThresh in dThreshList) {
  tempSub <- subset(simulationResults, dThresh == iDThresh) 
  benefitList <- c(benefitList, mean(tempSub$dendroBenefit))
  treeList <- c(treeList, mean(tempSub$countDBH))
  costList <- c(costList, mean(tempSub$dendroCost_dbhThresh))
  areaList <- c(areaList, mean(na.omit(tempSub$applicationAreaRequired)))
}

plot(dThreshList[2:6], areaList[2:6], xlab = 'dendrometer use threshold (dbh in cm)', ylab = 'project area required to break even (ha)')
plot(dThreshList[2:6], benefitList[2:6], xlab = 'dendrometer use threshold (dbh in cm)', ylab = 'additional carbon benefit (t C per ha per year)')
plot(dThreshList[2:6], treeList[2:6],  xlab = 'dendrometer use threshold (dbh in cm)', ylab = 'approximate number of trees requiring dendrometer (per ha)')
plot(dThreshList[2:6], costList[2:6],  xlab = 'dendrometer use threshold (dbh in cm)', ylab = 'total cost for purchase of dendrometers (USD)') 
