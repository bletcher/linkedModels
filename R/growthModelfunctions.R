#'Get data ready to run the growth model
#'
#'@param d a data frame created using getCoreData()
#'@return A list for importing into the jags run. Run with runGrowthModel()
#'@export

prepareDataForJags <- function(d){

  nInd <- n_distinct(d$tag, na.rm = T)
  nOcc <- n_distinct(d$sampleIndex, na.rm = T)
  minOcc <- min(d$sampleIndex)
  nRivers <- n_distinct(d$riverOrdered, na.rm = T)
  nSeasons <- n_distinct(d$season, na.rm = T)
  nYears <- n_distinct(d$year, na.rm = T)

  d <- d %>%
    group_by(tag) %>%
    mutate( minOcc = min(sampleNumber), fOcc = (sampleNumber == minOcc)*1,
            maxOcc = max(sampleNumber), lOcc = (sampleNumber == maxOcc)*1
    )


  evalRows <- which( d$lOcc == 0 )
  firstObsRows <- which( d$fOcc == 1 )
  lastObsRows <- which( d$lOcc == 1 )

  nEvalRows <- length(evalRows)
  nFirstObsRows <- length(firstObsRows)
  nLastObsRows <- length(lastObsRows)


  load(file = "./data/cutoffYOYInclSpring1DATA.RData")
  cutoffYOYDATA <- cutoffYOYInclSpring1DATA # will need to add years after 2012, I think

  d$riverN <- as.numeric(d$riverOrdered)

  means <- d %>%
    group_by(season,riverN) %>%
    summarize( meanLen = mean(observedLength, na.rm = T),
               sdLen = sd(observedLength, na.rm = T))

  data <- list( lengthDATA = d$observedLength,
                riverDATA = d$riverN,
                #ind = d$tagIndex,
                nRivers = nRivers,
                #nInd = nInd,
                #nOcc = nOcc,
                #occ = d$sampleIndex - minOcc + 1,
                season = d$season,
                year = d$year - min(d$year) + 1,
                nYears = nYears, #max(d$year) - min(d$year) + 1,
                nEvalRows = nEvalRows, evalRows = evalRows,
                nFirstObsRows = nFirstObsRows, firstObsRows = firstObsRows,
                nLastObsRows = nLastObsRows, lastObsRows = lastObsRows,
                #nSeasons = nSeasons,
                lengthMean = matrix(means$meanLen,c(nSeasons,nRivers),byrow = T),
                lengthSD = matrix(means$sdLen,c(nSeasons,nRivers),byrow = T),
                cutoffYOYDATA = cutoffYOYDATA
  )
  return(data)
}



#'Run the growth model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runGrowthModel <- function(d){   #iterToUse, firstNonBurnIter, chainToUse, simInfo, grDir, out){

  inits <- function(){
    #list(grBetaInt = array(rnorm(d$nSeasons*d$nRivers,0,2.25),c(d$nSeasons,d$nRivers)))
    list(grBetaInt = array(rnorm(2*4*4*6,0,2.25),c(2,4,4,6)))
  }

  params <- c("grBetaInt","muGrBetaInt","sigmaGrBetaInt","grBeta","grSigmaBeta","length")

  outGR <- jags(data = d,
                inits = inits,
                parameters.to.save = params,
                model.file = "./jags/grModel.jags",
                n.chains = 3,
                n.adapt = 1000, #1000
                n.iter = 2000,
                n.burnin = 1000,
                n.thin = 3,
                parallel = FALSE
  )

  #outGR$movementModelIterUsed <- iter
  return(outGR)
}

