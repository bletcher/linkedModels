#'Run the growth model using Nimble
#'
#'@param d, input dataframe
#'@param mcmcInfo, a list containing run info
#'@return a data frame
#'@export

runGrowthModel_Nimble <- function(d,mcmcInfo,code){

  ##
  nBetas <- 8
  nBetasSigma <- 3
  nBetasBNT <- 3
  nBetasATS <- 1

  constants <- list(riverDATA = d$riverDATA,
                    nRivers = d$nRivers,
                    nSpecies = d$nSpecies,
                    ind = d$ind,
                    nInd = d$nInd,
                    season = d$season,
                    nEvalRows = d$nEvalRows,
                    evalRows = d$evalRows,
                    nSeasons = d$nSeasons,
                  #  countPAllSppStd = d$countPAllSppStd,
                    cBKTStd = d$cBKTStd,
                    cBNTStd = d$cBNTStd,
                    cATSStd = d$cATSStd,
                    BKT01DATA = d$BKT01DATA,
                    BNT01DATA = d$BNT01DATA,
                    ATS01DATA = d$ATS01DATA,

                    isYOYDATA = d$isYOYDATA,
                    tempStd = d$tempStd,
                    flowStd = d$flowStd,
                    tempStd2 = d$tempStd2,
                    flowStd2 = d$flowStd2,
                    #lengthDATAStd = d$lengthDATAStd,
                    lengthMean = d$lengthMean,
                    lengthSD = d$lengthSD,
                    sampleInterval = d$sampleInterval[1:(max(d$evalRows)+1)],
                    nBetas = nBetas,
                    nBetasSigma = nBetasSigma,
                    nBetasATS = nBetasATS,
                    nBetasBNT = nBetasBNT
              #      predRange = seq(-15,15,5)
                    )
  ##
  data <- list(lengthDATA = d$lengthDATA[1:(max(constants$evalRows)+1)]) # so don't have trailing single obs fish at end of df
  print(paste("Trimmed", length(d$lengthDATA) - length(data$lengthDATA), "fish that had single observation(s) at end of df"))

  ##
  inits <- list(grInt = array(rnorm(2 * constants$nSeasons * constants$nRivers, 0.5, 0.25), c(2, constants$nSeasons, constants$nRivers)),
                grBeta = array(rnorm(nBetas * 2 * constants$nSeasons * constants$nRivers, 0, 0.1), c(nBetas, 2, constants$nSeasons, constants$nRivers)),
                ##grSigma = array(runif(2 * constants$nSeasons * constants$nRivers, 0, 0.05), c(2, constants$nSeasons, constants$nRivers)),
                sigmaBeta = array(rnorm(nBetasSigma * 2 * constants$nSeasons * constants$nRivers, 0, 0.05), c(nBetasSigma, 2, constants$nSeasons, constants$nRivers)),
                #grATSOffset = array(rnorm( 2 * constants$nSeasons, 0, 0.1), c(2,constants$nSeasons) ),
                grIndRE = rnorm(constants$nInd, 0, 0.1),
                sigmaIndRE = 1,
                grIntMu = rep(0, constants$nSeasons),
                grIntSigma = rep(1, constants$nSeasons),
                grBetaMu = array(0, c(nBetas, constants$nSeasons)),
                grBetaSigma = array(1, c(nBetas, constants$nSeasons)),
                sigmaInt = array(1, c(2, constants$nSeasons, constants$nRivers)),
                sigmaIntMu = rep(0, constants$nSeasons),
                sigmaIntSigma = rep(1, constants$nSeasons),
                sigmaBetaMu = array(0, c(nBetasSigma, constants$nSeasons)),
                sigmaBetaSigma = array(1, c(nBetasSigma, constants$nSeasons)),
                # mu = rnorm(max(constants$evalRows), 0.33, 0.05),
                # gr = rnorm(max(constants$evalRows), 0.33, 0.05),
                # expectedGRSigma = runif(max(constants$evalRows), 0, 1),
                # lengthExp = runif(max(constants$evalRows) + 1, 60, 200),
                # lengthStd = rnorm(max(constants$evalRows), 0, 1)
                #grBetaATS = array(rnorm( nBetasATS * 2 * constants$nSeasons, 0, 0.1), c(nBetasATS, 2, constants$nSeasons)),
                #grBetaATSMu = array(0, c(nBetasATS, constants$nSeasons)),
                #grBetaATSSigma = array(1, c(nBetasATS, constants$nSeasons)),
                grBetaATS = array(rnorm( nBetasATS * 2 * constants$nSeasons * constants$nRivers, 0, 0.1), c(nBetasATS, 2, constants$nSeasons, constants$nRivers)),
                grBetaATSMu = array(0, c(nBetasATS, constants$nSeasons)),
                grBetaATSSigma = array(1, c(nBetasATS, constants$nSeasons)),
                grBetaBNT = array(rnorm( nBetasBNT * 2 * constants$nSeasons * constants$nRivers, 0, 0.1), c(nBetasBNT, 2, constants$nSeasons, constants$nRivers)),
                grBetaBNTMu = array(0, c(nBetasBNT, constants$nSeasons)),
                grBetaBNTSigma = array(1, c(nBetasBNT, constants$nSeasons))

                #sigmaBetaATS = array(rnorm( nBetasATS * 2 * constants$nSeasons, 0, 0.1), c(nBetasATS, 2, constants$nSeasons)),
                #sigmaBetaATSMu = array(0, c(nBetasATS, constants$nSeasons)),
                #sigmaBetaATSSigma = array(1, c(nBetasATS, constants$nSeasons)),
                # sigmaBetaATS = array(rnorm( nBetasATS * 2 * constants$nSeasons * constants$nRivers, 0, 0.1), c(nBetasATS, 2, constants$nSeasons, constants$nRivers)),
                # sigmaBetaATSMu = array(0, c(nBetasATS, constants$nSeasons)),
                # sigmaBetaATSSigma = array(1, c(nBetasATS, constants$nSeasons)),
                # sigmaBetaBNT = array(rnorm( nBetasBNT * 2 * constants$nSeasons * constants$nRivers, 0, 0.1), c(nBetasBNT, 2, constants$nSeasons, constants$nRivers)),
                # sigmaBetaBNTMu = array(0, c(nBetasBNT, constants$nSeasons)),
                # sigmaBetaBNTSigma = array(1, c(nBetasBNT, constants$nSeasons))
  )
  ##
  params <- c('grInt', 'grIntMu', 'grIntSigma', 'grBeta', 'grBetaMu',
              'grBetaSigma', 'sigmaInt', 'sigmaIntMu', 'sigmaIntSigma',
              'sigmaBeta', 'sigmaBetaMu', 'sigmaBetaSigma', 'grIndRE', 'sigmaIndRE',
              'lengthExp',
              'grBetaATS','grBetaATSMu','grBetaATSSigma',
              'grBetaBNT','grBetaBNTMu','grBetaBNTSigma')

  return(list(code=code,data=data,constants=constants,inits=inits,params=params))
}


#'Run the growth model using jags
#'
#'@param d, the input dataframe
#'@param paralel, boolean for running chains in parallel
#'@return a data frame
#'@export

runGrowthModel <- function(d, parallel = FALSE){

  #grBetaOutside = array( runif( 11 *2*d$nSeasons*d$nRivers,-2,2),c( 11 ,2,d$nSeasons,d$nRivers))
  #grBetaOutside[ ,1,,2, ] <- 0

  nBetas <- 12
  nBetasSigma <- 6

  inits <- function(){
    list(
      grInt = array(rnorm(2*d$nSeasons*d$nRivers,0.5,0.25),c(2,d$nSeasons,d$nRivers)),
      #grBeta[x,1:2,1,1:4,1:4]
      grBeta = array(rnorm(nBetas*2*d$nSeasons*d$nRivers,0,0.1),c(nBetas,2,d$nSeasons,d$nRivers)),
      #grSigma[ yoy,spp,s,r ]
      grSigma = array(runif(2*d$nSeasons*d$nRivers,0,0.05),c(2,d$nSeasons,d$nRivers))
      # sigmaBeta[ b,yoy,spp,s,r ]
    #  sigmaBeta = array(rnorm(nBetasSigma*2*d$nSeasons*d$nRivers,0,0.05),c(nBetasSigma,2,d$nSeasons,d$nRivers))
   #   grIndRE = rnorm(d$nInd,0,0.001)
    )
  }

  # params <- c('grInt' ,'grIntMu','grIntSigma'
  #             ,'sigmaInt','sigmaIntMu','sigmaIntSigma'
  #             ,'grBeta','grBetaMu','grBetaSigma'
  #             , 'length','expectedGR'#,'expectedGR'
  #    #          , 'grIndRE','grIndREMean','grIndRETau'
  #             )

  params <- c('grInt','grIntMu','grIntSigma',
              'grBeta','grBetaMu', 'grBetaSigma',
              'sigmaInt','sigmaIntMu','sigmaIntSigma',
              'sigmaBeta','sigmaBetaMu', 'sigmaBetaSigma',
              'grIndRE', 'sigmaIndRE',
              'lengthExp')

  outGR <- jags(data = d,
                inits = inits,
                parameters.to.save = params,
                model.file = "./jags/grModel7.jags",
                n.chains = 3,
                n.adapt = 500, #1000
                n.iter = 2500,
                n.burnin = 1500,
                n.thin = 5,
                parallel = parallel
  )

  #outGR$movementModelIterUsed <- iter
  return(outGR)
}

#'Adjust counts
#'
#'@param cdIn core data
#'@param ddDIn results from the detection model
#'@param ddddGIn detection data
#'@param meanOrIterIn whether the growth model gets mean P from the detection model, or results from an iteration
#'@param sampleToUse if detection results are from an iteration, which iteration
#'@return a data frame with standardized counts, nAllFishBySpeciesPStd is standardized by species counts, nAllFishBySpeciesPAllSppStd is standardized by the sum of counts of all species in the analysis
#'@return Std n's for missing occasion samples (e.g. some winters) are interpolated
#'@export
#'


adjustCounts <- function( cdIn,ddDIn,ddddDIn,meanOrIterIn,sampleToUse ){

  if ( meanOrIterIn == "mean") {
    den <- getDensities(ddddDIn, ddDIn, meanOrIterIn, sampleToUse )
    #ggplot(filter(den,!is.na(nAllFishBySpeciesYOYP)),aes(yearN,nAllFishBySpeciesYOYP, color = factor(speciesN))) + geom_point() + geom_line() + facet_grid(riverN~isYOYN+seasonN)
    #ftable(den$yearN,den$speciesN,is.na(den$nAllFishBySpeciesP))
    #ggplot(filter(den,!is.na(nAllFishBySpeciesP)),aes(count,nAllFishBySpeciesP, color = species)) + geom_point()
    den$meanOrIter <- meanOrIterIn
    den$iterIn <- meanOrIterIn
    print("in mean")
  }

  if ( meanOrIterIn == "iter") {
    den <- getDensities(ddddDIn, ddDIn, meanOrIterIn, sampleToUse )
    den$meanOrIter <- meanOrIterIn
    den$iterIn <- sampleToUse

    print("in iter ")
    print(c(meanOrIterIn,sampleToUse))
  }

  denForMerge <- den %>%
    dplyr::select(isYOYN,species, season, riverOrdered, year,
                  nAllFishBySpeciesYOYP,
                  nAllFishBySpeciesYOY
                 # nAllFishP,
                  # massAllFishBySpeciesYOY, massAllFish,
                  #meanOrIter, iterIn
                  ) %>%
  #  filter( !is.na(nAllFishBySpeciesP) ) %>%
    arrange(species,riverOrdered,year,season) %>%
    distinct(isYOYN,species,riverOrdered,year,season,nAllFishBySpeciesYOYP,nAllFishBySpeciesYOY)

  #ftable(denForMerge$isYOYN,denForMerge$species,denForMerge$riverOrdered,denForMerge$season,denForMerge$year)

#######

  # stats for counts by yoy, species, average over years
  denForMergeSummaryBySpeciesYOY <- denForMerge %>%
    group_by( isYOYN,species,season,riverOrdered ) %>%
    summarize( nAllFishBySpeciesYOYPMean = mean( nAllFishBySpeciesYOYP, na.rm = T ),
               nAllFishBySpeciesYOYPSD =     sd( nAllFishBySpeciesYOYP, na.rm = T)
               #massAllFishBySpeciesMean = mean( massAllFishBySpecies, na.rm = T ),
               #massAllFishBySpeciesSD =     sd( massAllFishBySpecies, na.rm = T)
    )

  denForMerge2 <- left_join( denForMerge, denForMergeSummaryBySpeciesYOY ) %>%
    mutate( nAllFishBySpeciesYOYPStd = ( nAllFishBySpeciesYOYP - nAllFishBySpeciesYOYPMean ) / nAllFishBySpeciesYOYPSD )

  # create template for all possible occasions
  possibleOccasions <- data.frame(table(denForMerge$season,denForMerge$year,denForMerge$species,denForMerge$riverOrdered,denForMerge$isYOYN)) %>%
    mutate( FreqUp = lead(Freq), FreqDown = lag(Freq),
            occFilled = ifelse(FreqUp * FreqDown + Freq > 0,1,0) ) %>%
    filter(occFilled == 1) %>%
    rename( season=Var1, year=Var2, species=Var3, riverOrdered=Var4, isYOYN = Var5 ) %>%
    mutate( season = as.numeric(season), year = as.numeric(year) + min(denForMerge2$year, na.rm = TRUE) - 1, species = as.character(species),
            riverOrdered = factor(riverOrdered,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered = T),
            isYOYN = as.numeric(isYOYN)
          )

  # use this as merge template, then interpolate std values

  denForMerge_possibleOccasions <- left_join( possibleOccasions,denForMerge2, by = c("isYOYN","species","season", "year", "species", "riverOrdered") ) %>%
    arrange( isYOYN,species,riverOrdered,year,season) %>%
    mutate( nAllFishBySpeciesYOYPStd = zoo::na.approx(nAllFishBySpeciesYOYPStd),
            nAllFishBySpeciesYOYP = zoo::na.approx(nAllFishBySpeciesYOYP)) %>%
    dplyr::select(-c(Freq,FreqUp,FreqDown))

### isYOYN==1 and season == 2 is getting interpolated now - set values where no fish were caught to NA
  denForMerge_possibleOccasions$nAllFishBySpeciesYOYP <- ifelse(
    denForMerge_possibleOccasions$season == 2 &
    denForMerge_possibleOccasions$isYOYN == 1 &
    is.na(denForMerge_possibleOccasions$nAllFishBySpeciesYOY), NA, denForMerge_possibleOccasions$nAllFishBySpeciesYOYP )

  denForMerge_possibleOccasions$nAllFishBySpeciesYOYStd <- ifelse(
    denForMerge_possibleOccasions$season == 2 &
    denForMerge_possibleOccasions$isYOYN == 1 &
    is.na(denForMerge_possibleOccasions$nAllFishBySpeciesYOY), NA, denForMerge_possibleOccasions$nAllFishBySpeciesYOYPStd )

  cdIn <- left_join( cdIn,denForMerge_possibleOccasions, by = c("isYOY"="isYOYN","species", "year", "season", "riverOrdered") )

   # get counts of each species by averaging over nAllFishBySpeciesYOYP for each species
   # averaging over isYOY
  denForMerge_possibleOccasions2 <- denForMerge_possibleOccasions %>%
     group_by(species,season,riverOrdered,year) %>%
     summarize( nAllFishBySpeciesPStd_Mean = mean(nAllFishBySpeciesYOYPStd, na.rm=T),
                nAllFishBySpeciesPStd_SD = sd(nAllFishBySpeciesYOYPStd, na.rm=T))

  cdIn <- left_join(cdIn,denForMerge_possibleOccasions2)

   # Averaging over species
  denForMerge_possibleOccasions3 <- denForMerge_possibleOccasions %>%
    group_by(season,riverOrdered,year) %>%
    summarize( nAllFishPStd_Mean = mean(nAllFishBySpeciesYOYPStd, na.rm=T),
               nAllFishPStd_SD = sd(nAllFishBySpeciesYOYPStd, na.rm=T))

  cdIn <- left_join(cdIn,denForMerge_possibleOccasions3)

  #####
  # counts by species in separate columns
  nAllFishBySpeciesPStdBySpp <- denForMerge_possibleOccasions2 %>%
    dplyr::select(species, season, riverOrdered, year, nAllFishBySpeciesPStd_Mean) %>%
    distinct(species, season, riverOrdered, year, nAllFishBySpeciesPStd_Mean) %>%
    spread(key=species, value=nAllFishBySpeciesPStd_Mean, fill = -2.5) %>%
    rename(nAllFishBySpeciesPStdBKT = bkt,
           nAllFishBySpeciesPStdBNT = bnt,
           nAllFishBySpeciesPStdATS = ats
           )

  cdIn <- left_join( cdIn,nAllFishBySpeciesPStdBySpp )

 #  #####
 #  # counts by species and yoy in separate columns
 #  # yoy1 and yoy2 are too highly correlated - don't use
 #  nAllFishBySpeciesPStdBySppYOY <- denForMerge_possibleOccasions %>%
 #    dplyr::select(isYOYN,species, season, riverOrdered, year, nAllFishBySpeciesPStd) %>%
 #    unite(yoySpp, species, isYOYN, sep = "_") %>%
 #    spread(key=yoySpp, value=nAllFishBySpeciesPStd, fill = -2.5) %>%
 #    rename(nAllFishBySpeciesPStdBKT_yoy1 = bkt_1,
 #           nAllFishBySpeciesPStdBKT_yoy2 = bkt_2,
 #           nAllFishBySpeciesPStdBNT_yoy1 = bnt_1,
 #           nAllFishBySpeciesPStdBNT_yoy2 = bnt_2,
 #           nAllFishBySpeciesPStdATS_yoy1 = ats_1,
 #           nAllFishBySpeciesPStdATS_yoy2 = ats_2)
 #
 #  cdIn <- left_join( cdIn,nAllFishBySpeciesPStdBySppYOY )
 # # ggplot(nAllFishBySpeciesPStdBySppYOY , aes(nAllFishBySpeciesPStdBKT_yoy1,nAllFishBySpeciesPStdBKT_yoy2)) + geom_point()
 # # ggplot(nAllFishBySpeciesPStdBySppYOY , aes(nAllFishBySpeciesPStdBNT_yoy1,nAllFishBySpeciesPStdBNT_yoy2)) + geom_point()
 # #  ggplot(nAllFishBySpeciesPStdBySppYOY , aes(nAllFishBySpeciesPStdATS_yoy1,nAllFishBySpeciesPStdATS_yoy2)) + geom_point()
 #  nAllFishBySpeciesPStdBySppYOY[nAllFishBySpeciesPStdBySppYOY == -2.5] <- NA
 #  round(cor(nAllFishBySpeciesPStdBySppYOY[,c('nAllFishBySpeciesPStdBKT_yoy1','nAllFishBySpeciesPStdBKT_yoy2','nAllFishBySpeciesPStdBNT_yoy1','nAllFishBySpeciesPStdBNT_yoy2','nAllFishBySpeciesPStdATS_yoy1','nAllFishBySpeciesPStdATS_yoy2')], use="complete.obs"),2)

  return(cdIn)

}


#'Turn observedLength values to NA for a percentage of the observations
#'
#'@param d a dataframe
#'@param runCrossValidation boolean for running cross validation or not
#'@return a data frame with observedLength set to NA for percentLeftOut observations
#'@export
#'
crossValidate <- function(d, runCrossValidation_TF){
  d$observedLengthOriginal <- d$observedLength
  if ( runCrossValidation_TF ) {
    propFOcc <- sum(d$fOcc) / nrow(d)
    d$leftOut <- ifelse( (runif(nrow(d)) < percentLeftOut/propFOcc/100) & (d$fOcc == 0), T, F ) # adjust percentLeftOut higher to acct for the fOcc that can't be NA
    d$observedLength <- ifelse( d$leftOut, NA, d$observedLength )
  }
  return(d)
}





#'Adjust counts for all fish
#'
#'@param nAll counts of all fish
#'@param ddDIn results from the detection model
#'@param ddddGIn detection data
#'@param meanOrIterIn whether the growth model gets mean P from the detection model, or results from an iteration
#'@param sampleToUse if detection results are from an iteration, which iteration
#'@return a data frame with standardized counts, nAllFishBySpeciesPStd is standardized by species counts, nAllFishBySpeciesPAllSppStd is standardized by the sum of counts of all species in the analysis
#'@return Std n's for missing occasion samples (e.g. some winters) are interpolated
#'@export
#'


adjustCounts_allFish <- function( nAll,ddDIn,ddddDIn,meanOrIterIn,sampleToUse ){

  if ( meanOrIterIn == "mean") {
    den <- getDensities_allFish( nAll, ddDIn, meanOrIterIn, sampleToUse )
    #ggplot(filter(den,!is.na(nAllFishBySpeciesYOYP)),aes(yearN,nAllFishBySpeciesYOYP, color = factor(speciesN))) + geom_point() + geom_line() + facet_grid(riverN~isYOYN+seasonN)
    #ftable(den$yearN,den$speciesN,is.na(den$nAllFishBySpeciesP))
    #ggplot(filter(den,!is.na(nAllFishBySpeciesP)),aes(count,nAllFishBySpeciesP, color = species)) + geom_point()
    den$meanOrIter <- meanOrIterIn
    den$iterIn <- meanOrIterIn
    print("in mean")
  }

  if ( meanOrIterIn == "iter") {
    den <- getDensities(ddddDIn, ddDIn, meanOrIterIn, sampleToUse )
    den$meanOrIter <- meanOrIterIn
    den$iterIn <- sampleToUse

    print("in iter ")
    print(c(meanOrIterIn,sampleToUse))
  }

  denForMerge <- den %>%
    dplyr::select(isYOYN,species, season, riverOrdered, year,
                  nAllFishBySpeciesYOYP,
                  nAllFishBySpeciesYOY
                  # nAllFishP,
                  # massAllFishBySpeciesYOY, massAllFish,
                  #meanOrIter, iterIn
    ) %>%
    #  filter( !is.na(nAllFishBySpeciesP) ) %>%
    arrange(species,riverOrdered,year,season) %>%
    distinct(isYOYN,species,riverOrdered,year,season,nAllFishBySpeciesYOYP,nAllFishBySpeciesYOY)

  #ftable(denForMerge$isYOYN,denForMerge$species,denForMerge$riverOrdered,denForMerge$season,denForMerge$year)

  #######

  # stats for counts by yoy, species, average over years
  denForMergeSummaryBySpeciesYOY <- denForMerge %>%
    group_by( isYOYN,species,season,riverOrdered ) %>%
    summarize( nAllFishBySpeciesYOYPMean = mean( nAllFishBySpeciesYOYP, na.rm = T ),
               nAllFishBySpeciesYOYPSD =     sd( nAllFishBySpeciesYOYP, na.rm = T)
               #massAllFishBySpeciesMean = mean( massAllFishBySpecies, na.rm = T ),
               #massAllFishBySpeciesSD =     sd( massAllFishBySpecies, na.rm = T)
    )

  denForMerge2 <- left_join( denForMerge, denForMergeSummaryBySpeciesYOY ) %>%
    mutate( nAllFishBySpeciesYOYPStd = ( nAllFishBySpeciesYOYP - nAllFishBySpeciesYOYPMean ) / nAllFishBySpeciesYOYPSD )

  # create template for all possible occasions
  possibleOccasions <- data.frame(table(denForMerge$season,denForMerge$year,denForMerge$species,denForMerge$riverOrdered,denForMerge$isYOYN)) %>%
    mutate( FreqUp = lead(Freq), FreqDown = lag(Freq),
            occFilled = ifelse(FreqUp * FreqDown + Freq > 0,1,0) ) %>%
    filter(occFilled == 1) %>%
    rename( season=Var1, year=Var2, species=Var3, riverOrdered=Var4, isYOYN = Var5 ) %>%
    mutate( season = as.numeric(season), year = as.numeric(year) + min(denForMerge2$year, na.rm = TRUE) - 1, species = as.character(species),
            riverOrdered = factor(riverOrdered,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels = c("west brook","wb jimmy","wb mitchell","wb obear"), ordered = T),
            isYOYN = as.numeric(isYOYN)
    )

  # use this as merge template, then interpolate std values

  denForMerge_possibleOccasions <- left_join( possibleOccasions,denForMerge2, by = c("isYOYN","species","season", "year", "species", "riverOrdered") ) %>%
    arrange( isYOYN,species,riverOrdered,year,season) %>%
    mutate( nAllFishBySpeciesYOYPStd = zoo::na.approx(nAllFishBySpeciesYOYPStd),
            nAllFishBySpeciesYOYP = zoo::na.approx(nAllFishBySpeciesYOYP)) %>%
    dplyr::select(-c(Freq,FreqUp,FreqDown))

  ### isYOYN==1 and season == 2 is getting interpolated now - set values where no fish were caught to NA
  denForMerge_possibleOccasions$nAllFishBySpeciesYOYP <- ifelse(
    denForMerge_possibleOccasions$season == 2 &
      denForMerge_possibleOccasions$isYOYN == 1 &
      is.na(denForMerge_possibleOccasions$nAllFishBySpeciesYOY), NA, denForMerge_possibleOccasions$nAllFishBySpeciesYOYP )

  denForMerge_possibleOccasions$nAllFishBySpeciesYOYStd <- ifelse(
    denForMerge_possibleOccasions$season == 2 &
      denForMerge_possibleOccasions$isYOYN == 1 &
      is.na(denForMerge_possibleOccasions$nAllFishBySpeciesYOY), NA, denForMerge_possibleOccasions$nAllFishBySpeciesYOYPStd )

  nAll <- left_join( nAll,denForMerge_possibleOccasions, by = c("isYOY"="isYOYN","species", "year", "season", "riverOrdered","nAllFishBySpeciesYOY") )

  # get counts of each species by averaging over nAllFishBySpeciesYOYP for each species
  # averaging over isYOY
  denForMerge_possibleOccasions2 <- denForMerge_possibleOccasions %>%
    group_by(species,season,riverOrdered,year) %>%
    summarize( nAllFishBySpeciesPStd_Mean = mean(nAllFishBySpeciesYOYPStd, na.rm=T),
               nAllFishBySpeciesPStd_SD = sd(nAllFishBySpeciesYOYPStd, na.rm=T),
               nAllFishBySpeciesP_Mean = mean(nAllFishBySpeciesYOYP, na.rm=T),
               nAllFishBySpeciesP_SD = sd(nAllFishBySpeciesYOYP, na.rm=T))

  nAll <- left_join(nAll,denForMerge_possibleOccasions2)

  # Averaging over species
  denForMerge_possibleOccasions3 <- denForMerge_possibleOccasions %>%
    group_by(season,riverOrdered,year) %>%
    summarize( nAllFishPStd_Mean = mean(nAllFishBySpeciesYOYPStd, na.rm=T),
               nAllFishPStd_SD = sd(nAllFishBySpeciesYOYPStd, na.rm=T),
               nAllFishP_Mean = mean(nAllFishBySpeciesYOYP, na.rm=T),
               nAllFishP_SD = sd(nAllFishBySpeciesYOYP, na.rm=T))

  nAll <- left_join(nAll,denForMerge_possibleOccasions3)

  #####
  # counts by species in separate columns
  nAllFishBySpeciesPStdBySpp <- denForMerge_possibleOccasions2 %>%
    dplyr::select(species, season, riverOrdered, year, nAllFishBySpeciesPStd_Mean) %>%
    distinct(species, season, riverOrdered, year, nAllFishBySpeciesPStd_Mean) %>%
    spread(key=species, value=nAllFishBySpeciesPStd_Mean, fill = -2.5) %>%
    rename(nAllFishBySpeciesPStdBKT = bkt,
           nAllFishBySpeciesPStdBNT = bnt,
           nAllFishBySpeciesPStdATS = ats
    )

  nAll <- left_join( nAll,nAllFishBySpeciesPStdBySpp )

  #####
  # unstandardized counts by species in separate columns
  # maybe better than above because fills in 0s for no data rather than arbitrary -2.5

  # 05/03/2018 took out fill because don't want value entered when no fish possible, e.g. bnt in obear
  nAllFishBySpeciesPBySpp <- denForMerge_possibleOccasions2 %>%
    dplyr::select(species, season, riverOrdered, year, nAllFishBySpeciesP_Mean) %>%
    distinct(species, season, riverOrdered, year, nAllFishBySpeciesP_Mean) %>%
    spread(key=species, value=nAllFishBySpeciesP_Mean) %>%
    rename(nAllFishBySpeciesPBKT = bkt,
           nAllFishBySpeciesPBNT = bnt,
           nAllFishBySpeciesPATS = ats
    )

  nAll <- left_join( nAll,nAllFishBySpeciesPBySpp )


  return(nAll)
}





#'Remove fish with many intermediate NAs - these cause problems with indRE estimates
#'
#'@param d a dataframe
#'@return a data frame with problem fish removed
#'@export
#'

removeFishWithManyIntermediateNAs <- function(d){

  tagsToRemove <- c('00088cc02a','1c2d6c51d3','00088cc364','1c2c582218')
  d <- d %>% filter( !(tag %in% tagsToRemove) )

  return(d)
}
