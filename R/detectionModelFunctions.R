#'Run the detection model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runDetectionModel <- function(d, parallel = FALSE){   #iterToUse, firstNonBurnIter, chainToUse, simInfo, grDir, out){

  inits <- function(){
    list(#pBetaInt = array(rnorm(d$nSeasons*d$nRivers*d$nYears,0,2.25),c(d$nSeasons,d$nRivers,d$nYears)),
         pBetaInt = array(runif(d$nSpecies*d$nSeasons*d$nRivers*d$nYears, -2.5, 2),c(d$nSpecies,d$nSeasons,d$nRivers,d$nYears))
#        z = d$zForInit
         )
  }

  params <- c("pBetaInt","pBetaIntMean","pBetaIntSigma","phiBetaInt")

  outDet <- jags(data = d,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "./jags/detModel.jags",
                 n.chains = 3,
                 n.adapt = 50, #1000
                 n.iter = 200,
                 n.burnin = 50,
                 n.thin = 4,
                 parallel = parallel
                )

  #outDet$movementModelIterUsed <- iter
  return(outDet)
}


#'Get (estimate) densities by section, river, year based on num sammples and p from the density model
#'
#'@param dddd the input data frame for the detection model run
#'@param dd the output data frame for the detection model run
#'@param meanOrIter Whether model run output is the mean ('mean') of all iterations or a single iteration ('iter'). Default is 'mean'.
#'@param sampleToUse if meanOrIter == 'iter', which iteration to use. Assumes use of chain 1.
#'@param chainToUse if meanOrIter == 'iter', which chain to use. Default is chain 1.
#'@return a data frame
#'@export

getDensities <- function(dddd,dd, meanOrIter = "mean", sampleToUse = sampleToUse){

  try( if (sampleToUse > dd$mcmc.info$n.samples) stop("requested iteration beyond max # of iterations") )

  counts <- dddd %>%
    filter( enc == 1 ) %>%
    group_by( species,season,riverOrdered,year ) %>%
    dplyr::summarise( count = n() ) %>%
    mutate( speciesN = as.numeric(as.factor(species)),
            seasonN = season,
            riverN = as.numeric(riverOrdered),
            yearN = year - min(year) + 1
          ) %>%
    ungroup()

  if ( meanOrIter == 'mean') ddIn <- dd$q50$pBetaInt
  if ( meanOrIter == 'iter') ddIn <- dd$sims.list$pBetaInt[ sampleToUse,,,, ]

  #print(c("in getDensities()",meanOrIter,ddIn))

  #pBetaInt[ species,season,riverDATA,year ]
  # convert array to data frame
  p <- as.data.frame.table(ddIn) %>%
    rename( speciesN = Var1,
            seasonN = Var2,
            riverN = Var3,
            yearN = Var4,
            logitP = Freq
          ) %>%
    mutate( speciesN = as.numeric(speciesN),
            seasonN = as.numeric(seasonN),
            riverN = as.numeric(riverN),
            yearN = as.numeric(yearN),
            p = invlogit( logitP )
          )

  p <- left_join( p, counts ) %>%
    mutate( countP = count/p )

  return(p)
}


#'Add survival estimates from detection model by section, river, year b
#'
#'@param dddd the input data frame for the detection model run
#'@param dd the output data frame for the detection model run
#'@param meanOrIter Whether model run output is the mean ('mean') of all iterations or a single iteration ('iter'). Default is 'mean'.
#'@param sampleToUse if meanOrIter == 'iter', which iteration to use. Assumes use of chain 1.
#'@param chainToUse if meanOrIter == 'iter', which chain to use. Default is chain 1.
#'@return a data frame with merged phi's
#'@export

addSurvivals <- function(dddd,dd, meanOrIter = "mean", sampleToUse = sampleToUse){

  try( if (sampleToUse > dd$mcmc.info$n.samples) stop("requested iteration beyond max # of iterations") )

  if ( meanOrIter == 'mean') ddIn <- dd$q50$phiBetaInt
  if ( meanOrIter == 'iter') ddIn <- dd$sims.list$phiBetaInt[ sampleToUse,,,, ]

  #print(c("in getDensities()",meanOrIter,ddIn))

  #phiBetaInt[ species,season,riverDATA,year ]
  # convert array to data frame
  phi <- as.data.frame.table(ddIn) %>%
    rename( speciesN = Var1,
            seasonN = Var2,
            riverN = Var3,
            yearN = Var4,
            logitPhi = Freq
    ) %>%
    mutate( speciesN = as.numeric(speciesN),
            seasonN = as.numeric(seasonN),
            riverN = as.numeric(riverN),
            yearN = as.numeric(yearN),
            phi = invlogit( logitPhi ),
            logitPhiStd = ( logitPhi - mean(logitPhi,na.rm=T) ) / sd(logitPhi,na.rm=T)
    )



  # add temporary variables for merging
  dddd <- dddd %>%
    mutate( speciesN = as.numeric(as.factor(species)),
            seasonN = season,
            riverN = as.numeric(riverOrdered),
            yearN = year - min(year) + 1
          )

  dddd <- left_join( dddd, phi ) %>%
    dplyr::select(-speciesN,-seasonN,-riverN,-yearN)

  return(dddd)
}
