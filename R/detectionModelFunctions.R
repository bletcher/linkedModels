#'Run the detection model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runDetectionModel <- function(d, parallel = FALSE){   #iterToUse, firstNonBurnIter, chainToUse, simInfo, grDir, out){

  inits <- function(){
    list(#pBetaInt = array(rnorm(dddD$nSeasons*dddD$nRivers*dddD$nYears,0,2.25),c(dddD$nSeasons,dddD$nRivers,dddD$nYears)),
         pBetaInt = array(runif(dddD$nSpecies*dddD$nSeasons*dddD$nRivers*dddD$nYears, -2.5, 2),c(dddD$nSpecies,dddD$nSeasons,dddD$nRivers,dddD$nYears)),
         z = dddD$zForInit
         )
  }

  params <- c("pBetaInt","phiBetaInt")

  outDet <- jags(data = d,
                 inits = inits,
                 parameters.to.save = params,
                 model.file = "./jags/detModel.jags",
                 n.chains = 3,
                 n.adapt = 1000, #1000
                 n.iter = 2000,
                 n.burnin = 1000,
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

getDensities <- function(dddd,dd, meanOrIter = "mean", sampleToUse = 1, chainToUse = 1){

  try( if (sampleToUse > dd$mcmc.info$n.samples/dd$mcmc.info$n.chains) stop("requested iteration beyond max # of iterations") )
  try( if (chainToUse > dd$mcmc.info$n.chains) stop("requested chain beyond max # of chains") )


  counts <- dddd %>%
    filter( enc == 1 ) %>%
    group_by( species,season,riverOrdered,year ) %>%
    summarise( count = n() ) %>%
    mutate( speciesN = as.numeric(as.factor(species)),
            seasonN = season,
            riverN = as.numeric(riverOrdered),
            yearN = year - min(year) + 1
          ) %>%
    ungroup() #%>%
#    dplyr::select( speciesN, seasonN, riverN, yearN, count )


  if ( meanOrIter == 'mean') ddIn <- dd$q50$pBetaInt

  if ( meanOrIter == 'iter') {
    sampleN <- sampleToUse + dd$mcmc.info$n.samples/dd$mcmc.info$n.chains * (chainToUse - 1)
    ddIn <- dd$sims.list$pBetaInt[ sampleN,,,, ]
  }

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
