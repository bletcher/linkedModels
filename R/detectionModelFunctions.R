#'Run the detection model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runDetectionModel <- function(d, parallel = FALSE){   #iterToUse, firstNonBurnIter, chainToUse, simInfo, grDir, out){

  inits <- function(){
    list(#pBetaInt = array(rnorm(dddD$nSeasons*dddD$nRivers*dddD$nYears,0,2.25),c(dddD$nSeasons,dddD$nRivers,dddD$nYears)),
         pBetaInt = array(runif(dddD$nSeasons*dddD$nRivers*dddD$nYears, -2.5, 2),c(dddD$nSeasons,dddD$nRivers,dddD$nYears)),
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
