#'Run the detection model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runDetectionModel <- function(d, parallel = FALSE){   #iterToUse, firstNonBurnIter, chainToUse, simInfo, grDir, out){

  inits <- function(){
    list(pBetaInt = array(rnorm(2*ddd$nSeasons*ddd$nRivers*ddd$nYears,0,2.25),c(2,ddd$nSeasons,ddd$nRivers,ddd$nYears)))
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
