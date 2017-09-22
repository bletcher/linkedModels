

#'Run the growth model
#'
#'@param drainage Which drainage, "west" or "stanley"
#'@return a data frame
#'@export

runGrowthModel <- function(d, parallel = FALSE){   #iterToUse, firstNonBurnIter, chainToUse, simInfo, grDir, out){

  inits <- function(){
    #list(grBetaInt = array(rnorm(d$nSeasons*d$nRivers,0,2.25),c(d$nSeasons,d$nRivers)))
    #list(grBetaInt = array(rnorm(2*4*4*6,0,2.25),c(2,4,4,6)))
    list(grBetaInt = array(rnorm(2*dddG$nSeasons*dddG$nRivers*dddG$nYears,0,2.25),c(2,dddG$nSeasons,dddG$nRivers,dddG$nYears)))
  }

  params <- c("grBetaInt","muGrBetaInt","sigmaGrBetaInt","grBeta","muGrBeta","grSigmaBeta","length")

  outGR <- jags(data = d,
                inits = inits,
                parameters.to.save = params,
                model.file = "./jags/grModel.jags",
                n.chains = 3,
                n.adapt = 1000, #1000
                n.iter = 2000,
                n.burnin = 1000,
                n.thin = 4,
                parallel = parallel
  )

  #outGR$movementModelIterUsed <- iter
  return(outGR)
}
