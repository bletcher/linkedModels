#'Estimate proportion overlap for betas
#'
#'@param mcmcProcessed a list generated using process.output()
#'@param nB a data frame containing the number of betas for the run
#'@return Prints the proportion of [,isYOY,season,river] beta means that overlap 0
#'@export

getPropOverlap0Betas <- function(mcmcProcessed,nB){

  getProp <- function(i) {
    sum(i*1)/length(i)
  }

#  lapply( mcmcProcessed$overlap0, getProp )
#  lapply( mcmcProcessed$mean, mean )

  # for each beta
  print("grBeta")
  for ( i in 1:nB$nBetas ){
    prop <- lapply( mcmcProcessed$overlap0$grBeta[i,,,], getProp )
    print( sum(data.frame(prop)) / length(data.frame(prop)) )
  }
  print("")

  # for each BNT beta
  print("grBetaBNT")
  for ( i in 1:nB$nBetasBNT ){
    prop <- lapply( mcmcProcessed$overlap0$grBetaBNT[i,,,1:2], getProp )
    print( sum(data.frame(prop)) / length(data.frame(prop)) )
  }
  print("")

  # for each ATS beta
  print("grBetaATS")
  for ( i in 1:nB$nBetasATS ){
    prop <- lapply( mcmcProcessed$overlap0$grBetaATS[i,,,1], getProp )
    print( sum(data.frame(prop)) / length(data.frame(prop)) )
  }
}
