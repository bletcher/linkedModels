plotInt <- function(){
  ggGrInt <- array2df(dG[[1]]$sims.list$grInt, label.x = "est") %>%
    rename(yoy=d2,species=d3,season=d4,river=d5)

  ggGrInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggGrInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  gg <- ggplot(filter(ggGrInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + ylim(-1.5,1.5) +
    #   facet_grid(d2+d4~d5+d3)
    facet_grid(yoy + river ~ season + species)
  print(gg)
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotBetas <- function(b){
  ggGrBeta <- array2df(dG[[1]]$sims.list$grBeta, label.x = "est")

  ggGrBeta$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggGrBeta$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)

  gg <- list()
  numBetas <- 18
  for (i in 1:numBetas){
    gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) + ylim(-1,1) + facet_grid(d3+d5~d6+d4) + ggtitle(paste("beta =", i))
    if(i %in% b) print(gg[[i]])
  }
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotSigma <- function(){
  ggSigmaInt <- array2df(dG[[1]]$sims.list$grSigma, label.x = "est")

  ggSigmaInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggSigmaInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  gg <- ggplot(filter(ggSigmaInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + ylim(-5,5) + facet_grid(d2+d4~d5+d3)
  print(gg)
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotSigmaBetas <- function(b){
  ggGrBeta <- array2df(dG[[1]]$sims.list$sigmaBeta, label.x = "est")

  ggGrBeta$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggGrBeta$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)

  gg <- list()
  numBetas <- 4
  for (i in 1:numBetas){
    gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) + ylim(-1,1) + facet_grid(d3+d5~d6+d4) + ggtitle(paste("beta =", i))
    if(i %in% b) print(gg[[i]])
  }
}
