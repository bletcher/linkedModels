##################
# Nimble

plotInt_Nimble <- function(d){
  ggGrInt <- array2df(d$sims.list$grInt, label.x = "est") %>%
    #  rename(yoy=d2,species=d3,season=d4,river=d5)
    rename(yoy=d2,season=d3,river=d4)

  ggGrInt$chain <- rep(1:mcmcInfo$nChains, each = mcmcInfo$nSamples/mcmcInfo$nChains)
  ggGrInt$iter <- 1:as.numeric(mcmcInfo$nSamples/mcmcInfo$nChains)

  gg <- ggplot(filter(ggGrInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + ylim(-1.5,1.5) +
    #   facet_grid(d2+d4~d5+d3)
    facet_grid(yoy + river ~ season )
  print(gg)
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotBetas_Nimble <- function(d,b){
  ggGrBeta <- array2df(d$sims.list$grBeta, label.x = "est")

  ggGrBeta$chain <- rep(1:mcmcInfo$nChains, each = mcmcInfo$nSamples/mcmcInfo$nChains)
  ggGrBeta$iter <- 1:as.numeric(mcmcInfo$nSamples/mcmcInfo$nChains)

  gg <- list()
  numBetas <- 18
  for (i in 1:numBetas){
    gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) +
      geom_hline(yintercept = 0) +
      geom_point( aes(color = factor(chain)), size = 0.1 ) +
      ylim(-1,1) +
      facet_grid(d3+d5 ~ d4) +
      ggtitle(paste("beta =", i))
    if(i %in% b) print(gg[[i]])
  }
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotSigmaInt_Nimble <- function(d){
  ggSigmaInt <- array2df(d$sims.list$sigmaInt, label.x = "est")

  ggSigmaInt$chain <- rep(1:mcmcInfo$nChains, each = mcmcInfo$nSamples/mcmcInfo$nChains)
  ggSigmaInt$iter <- 1:as.numeric(mcmcInfo$nSamples/mcmcInfo$nChains)

  gg <- ggplot(filter(ggSigmaInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + ylim(-5,5) + facet_grid(d2+d4~d3)
  print(gg)
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotSigmaBetas_Nimble <- function(d,b){
  ggSigmaBeta <- array2df(d$sims.list$sigmaBeta, label.x = "est")

  ggSigmaBeta$chain <- rep(1:mcmcInfo$nChains, each = mcmcInfo$nSamples/mcmcInfo$nChains)
  ggSigmaBeta$iter <- 1:as.numeric(mcmcInfo$nSamples/mcmcInfo$nChains)

  gg <- list()
  numBetas <- 6
  for (i in 1:numBetas){
    gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) +
      ylim(-1,1) + facet_grid(d3+d5~d4) + ggtitle(paste("beta =", i))
    if(i %in% b) print(gg[[i]])
  }
}


##################
# jags

plotInt <- function(){
  ggGrInt <- array2df(dG[[1]]$sims.list$grInt, label.x = "est") %>%
  #  rename(yoy=d2,species=d3,season=d4,river=d5)
    rename(yoy=d2,season=d3,river=d4)

  ggGrInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggGrInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  gg <- ggplot(filter(ggGrInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + ylim(-1.5,1.5) +
    #   facet_grid(d2+d4~d5+d3)
    facet_grid(yoy + river ~ season )
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
    gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) +
      geom_hline(yintercept = 0) +
      geom_point( aes(color = factor(chain)), size = 0.1 ) +
      ylim(-1,1) +
      facet_grid(d3+d5 ~ d4) +
      ggtitle(paste("beta =", i))
    if(i %in% b) print(gg[[i]])
  }
}

# isYOY[ evalRows[i] ],species[ evalRows[i]],season[ evalRows[i] ],riverDATA[ evalRows[i] ]
# [1:675, 1, 1:2, 1:2, 1:4, 1:4]
plotSigma <- function(){
  ggSigmaInt <- array2df(dG[[1]]$sims.list$grSigma, label.x = "est")

  ggSigmaInt$chain <- rep(1:dG[[1]]$mcmc.info$n.chains, each = dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  ggSigmaInt$iter <- 1:as.numeric(dG[[1]]$mcmc.info$n.samples/dG[[1]]$mcmc.info$n.chains)
  gg <- ggplot(filter(ggSigmaInt), aes(iter,est)) + geom_point( aes(color=factor(chain)), size = 0.1 ) + ylim(-5,5) + facet_grid(d2+d4~d3)
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
    gg[[i]] <- ggplot(filter(ggGrBeta,d2 == i), aes(iter,est)) + geom_hline(yintercept = 0) + geom_point( aes(color = factor(chain)), size = 0.1 ) +
      ylim(-1,1) + facet_grid(d3+d5~d4) + ggtitle(paste("beta =", i))
    if(i %in% b) print(gg[[i]])
  }
}


# from Daniel Turek
samplesPlot <- function(samples, var=colnames(samples), ind=NULL, burnin=NULL, width=7, height=4, legend=TRUE, legend.location='topright', traceplot=TRUE, densityplot=TRUE, file=NULL) {
  if(!is.null(file)) pdf(file, width=width, height=height) else
    if(inherits(try(knitr::opts_chunk$get('dev'), silent=TRUE), 'try-error') || is.null(knitr::opts_chunk$get('dev')))   ## if called from Rmarkdown/knitr
      dev.new(height=height, width=width)
  par.save <- par(no.readonly = TRUE)
  par(mfrow=c(1,traceplot+densityplot), cex=0.7, cex.main=1.5, cex.axis=0.9, lab=c(3,3,7), mgp=c(0,0.4,0), mar=c(1.6,1.6,2,0.6), oma=c(0,0,0,0), tcl=-0.3, bty='l')
  ## process samples
  var <- gsub('\\[', '\\\\\\[', gsub('\\]', '\\\\\\]', var))   ## add \\ before any '[' or ']' appearing in var
  var <- unlist(lapply(var, function(n) grep(paste0('^', n,'(\\[.+\\])?$'), colnames(samples), value=TRUE)))  ## expanded any indexing
  samples <- samples[, var, drop=FALSE]
  if(!is.null(ind) && !is.null(burnin)) stop('only specify either ind or burnin')
  if(!is.null(ind))     samples <- samples[ind, , drop=FALSE]
  if(!is.null(burnin))  samples <- samples[(burnin+1):dim(samples)[1], , drop=FALSE]
  nparam <- ncol(samples)
  rng <- range(samples)
  if(!traceplot & !densityplot) stop('both traceplot and densityplot are false')
  if(traceplot) {  ## traceplot
    plot(1:nrow(samples), ylim=rng, type='n', main='Traceplots', xlab='', ylab='')
    for(i in 1:nparam)
      lines(samples[,i], col=rainbow(nparam, alpha=0.75)[i])
    if(legend & !densityplot & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
      legend(legend=dimnames(samples)[[2]], fill=rainbow(nparam, alpha=0.5), bty='n', x=legend.location)
  }  ## finish traceplot
  if(densityplot) {  ## denstyplot
    xMin <- xMax <- yMax <- NULL
    for(i in 1:nparam) {
      d <- density(samples[,i])
      xMin <- min(xMin,d$x); xMax <- max(xMax,d$x); yMax <- max(yMax, d$y) }
    plot(1, xlim=c(xMin,xMax), ylim=c(0,yMax), type='n', main='Posterior Densities', xlab='', ylab='', yaxt='n')
    for(i in 1:nparam)
      polygon(density(samples[,i]), col=rainbow(nparam, alpha=0.2)[i], border=rainbow(nparam, alpha=0.2)[i])
    if(legend & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
      legend(legend=dimnames(samples)[[2]], fill=rainbow(nparam, alpha=0.5), bty='n', x=legend.location)
  }  ## finish densityplot
  if(!is.null(file)) dev.off()
  invisible(par(par.save))
}
