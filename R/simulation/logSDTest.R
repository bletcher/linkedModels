library(nimble)
library(tidyverse)

mc <- nimbleCode({

  b ~ dnorm(0, sd = 20)
  m ~ dnorm(0, sd = 20)
  sigmaInt ~ dnorm(0, sd = 20)
  sigmaSlope ~ dnorm(0, sd = 20)

  for(i in 1:nrow){
    y[i] ~ dnorm(mean[i], sd = sigmaExp[i])
    mean[i] <- b + m * x[i]
    log(sigmaExp[i]) <- sigmaInt + sigmaSlope * x[i]
  }

})

m <- 5
b <- 10
sigmaInt <- 0.2
sigmaSlope <- 2
nReps <- 1000
xes <- seq(0,4,by=0.5)
d <- data.frame( xIn = rep(xes,nReps) )
d$y <- rnorm(length(xes)*nReps, b + m * d$xIn, sigmaInt + sigmaSlope * d$xIn)
offset <- 0
d$x <- d$xIn + offset

model <- nimbleModel(mc, data = list( y = d$y, x = d$x), constants = list(nrow = nrow(d)), inits = list(b=rnorm(1,b),m=rnorm(1,m),sigmaInt=rnorm(1,sigmaInt),sigmaSlope=rnorm(1,sigmaSlope)))

out <- nimbleMCMC(model=model, monitors = c( 'b','m','sigmaInt','sigmaSlope'),
                       nchains = 2, niter = 2000, nburnin = 1000,
                       summary = TRUE, WAIC = TRUE)
out$summary

means <- data.frame(x = b+m*(unique(d$x)), y = (out$summary$all.chains['b','Mean']+out$summary$all.chains['m','Mean']*(unique(d$x) + offset)))
ggplot(means, aes(x,y)) + geom_point() + geom_abline(slope=1,intercept=0)

sds <- data.frame( x = sigmaInt+sigmaSlope*(unique(d$x)), y = exp(out$summary$all.chains['sigmaInt','Mean'] + out$summary$all.chains['sigmaSlope','Mean']*(unique(d$x) + offset)) )
ggplot(sds, aes(x,y)) + geom_point() + geom_abline(slope=1,intercept=0)


