#From Evan's Ecology paper'

d <- data.frame(temp = seq(1,20.5,0.1))

# bkt first, bnt second

tOpt <- c(11.9,13.3)
cTMax <- c(19.5,19.3)
sigma <- c(3.7,4.1)

# lInf <- 354
# k <- 0.0013
# flowBeta <- 0.098
# bktBeta <- -0.051
# bntBeta <- -0.039
#
# bktDen <- 7.01
# bntDen <- 13.00
# flow <- 0.14
#
# gOpt <- k*(lInf-len) + flowBeta*flow + bktBeta*bktDen + bntBeta*bntDen

i <- 1
d$pT_BKT <- ifelse( d$temp <= tOpt[i],
                 exp(-((d$temp - tOpt[i])/(2*sigma[i]))^2), #there's an error in the paper
                 1 - ((d$temp - tOpt[i])/(tOpt[i] - cTMax[i]))^2
              )

i <- 2
d$pT_BNT <- ifelse( d$temp <= tOpt[i],
                    exp(-((d$temp - tOpt[i])/(2*sigma[i]))^2), #there's an error in the paper
                    1 - ((d$temp - tOpt[i])/(tOpt[i] - cTMax[i]))^2
)

ggplot(d, aes(temp,pT_BKT)) + geom_line() + geom_line(aes(temp,pT_BNT), color = "brown")
