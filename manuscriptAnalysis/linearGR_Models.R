library(lme4)
library(tidyverse)

# lm to compare against the jags model
# patterns generally match up well

load(file = paste0('./data/out/dG_bktbntats1997_forLmer.RData')) # made by running mainAnalysis.R with speciesGr = c("bkt", "bnt","ats") through ddG[[ii]] <- dddG[[ii]] %>% prepareDataForJags("growth")

d2 <- data.frame(
  #sampleNumber=ddG[[1]][[2]]$sampleNumber,
  ind=ddG[[1]][[1]]$ind,
  river=ddG[[1]][[1]]$riverDATA,
  season=ddG[[1]][[1]]$season,
  species=ddG[[1]][[1]]$species,
  isYOY=ddG[[1]][[1]]$isYOYDATA,
  length=ddG[[1]][[1]]$lengthDATAStd,
  lengthByYear=ddG[[1]][[1]]$lengthDATAStd_ByYOY_Year,
  grLength=ddG[[1]][[1]]$grNotUse,
  tempStd=ddG[[1]][[1]]$tempStd,
  flowStd=ddG[[1]][[1]]$flowStd,
  countPStd=ddG[[1]][[1]]$countPAllSppStd,
  countPStdBKT=ddG[[1]][[1]]$countPStdBKT,
  countPStdBNT=ddG[[1]][[1]]$countPStdBNT,
  countPStdATS=ddG[[1]][[1]]$countPStdATS
)

d <- d2 %>% filter(
  #  species == as.numeric(speciesInGr),
    !is.na(grLength)
  )
d <- d %>% mutate( tempStd2 = tempStd^2, flowStd2 = flowStd^2, countPStd2 = countPStd^2, length2 = length^2,
                   isYOY = as.factor(isYOY), season = as.factor(season), river = as.factor(river), species = as.factor(species)
                 )
d$grLength <- ifelse(d$species == 2 & d$river == 3,NA,d$grLength)

#################################################################
# ?s use length std by sample# and isYOY?
#    use abundance std by isYOY?

# model with length is a lot better. Stick with that to start
modl0 <- lmer( grLength ~ isYOY * species * river * season * length + (1|ind), data = d )
modl1 <- lmer( grLength ~ isYOY * species * river * season * lengthByYear + (1|ind), data = d )
ggplot(d,aes(length,lengthByYear,color=isYOY)) + geom_point()
AIC(modl0,modl1)

# model with spp-specific abundances is a lot better. Stick with that to start
moda0 <- lmer( grLength ~ isYOY * species * river * season * countPStd + (1|ind), data = d )
moda1 <- lmer( grLength ~ isYOY * species * river * season * countPStdBKT * countPStdBNT * countPStdATS + (1|ind), data = d )
AIC(moda0,moda1)


# model selection
# main effects
mod0 <- lmer( grLength ~ isYOY * species * river * season + (1|ind), data = d )
mod1 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd ) + (1|ind), data = d )
mod2 <- lmer( grLength ~ isYOY * species * river * season * ( flowStd ) + (1|ind), data = d )
mod3 <- lmer( grLength ~ isYOY * species * river * season * ( countPStdBKT * countPStdBNT * countPStdATS ) + (1|ind), data = d )
mod4 <- lmer( grLength ~ isYOY * species * river * season * ( length ) + (1|ind), data = d )

mod100 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS * length ) + (1|ind), data = d )
mod102 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + (1|ind), data = d )
mod103 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd ) + (1|ind), data = d )
mod104 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * countPStdBKT * countPStdBNT * countPStdATS ) + (1|ind), data = d )
mod105 <- lmer( grLength ~ isYOY * species * river * season * ( flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + (1|ind), data = d )

#additive better?, no
mod102a <- lmer( grLength ~ isYOY * species * river * season * ( tempStd + flowStd + countPStdBKT + countPStdBNT + countPStdATS ) + (1|ind), data = d )

# length has clear effect on graphs. what structure is needed - full interactive (mod100) worse than model w/o length (mod102)
mod102b <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + isYOY * length + (1|ind), data = d )# bkt ats, best model #
mod102c <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + isYOY * species * length + (1|ind), data = d )# bkt ats, best model #
mod102d <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + isYOY * species * river * length + (1|ind), data = d )# bkt ats, best model #
mod102e <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + isYOY * species * river * season * length + (1|ind), data = d )# bkt ats, best model #

# which interactions with length?
mod102e1 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + isYOY * species * river * season * length * tempStd + (1|ind), data = d )# bkt ats, best model #
mod102e2 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + isYOY * species * river * season * length * flowStd + (1|ind), data = d )# bkt ats, best model #
mod102e3 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStdBKT * countPStdBNT * countPStdATS ) + isYOY * species * river * season * length * countPStdBKT * countPStdBNT * countPStdATS + (1|ind), data = d )# bkt ats, best model #

aicLinear <- AIC(mod102e1,mod102e2,mod102e3,mod102a,mod102b,mod102c,mod102d,mod102e,mod0,mod1,mod2,mod3,mod4,mod100,mod102,mod103,mod104,mod105) %>% rownames_to_column() %>% arrange(AIC)

##########################
# add in squared terms

mod2_1 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
                           isYOY * species * river * season * tempStd2+
                           isYOY * species * river * season * flowStd2+
                           isYOY * species * river * season * countPStd2 +
                           (1|ind), data = d )

mod2_1a <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
                  isYOY * species * river * season * tempStd2+
                  isYOY * species * river * season * flowStd2+
                  isYOY * species * river * season * countPStd2 +
                  isYOY * species * river * season * length2 +
                  (1|ind), data = d )

mod2_2 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
             #   isYOY * species * river * season * tempStd2+
            #    isYOY * species * river * season * flowStd2+
            #    isYOY * species * river * season * countPStd2 +
                (1|ind), data = d )

mod2_3 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
                   isYOY * species * river * season * tempStd2+
                #    isYOY * species * river * season * flowStd2+
                #    isYOY * species * river * season * countPStd2 +
                (1|ind), data = d )

mod2_4 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
                #   isYOY * species * river * season * tempStd2+
                   isYOY * species * river * season * flowStd2+
                #    isYOY * species * river * season * countPStd2 +
                (1|ind), data = d )

mod2_5 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
                #isYOY * species * river * season * tempStd2+
                #   isYOY * species * river * season * flowStd2+
                    isYOY * species * river * season * countPStd2 +
                (1|ind), data = d )

mod2_6 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
                isYOY * species * river * season * tempStd2+
                isYOY * species * river * season * flowStd2+
                #isYOY * species * river * season * countPStd2 +
                (1|ind), data = d )

mod2_7 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
                isYOY * species * river * season * tempStd2+
                #isYOY * species * river * season * flowStd2+
                isYOY * species * river * season * countPStd2 +
                (1|ind), data = d )

mod2_8 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd * flowStd * countPStd ) + isYOY * species * river * season * length +
                #isYOY * species * river * season * tempStd2+
                isYOY * species * river * season * flowStd2+
                isYOY * species * river * season * countPStd2 +
                (1|ind), data = d )


mod2_100 <- lmer( grLength ~ isYOY * species * river * season * ( tempStd + flowStd + countPStd * length) +
                #   isYOY * species * river * season * tempStd2+
                #    isYOY * species * river * season * flowStd2+
                #    isYOY * species * river * season * countPStd2 +
                (1|ind), data = d )

aicSquared <- AIC(mod2_1,mod2_1a,mod2_2,mod2_3,mod2_4,mod2_5,mod2_6,mod2_7,mod2_8,mod2_100) %>% rownames_to_column() %>% arrange(AIC)
aicSquared <- AIC(mod2_1,mod2_1a,mod2_2,mod2_3,mod2_5,mod2_6,mod2_7,mod2_100) %>% rownames_to_column() %>% arrange(AIC)

aicLinear
aicSquared

ggplot(d, aes(tempStd, grLength)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(season ~ river)

ggplot(d, aes(countPStd, grLength)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(season ~ species)

ggplot(d, aes(flowStd, grLength)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(season ~ species)

ggplot(d, aes(tempStd, grLength, color=countPStd )) +
  geom_point() +
  geom_smooth(method='lm') + ylim(-0.2,1) +
  facet_grid(season ~ river)

ggplot(d, aes(length, grLength, color=tempStd )) +
  geom_point() +
  geom_smooth(method='lm') + ylim(-0.2,1) +
  facet_grid(season ~ river)

# simple lm
#d <- d %>% mutate( tempStd2  =tempStd^2, flowStd2 = flowStd^2, countPStd2 = countPStd^2)
#mod1 <- lm( grLength ~ isYOY * river * season * tempStd * flowStd * countPStd + tempStd2+flowStd2+countPStd2, data = d )

# predictions for the data
d$predicted <- predict( mod2_1 )
ggplot(d, aes(tempStd, predicted, color=countPStd )) +
  geom_point() +
  facet_grid(season ~ river)

ggplot(d, aes(tempStd, predicted, color=flowStd)) +
  geom_point() +
  facet_grid(season ~ river)


# predict over a grid
nPoints <- 4
x <- seq( -2,2,length.out = 5 )

  lengthData <- rep(x, each = nPoints ^ 6)
  isYOYData <- rep(1:4, each = nPoints ^ 5)
  seasonData <- rep(1:4, each = nPoints ^ 4)
  riverData <- rep(unique(d$river), each = nPoints ^ 3)
  tempData <- rep(x, each = nPoints ^ 2)
  flowData <- rep(x, each = nPoints ^ 1)
  countData <- rep(x, each = nPoints ^ 0)


predTemplate <- data.frame( length = lengthData,
                            isYOY = isYOYData,
                            season = seasonData,
                            river = riverData,
                            countPStd = countData,
                            flowStd =  flowData,
                            tempStd =  tempData
) %>% filter(isYOY %in% 1:2)
predTemplate <- predTemplate %>% mutate( tempStd2  = tempStd^2, flowStd2 = flowStd^2, countPStd2 = countPStd^2, length2 = length^2)

predTemplate$predicted <- predict( mod2_15, newdata = predTemplate )

countIn <- 0
isYOYIn <- 2
ggplot(predTemplate %>% filter(countPStd == countIn, isYOY == isYOYIn), aes(tempStd, predicted, color=flowStd)) +
  geom_line(aes(group=flowStd)) +
  geom_point() +
  ylim(-0.2,1) +
  ggtitle(paste(countIn, isYOYIn,speciesIn)) +
  facet_grid( river ~ season )

tempIn <- 0
isYOYIn <- 2
ggplot(predTemplate %>% filter(tempStd == tempIn, isYOY == isYOYIn), aes(countPStd, predicted, color=flowStd)) +
  geom_line(aes(group=flowStd)) +
  geom_point() +
  ylim(-0.2,1) +
  ggtitle(paste(tempIn, isYOYIn,speciesIn)) +
  facet_grid( river ~ season )

tempIn <- 0
flowIn=0
isYOYIn <- 2
ggplot(predTemplate %>% filter(tempStd == tempIn,flowStd==flowIn, isYOY == isYOYIn), aes(countPStd, predicted)) +
  geom_line() +
  geom_point() +
  ylim(-0.2,1) +
  ggtitle(paste(tempIn, isYOYIn,speciesIn)) +
  facet_grid( river ~ season )


