library(lme4)

# lm to compare against the jags model
# patterns generally match up well

# d <- data.frame(
#   sampleNumber=dddG[[1]]$sampleNumber,
#   ind=ddG[[1]]$ind,
#   river=dddG[[1]]$riverOrdered,
#   season=dddG[[1]]$season,
#   species=dddG[[1]]$species,
#   isYOY=ddG[[1]]$isYOYDATA,
#   length=ddG[[1]]$lengthDATA,
#   grLength=dddG[[1]]$grLength,
#   tempStd=ddG[[1]]$tempStd,
#   flowStd=ddG[[1]]$flowStd,
#   countPStd=ddG[[1]]$countPStd)

d2 <- data.frame(
  #sampleNumber=ddG[[1]][[2]]$sampleNumber,
  ind=ddG[[1]][[1]]$ind,
  river=ddG[[1]][[1]]$riverDATA,
  season=ddG[[1]][[1]]$season,
  species=ddG[[1]][[1]]$species,
  isYOY=ddG[[1]][[1]]$isYOYDATA,
  length=ddG[[1]][[1]]$lengthDATA,
  grLength=ddG[[1]][[1]]$grNotUse,
  tempStd=ddG[[1]][[1]]$tempStd,
  flowStd=ddG[[1]][[1]]$flowStd,
  countPStd=ddG[[1]][[1]]$countPAllSppStd
)

# to run a new spp, need to rerun mainAnalysis.R from beginning to 'ddG[[ii]] <- dddG[[ii]] %>% prepareDataForJags("growth")'
d <- d2 %>% filter(species == as.numeric(speciesInGr), !is.na(grLength))
d <- d %>% mutate( tempStd2  =tempStd^2, flowStd2 = flowStd^2, countPStd2 = countPStd^2,
                   isYOY = as.factor(isYOY),season = season)
# model selection
# main effects

mod1 <- lmer( grLength ~ isYOY * river * season * ( tempStd ) + (1|ind), data = d )
mod2 <- lmer( grLength ~ isYOY * river * season * ( flowStd ) + (1|ind), data = d )
mod3 <- lmer( grLength ~ isYOY * river * season * ( countPStd ) + (1|ind), data = d )
mod4 <- lmer( grLength ~ isYOY * river * season * ( length ) + (1|ind), data = d )

mod100 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd * length ) + (1|ind), data = d )
mod102 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) + (1|ind), data = d ) # bkt, best model #
mod103 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd ) + (1|ind), data = d ) # bnt, best model #
mod104 <- lmer( grLength ~ isYOY * river * season * ( tempStd * countPStd ) + (1|ind), data = d )
mod105 <- lmer( grLength ~ isYOY * river * season * ( flowStd * countPStd) + (1|ind), data = d )

aicLinear <- AIC(mod1,mod2,mod3,mod4,mod100,mod102,mod103,mod104,mod105) %>% rownames_to_column() %>% arrange(AIC)

# add in squared term

mod2_1 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) +
                         isYOY * river * season * tempStd2+
                         isYOY * river * season * flowStd2+
                         isYOY * river * season * countPStd2 +
                         (1|ind), data = d )

mod2_2 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) +
             #   isYOY * river * season * tempStd2+
            #    isYOY * river * season * flowStd2+
            #    isYOY * river * season * countPStd2 +
                (1|ind), data = d )

mod2_3 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) +
                   isYOY * river * season * tempStd2+
                #    isYOY * river * season * flowStd2+
                #    isYOY * river * season * countPStd2 +
                (1|ind), data = d )
# bkt, best model #
mod2_4 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) +
                #   isYOY * river * season * tempStd2+
                   isYOY * river * season * flowStd2+
                #    isYOY * river * season * countPStd2 +
                (1|ind), data = d )

mod2_5 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) +
                #isYOY * river * season * tempStd2+
                #   isYOY * river * season * flowStd2+
                    isYOY * river * season * countPStd2 +
                (1|ind), data = d )

mod2_6 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) +
                isYOY * river * season * tempStd2+
                isYOY * river * season * flowStd2+
                #isYOY * river * season * countPStd2 +
                (1|ind), data = d )

mod2_7 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) +
                isYOY * river * season * tempStd2+
                #isYOY * river * season * flowStd2+
                isYOY * river * season * countPStd2 +
                (1|ind), data = d )
# bnt, best model #
mod2_8 <- lmer( grLength ~ isYOY * river * season * ( tempStd * flowStd * countPStd ) +
                #isYOY * river * season * tempStd2+
                isYOY * river * season * flowStd2+
                isYOY * river * season * countPStd2 +
                (1|ind), data = d )

mod2_100 <- lmer( grLength ~ isYOY * river * season * ( tempStd + flowStd + countPStd) +

                #   isYOY * river * season * tempStd2+
                #    isYOY * river * season * flowStd2+
                #    isYOY * river * season * countPStd2 +

                (1|ind), data = d )

aicSquared <- AIC(mod2_1,mod2_2,mod2_3,mod2_4,mod2_5,mod2_6,mod2_7,mod2_8,mod2_100) %>% rownames_to_column() %>% arrange(AIC)

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

ggplot(d, aes(tempStd, grLength, color=countPStd)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(season ~ river)

# simple lm
d <- d %>% mutate( tempStd2  =tempStd^2, flowStd2 = flowStd^2, countPStd2 = countPStd^2)
mod1 <- lm( grLength ~ isYOY * river * season * tempStd * flowStd * countPStd + tempStd2+flowStd2+countPStd2, data = d )

# predictions for the data
d$predicted <- predict( mod1 )
ggplot(d, aes(tempStd, predicted, color=countPStd)) +
  geom_point() +
  facet_grid(season ~ river)

ggplot(d, aes(tempStd, predicted, color=flowStd)) +
  geom_point() +
  facet_grid(season ~ river)


# predict over a grid
x <- seq( -2,2,length.out = 5 )

  isYOYData <- rep(1:4, each = nPoints ^ 5)
  seasonData <- rep(1:4, each = nPoints ^ 4)
  riverData <- rep(unique(d$river), each = nPoints ^ 3)
  tempData <- rep(x, each = nPoints ^ 2)
  flowData <- rep(x, each = nPoints ^ 1)
  countData <- rep(x, each = nPoints ^ 0)


predTemplate <- data.frame( isYOY = isYOYData,
                            season = seasonData,
                            river = riverData,
                            countPStd = countData,
                            flowStd =  flowData,
                            tempStd =  tempData
) %>% filter(isYOY %in% 1:2)
predTemplate <- predTemplate %>% mutate( tempStd2  = tempStd^2, flowStd2 = flowStd^2, countPStd2 = countPStd^2)

predTemplate$predicted <- predict( mod1, newdata=predTemplate )

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


