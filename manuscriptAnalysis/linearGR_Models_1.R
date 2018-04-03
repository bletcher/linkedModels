# library(tidyverse)
# #d$ATS01 <- ifelse( d$sampleNumber <= 58,1,0 )
# d <- d %>% mutate( tempStd2  =tempStd^2, flowStd2 = flowStd^2, countPStd2 = countPStd^2, length2 = length^2)
#
# #speciesIn <- 'bkt'
# #d <- d %>% filter(species == speciesIn, !is.na(grLength))
#
# ggplot(d, aes(tempStd, grLength)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   facet_grid(season ~ river)
#
# ggplot(d, aes(countPStd, grLength)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   facet_grid(season ~ species)
#
#
# pairs(d[,c("nBKT","nBNT","nATS","countPStd")])
#
# ggplot(d, aes(nBKT, grLength)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   facet_grid(season ~ species)
# ggplot(d, aes(nBNT, grLength)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   facet_grid(season ~ species)
# ggplot(d, aes(nATS, grLength)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   facet_grid(season ~ species)
#
#
# ggplot(d, aes(flowStd, grLength)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   facet_grid(season ~ species)
#
# ggplot(d, aes(tempStd, grLength, color=countPStd)) +
#   geom_point() +
#   geom_smooth(method='lm') +
#   facet_grid(season ~ river)
#
# # simple lm
# mod0 <- lm( grLength ~ as.factor(isYOY) * river * as.factor(season) * tempStd * flowStd * countPStd + tempStd2+flowStd2+countPStd2, data = d )
#
# # lmer"
#
# mod1 <- lmer( grLength ~ as.factor(isYOY) + river + as.factor(season) +
#                          tempStd + flowStd + countPStd +
#                          tempStd2 + flowStd2 + countPStd2 + (1|ind), data = d)
#
# mod2 <- lmer( grLength ~ as.factor(isYOY) + river + as.factor(season) +
#                 length + tempStd + flowStd + countPStd +
#                 tempStd2 + flowStd2 + countPStd2 + length2 +
#                 tempStd * flowStd +
#                 flowStd * countPStd +
#                 tempStd * countPStd +
#                 tempStd * flowStd * countPStd +
#                 (1|ind), data = d)
#
# mod3 <- lmer( grLength ~ as.factor(isYOY) + river + as.factor(season) +
#                 length + tempStd + flowStd + countPStd +
#
#                 tempStd2 + flowStd2 + countPStd2 + length2 +
#                 tempStd * flowStd +
#                 flowStd * countPStd +
#                 tempStd * countPStd +
#                 tempStd * flowStd * countPStd +
#
#                 (1|ind), data = d)
#
# mod4 <- lmer( grLength ~ as.factor(isYOY) + river + as.factor(season) +
#                 length + tempStd + flowStd + nBKT + nBNT + nATS +
#                 tempStd2 + flowStd2 + countPStd2 + length2 +
#                 tempStd * flowStd +
#                 flowStd * countPStd +
#                 tempStd * countPStd +
#                 tempStd * flowStd * countPStd +
#
#                 (1|ind), data = d)
#
# mod5 <- lmer( grLength ~ as.factor(isYOY) + river + as.factor(season) +
#                 length + tempStd + flowStd + nBKT + nBNT + nATS +
#                 tempStd2 + flowStd2 + countPStd2 + length2 +
#                 tempStd * flowStd +
#                 flowStd * nBKT +
#                 flowStd * nBNT +
#                 flowStd * nATS +
#                 tempStd * nBKT +
#                 tempStd * nBNT +
#                 tempStd * nATS +
#                 tempStd * flowStd * nBKT*nBNT*nATS +
#
#                 (1|ind), data = d)
#
# mod5a <- lmer( grLength ~ as.factor(isYOY) + river + as.factor(season) +
#                 length + tempStd + flowStd + nBKT + nBNT + nATS +
#                 tempStd2 + flowStd2 + countPStd2 + length2 +
#                 tempStd * flowStd * nBKT*nBNT*nATS +
#
#                 (1|ind), data = d)
#
#
# AIC(mod1,mod2,mod3,mod4,mod5,mod5a)
#
#
#
#
#
#
#
#
#
#
#
#
#
# # predictions for the data
# d$predicted <- predict( mod1 )
# ggplot(d, aes(tempStd, predicted, color=countPStd)) +
#   geom_point() +
#   facet_grid(season ~ river)
#
# ggplot(d, aes(tempStd, predicted, color=flowStd)) +
#   geom_point() +
#   facet_grid(season ~ river)
#
#
# # predict over a grid
# x <- seq( -2,2,length.out = 5 )
#
#   isYOYData <- rep(1:4, each = nPoints ^ 5)
#   seasonData <- rep(1:4, each = nPoints ^ 4)
#   riverData <- rep(unique(d$river), each = nPoints ^ 3)
#   tempData <- rep(x, each = nPoints ^ 2)
#   flowData <- rep(x, each = nPoints ^ 1)
#   countData <- rep(x, each = nPoints ^ 0)
#
#
# predTemplate <- data.frame( isYOY = isYOYData,
#                             season = seasonData,
#                             river = riverData,
#                             countPStd = countData,
#                             flowStd =  flowData,
#                             tempStd =  tempData
# ) %>% filter(isYOY %in% 1:2)
# predTemplate <- predTemplate %>% mutate( tempStd2  = tempStd^2, flowStd2 = flowStd^2, countPStd2 = countPStd^2)
#
# predTemplate$predicted <- predict( mod1, newdata=predTemplate )
#
# countIn <- 0
# isYOYIn <- 2
# ggplot(predTemplate %>% filter(countPStd == countIn, isYOY == isYOYIn), aes(tempStd, predicted, color=flowStd)) +
#   geom_line(aes(group=flowStd)) +
#   geom_point() +
#   ylim(-0.2,1) +
#   ggtitle(paste(countIn, isYOYIn,speciesIn)) +
#   facet_grid( river ~ season )
#
# tempIn <- 0
# isYOYIn <- 2
# ggplot(predTemplate %>% filter(tempStd == tempIn, isYOY == isYOYIn), aes(countPStd, predicted, color=flowStd)) +
#   geom_line(aes(group=flowStd)) +
#   geom_point() +
#   ylim(-0.2,1) +
#   ggtitle(paste(tempIn, isYOYIn,speciesIn)) +
#   facet_grid( river ~ season )
#
# tempIn <- 0
# flowIn=0
# isYOYIn <- 2
# ggplot(predTemplate %>% filter(tempStd == tempIn,flowStd==flowIn, isYOY == isYOYIn), aes(countPStd, predicted)) +
#   geom_line() +
#   geom_point() +
#   ylim(-0.2,1) +
#   ggtitle(paste(tempIn, isYOYIn,speciesIn)) +
#   facet_grid( river ~ season )
#
#
