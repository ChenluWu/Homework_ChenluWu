rm(list = ls())
gc()

library(foreign)
library(dplyr)

setwd('/Users/wuchenlu/Desktop/')

UGAdata <- read.dta('dataUGA.dta')
dataNeed <- subset(UGAdata,select = c("hh","ctotal","inctotal","age","age_sq","ethnic","familysize","female","urban","wave","year"))

dataNeed$t[dataNeed$wave=='2009-2010']=1
dataNeed$t[dataNeed$wave=='2010-2011']=2
dataNeed$t[dataNeed$wave=='2011-2012']=3
dataNeed$t[dataNeed$wave=='2013-2014']=4

dataNeed <- dataNeed[!is.na(dataNeed$ctotal),]
dataNeed <- dataNeed[!is.na(dataNeed$inctotal),]

dataNeed$lnc <- log(dataNeed$ctotal)
dataNeed$lninc <- log(dataNeed$inctotal)

dataNeed <- dataNeed[!is.infinite(dataNeed$lnc),]
dataNeed <- dataNeed[!is.infinite(dataNeed$lninc),]

dataNeed <- na.omit(dataNeed)

resultc <- lm(lnc~1+age+age_sq+familysize+factor(t)+factor(ethnic)+factor(female)+factor(urban),data = dataNeed)
resultinc <- lm(lninc~1+age+age_sq+familysize+factor(t)+factor(ethnic)+factor(female)+factor(urban),data = dataNeed)

dataNeed$resc <- resultc$residuals
dataNeed$resinc <- resultinc$residuals

dataNeed$aggc[dataNeed$t==1] = log(mean(dataNeed$ctotal[dataNeed$t==1]))
dataNeed$aggc[dataNeed$t==2] = log(mean(dataNeed$ctotal[dataNeed$t==2]))
dataNeed$aggc[dataNeed$t==3] = log(mean(dataNeed$ctotal[dataNeed$t==3]))
dataNeed$aggc[dataNeed$t==4] = log(mean(dataNeed$ctotal[dataNeed$t==4]))

# dataNeed <- dataNeed[dataNeed$urban==0,]

hhi <- unique(dataNeed$hh)

res_question1 <- matrix(0,length(hhi),3)
alldlnc <- vector(length=0)
alldlninc <- vector(length=0)
alldlnaggc <- vector(length=0)


for (i in 1:length(hhi)) {
  id <- hhi[i]
  regi <- dataNeed[dataNeed$hh==id,]
  regi <- regi[order(regi$t),,]
  cnt <- nrow(regi)
  if (cnt>2) {
  dlnc <- vector(length = cnt-1)
  dlninc <- vector(length = cnt-1)
  dlnaggc <- vector(length = cnt-1)
  for (j in 1:cnt-1) {
    dlnc[j] <- (regi$resc[j+1]-regi$resc[j])/(regi$t[j+1]-regi$t[j])
    dlninc[j] <- (regi$resinc[j+1]-regi$resinc[j])/(regi$t[j+1]-regi$t[j])
    dlnaggc[j] <-(regi$aggc[j+1]-regi$aggc[j])/(regi$t[j+1]-regi$t[j])
  }
  resulti <- lm(dlnc ~ 0 + dlninc + dlnaggc)
  res_question1[i,] <-c(id,resulti$coefficients[1],resulti$coefficients[2])
  dataNeed$beta[dataNeed$hh==id] = resulti$coefficients[1]
  dataNeed$phi[dataNeed$hh==id] = resulti$coefficients[2]
  
  alldlnc <- c(alldlnc,dlnc)
  alldlninc <- c(alldlninc,dlninc)
  alldlnaggc <- c(alldlnaggc,dlnaggc)
  
  } else {
    res_question1[i,] <- c(id,NaN,NaN)
  }
}

final_q1 <- as.data.frame(na.omit(res_question1))

final_q1 <- filter(final_q1, V2>-10 & V2<10 & V3>-50 & V3<50)
colnames(final_q1) <- c('hh','beta','phi')

colMeans(final_q1)
quantile(final_q1[,2],probs = c(0.25,0.5,0.75,1))
quantile(final_q1[,3],probs = c(0.25,0.5,0.75,1))

hist(final_q1[,2], breaks = 100, xlim = c(-10,10), probability = TRUE, main = "histogram of beta")
hist(final_q1[,3], breaks = 100, xlim = c(-50,50), probability = TRUE, main = "histogram of phi")

data_group <- dataNeed %>% group_by(hh) %>% 
  summarise(aggY=mean(inctotal))%>% arrange(hh)

res_q2 <- left_join(final_q1,data_group,by='hh')
res_q2$absbeta <- abs(res_q2$beta)

incomeQuantile <- quantile(res_q2$aggY, seq(from = 0, to = 1, by = 0.2))
betaMeanQ <- tapply(res_q2$beta, findInterval(res_q2$aggY, incomeQuantile), mean)
betaMedQ <- tapply(res_q2$beta, findInterval(res_q2$aggY, incomeQuantile), median)

betaQuantile <- quantile(res_q2$absbeta, seq(from = 0, to = 1, by = 0.2))
incMeanQ <- tapply(res_q2$aggY, findInterval(res_q2$absbeta, betaQuantile), mean)
incMedQ <- tapply(res_q2$aggY, findInterval(res_q2$absbeta, betaQuantile), median)

result_q3 <- lm(alldlnc~0 + alldlninc +alldlnaggc)
betaAgg <- result_q3$coefficients[1]
phiAgg <- result_q3$coefficients[2]

final_q1all <- final_q1
betaAgg_all <- betaAgg
phiAgg_all <- phiAgg

