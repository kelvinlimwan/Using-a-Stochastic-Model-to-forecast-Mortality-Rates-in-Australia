library(StMoMo)
library(demography)
library(ggplot2)
library(fanplot)
library(tidyverse)
library(janitor)
library(showtext)

#General overview plot
death.rates <- subset.data.frame(read.csv("Mx.csv"), Age>=20 & Age <=89 & Year>=1960 & Year<=2016)
death.rates.1960 <- subset.data.frame(death.rates, Age >= 40 & Year==1960, select = c(Age, Total))
death.rates.1970 <- subset.data.frame(death.rates, Age >= 40 & Year==1970, select = c(Age, Total))
death.rates.1980 <- subset.data.frame(death.rates, Age >= 40 & Year==1980, select = c(Age, Total))
death.rates.1990 <- subset.data.frame(death.rates, Age >= 40 & Year==1990, select = c(Age, Total))
death.rates.2000 <- subset.data.frame(death.rates, Age >= 40 & Year==2000, select = c(Age, Total))
death.rates.2010 <- subset.data.frame(death.rates, Age >= 40 & Year==2010, select = c(Age, Total))

par(mfrow=c(1,1))
plot(death.rates.1960, xlab = "Age", ylab = "Central rate of death", main = "Death rates across decades", type = "l", col = "blue4", lwd = 2)
lines(death.rates.1970, col = "royalblue", lwd = 2)
lines(death.rates.1980, col = "dodgerblue", lwd = 2)
lines(death.rates.1990, col = "deepskyblue", lwd = 2)
lines(death.rates.2000, col = "skyblue1", lwd = 2)
lines(death.rates.2010, col = "lightskyblue1", lwd = 2)
legend("topleft", legend = c("1960", "1970", "1980", "1990", "2000", "2010"), col = c("blue4", "royalblue", "dodgerblue", "deepskyblue", "skyblue1", "lightskyblue1"), lty = 1, lwd = 2, title = "Legend", cex = 0.8)

#Heat map
showtext_auto()
aussie <- read_table(paste0("Deaths_1x1.txt"), skip = 2, na = ".") %>% clean_names()
aussie$age <- as.integer(recode(aussie$age,"110+"="110"))
aussie <- aussie %>% mutate(ratio = male / female,
                            deciles = cut(ratio, breaks = quantile(ratio, probs = seq(0, 1, 0.1), na.rm = TRUE)),
                            pct_diff = ((male - female) / (male + female))*100,
                            bin_ratio = ntile(ratio, 100))

p <- ggplot(aussie, aes(x = year, y = age, fill = ntile(total, 100)))
p_out <- p + geom_raster() +
  scale_fill_viridis_c(option = "A", direction = -1) +
  scale_x_continuous(breaks = seq(1960, 2020, by = 10)) +
  scale_y_continuous(breaks = seq(20, 90, by = 10)) +
  guides(fill = guide_legend(nrow = 1, title.position = "top", label.position = "bottom")) +
  labs(x = "Year", y = "Age", fill = "Death Rate Percentile",
       title = "Australian Mortality 1960-2016")+
  theme(legend.position = "top",
        legend.title = element_text(size = 8))
p_out

AU<-read.demogdata("Mx_1x1.txt","Exposures_1x1.txt",type="mortality",label="Australia")
AUtotal<-StMoMoData(AU,series="total")

#Fitting our models
## M1
LC1 <- lc()

LCfit1<-fit(LC1, data = AUtotal, ages.fit = 20:89,years.fit=1960:2016)
LCfit2<-fit(LC1, data = AUtotal, ages.fit = 20:89,years.fit=1980:2016)

plot(LCfit1)
plot(LCfit2)

logLik(LCfit1)

resLC1<-residuals(LCfit1)
plot(resLC1,type="signplot")
plot(resLC1,type="colourmap")

par(mfrow = c(2,2))
plot(LCfit1$ages, LCfit1$ax, xlab="Age", ylab="", main="M1: beta_1(x)", pch = 18, cex = 0.5)
lines(LCfit2$ages, LCfit2$ax, col = "blue")
plot(LCfit1$ages, LCfit1$bx, xlab="Age", ylab="", main="M1: beta_2(x)", pch = 18, cex = 0.5)
lines(LCfit2$ages, LCfit2$bx, col = "blue")
plot(LCfit1$years, LCfit1$kt, xlab="Year", ylab="", main="M1: kappa(t)", pch = 18, cex = 0.5)
lines(LCfit2$years, LCfit2$kt, col = "blue")

## M2
RH1<-rh(cohortAgeFun="NP",approxConst = TRUE)
RH2<-rh(approxConst = TRUE)
RHfit1<- fit(RH2, data = AUtotal, ages.fit = 20:89, years.fit=1960:2016, start.ax = LCfit1$ax, 
             start.bx = LCfit1$bx, start.kt = LCfit1$kt)
RHfit2<- fit(RH2, data = AUtotal, ages.fit = 20:89, years.fit=1980:2016, start.ax = LCfit2$ax, 
             start.bx = LCfit2$bx, start.kt = LCfit2$kt)

plot(RHfit1)
plot(RHfit2)

logLik(RHfit1)

resRH1<-residuals(RHfit1)
plot(resRH1,type="signplot")
plot(resRH1,type="colourmap")

par(mfrow = c(2,2))
plot(RHfit1$ages, RHfit1$ax, xlab="Age", ylab="", main = "M2: beta_1(x)", pch = 18, cex = 0.5)
lines(RHfit2$ages, RHfit2$ax, col = "blue")
plot(RHfit1$ages, RHfit1$bx, ylim = c(0.003, 0.026), xlab="Age", ylab="",
     main = "M2: beta_2(x)", pch = 18, cex = 0.5)
lines(RHfit2$ages, RHfit2$bx, col = "blue")
plot(RHfit1$years, RHfit1$kt, xlab="Year", ylab="", main="M2: kappa(t)", pch = 18, cex = 0.5)
lines(RHfit2$years, RHfit2$kt, col = "blue")
plot(RHfit1$cohorts, RHfit1$gc, ylim = c(-0.8, 0.3), xlab="Cohort", ylab="", main="M2: gamma(t-x)", pch = 18, cex = 0.5)
lines(RHfit2$cohorts, RHfit2$gc, col = "blue")


## M3
APC1 <- apc()
APCfit1<-fit(APC1, data = AUtotal, ages.fit= 20:89, years.fit = 1960:2016)
APCfit2<-fit(APC1, data = AUtotal, ages.fit= 20:89, years.fit = 1980:2016)

plot(APCfit1,parametricbx = FALSE)
plot(APCfit2,parametricbx = FALSE)

logLik(APCfit1)

resAPC<-residuals(APCfit1)
plot(resAPC,type="signplot")
plot(resAPC,type="colourmap")

par(mfrow = c(2,2))
plot(APCfit1$ages,APCfit1$ax,xlab="Age",ylab="",main="M3: beta_1(x)", pch = 18, cex = 0.5)
lines(APCfit2$ages,APCfit2$ax, col = "blue")
plot(APCfit1$years,APCfit1$kt,xlab="Year",ylab="",main="M3: kappa(t)", pch = 18, cex = 0.5)
lines(APCfit2$years,APCfit2$kt, col = "blue")
plot(APCfit1$cohorts,APCfit1$gc,xlab="Year of birth",ylab="",main="M3: gamma(t-x)", pch = 18, cex = 0.5)
lines(APCfit2$cohorts,APCfit2$gc, col = "blue")

## M5
CBD1 <- cbd(link="logit")
CBDfit1<-fit(CBD1, data = central2initial(AUtotal), ages.fit= 20:89, years.fit = 1960:2016)
CBDfit2<-fit(CBD1, data = central2initial(AUtotal), ages.fit= 20:89, years.fit = 1980:2016)

plot(CBDfit1, parametricbx = FALSE)
plot(CBDfit2, parametricbx = FALSE)

logLik(CBDfit1)

resCBD<-residuals(CBDfit1)
plot(resCBD,type="signplot")
plot(resCBD,type="colourmap")

plot(CBDfit1$years, CBDfit1$kt[1,], xlab="Year", ylab="", main="M5: kappa_1(t)", pch = 18, cex = 0.5)
lines(CBDfit2$years, CBDfit2$kt[1,], col = "blue")
plot(CBDfit1$years, CBDfit1$kt[2,], xlab="Year", ylab="", main="M5: kappa_2(t)", pch = 18, cex = 0.5)
lines(CBDfit2$years, CBDfit2$kt[2,], col = "blue")

## M6
M6.1 <- m6(link ="logit")
M6fit1<-fit(M6.1, data = central2initial(AUtotal), ages.fit= 20:89, years.fit = 1960:2016)
M6fit2<-fit(M6.1, data = central2initial(AUtotal), ages.fit= 20:89, years.fit = 1980:2016)

plot(M6fit1,parametricbx = FALSE)
plot(M6fit2,parametricbx = FALSE)

logLik(M6fit1)

resM6<-residuals(M6fit1)
plot(resM6,type="signplot")
plot(resM6,type="colourmap")## mapping residuals

par(mfrow=c(2,2))
plot(M6fit1$years, M6fit1$kt[1,], xlab="Year", ylab="", main="M6: kappa_1(t)", pch = 18, cex = 0.5)
lines(M6fit2$years, M6fit2$kt[1,], col = "blue")
plot(M6fit1$years, M6fit1$kt[2,], xlab="Year", ylab="", main="M6: kappa_2(t)", pch = 18, cex = 0.5)
lines(M6fit2$years,M6fit2$kt[2,], col = "blue")
plot(M6fit1$cohorts, M6fit1$gc, xlab="Year of birth", ylab="", main="M6: gamma(t-x)", pch = 18, cex = 0.5)
lines(M6fit2$cohorts,M6fit2$gc, col = "blue")

## M7
M7.1 <- m7(link="logit")
M7fit1<-fit(M7.1, data = central2initial(AUtotal), ages.fit= 20:89, years.fit = 1960:2016)
M7fit2<-fit(M7.1, data = central2initial(AUtotal), ages.fit= 20:89, years.fit = 1980:2016)

plot(M7fit1,parametricbx = FALSE)
plot(M7fit2,parametricbx = FALSE)

logLik(M7fit1)

resM7<-residuals(M7fit1)
plot(resM7,type="signplot")
plot(resM7,type="colourmap")

par(mfrow = c(2,2))
plot(M7fit1$years, M7fit1$kt[1,], xlab="Year", ylab="", main="M7: kappa_1(t)", pch = 18, cex = 0.5)
lines(M7fit2$years,M7fit2$kt[1,], col = "blue")
plot(M7fit1$years, M7fit1$kt[2,], xlab="Year", ylab="", main="M7:kappa_2(t)", pch = 18, cex = 0.5)
lines(M7fit2$years, M7fit2$kt[2,], col = "blue")
plot(M7fit1$years, M7fit1$kt[3,], xlab="Year", ylab="", main="M7:kappa_3(t)", pch = 18, cex = 0.5)
lines(M7fit2$years, M7fit2$kt[3,], col = "blue")
plot(M7fit1$cohorts, M7fit1$gc, xlab="Year of birth", ylab="", main="M7: gamma(t-x)", pch = 18, cex = 0.5)
lines(M7fit2$cohorts, M7fit2$gc, col = "blue")


## M8
M8.1 <- m8(xc=89)
M8fit1<-fit(M8.1, data = central2initial(AUtotal), ages.fit= 20:89, years.fit = 1960:2016)
M8fit2<-fit(M8.1, data = central2initial(AUtotal), ages.fit= 20:89, years.fit = 1980:2016)

plot(M8fit1,parametricbx = FALSE)
plot(M8fit2,parametricbx = FALSE)

logLik(M8fit1)

resM8<-residuals(M8fit1)
plot(resM8,type="signplot")
plot(resM8,type="colourmap")

par(mfrow = c(2,2))
plot(M8fit1$years,M8fit1$kt[1,], xlab="Year", ylab="", main="M8: kappa_1(t)", pch = 18, cex = 0.5)
lines(M8fit2$years, M8fit2$kt[1,], col = "blue")
plot(M8fit1$years,M8fit1$kt[2,], xlab="Year", ylab="", main="M8: kappa_2(t)", pch = 18, cex = 0.5)
lines(M8fit2$years, M8fit2$kt[2,], col = "blue")
plot(M8fit1$cohorts,M8fit1$gc, xlab="Year of birth", ylab="", main="M8: gamma(t-x)", pch = 18, cex = 0.5)
lines(M8fit2$cohorts, M8fit2$gc, col = "blue")


#Forecasts using M3
APCfor1 <- forecast(APCfit1,h=20,level=c(50,80,95))
plot(APCfor1,parametricbx = FALSE, col = c("skyblue3", "lightblue"))

death.rates.20 <- subset.data.frame(death.rates, Age == 20, select = c(Year, Total))
death.rates.40 <- subset.data.frame(death.rates, Age == 40, select = c(Year, Total))
death.rates.60 <- subset.data.frame(death.rates, Age == 60, select = c(Year, Total))
death.rates.80 <- subset.data.frame(death.rates, Age == 80, select = c(Year, Total))

#Plots for given ages forecast, also showing the fitted model for reference
par(mfrow = c(2,2))
plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["20", ], 
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 20", lwd = 2)
lines(APCfor1$years, APCfor1$rates["20", ], col = "red", lwd = 2)
points(death.rates.20, pch = 18, col = "blue", cex = 0.8) 

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["40", ], 
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 40", lwd = 2)
lines(APCfor1$years, APCfor1$rates["40", ], col = "red", lwd = 2)
points(death.rates.40, pch = 18, col = "blue", cex = 0.8) 

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["60", ], 
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 60", lwd = 2)
lines(APCfor1$years, APCfor1$rates["60", ], col = "red", lwd = 2)
points(death.rates.60, pch = 18, col = "blue", cex = 0.8) 

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["80", ], 
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 80", lwd = 2)
lines(APCfor1$years, APCfor1$rates["80", ], col = "red", lwd = 2)
points(death.rates.80, pch = 18, col = "blue", cex = 0.8) 

#Simulate data
APCsim<-simulate(APCfit1,nsim=10000,seed=1)

par(mfrow = c(2,2))
plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["20", ],
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 20", 
     xlim = c(1980, 2036), ylim = c(0.0002, 0.0014), lwd = 2)
points(death.rates.20, pch = 18, col = "blue", cex = 0.8) 
abline(h=c(0.00035, 0.00135), col="gray", lty=3)
fan(t(APCsim$rates["20",,]),start=2017,probs=c(2.5,10,25,75,90,97.5),n.fan=4,ln=NULL,fan.col = colorRampPalette(c("royalblue", "white")), med.ln = TRUE, med.col = "navy")

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["40", ],
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 40", 
     xlim = c(1980, 2036), ylim = c(0.0004, 0.0017), lwd = 2)
points(death.rates.40, pch = 18, col = "blue", cex = 0.8)
abline(h=c(0.0005, 0.0015), col="gray", lty=3)
fan(t(APCsim$rates["40",,]),start=2017,probs=c(2.5,10,25,75,90,97.5),n.fan=4,ln=NULL,fan.col = colorRampPalette(c("royalblue", "white")), med.ln = TRUE, med.col = "navy")

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["60", ],
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 60", 
     xlim = c(1980, 2036), ylim = c(0.0035, 0.014), lwd = 2)
points(death.rates.60, pch = 18, col = "blue", cex = 0.8) 
abline(h=seq(0.0035, 0.014, by=0.001), col="gray", lty=3)
fan(t(APCsim$rates["60",,]),start=2017,probs=c(2.5,10,25,75,90,97.5),n.fan=4,ln=NULL,fan.col = colorRampPalette(c("royalblue", "white")), med.ln = TRUE, med.col = "navy")

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["80", ],
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 80", 
     xlim = c(1980, 2036), ylim = c(0.02, 0.085), lwd = 2)
points(death.rates.80, pch = 18, col = "blue", cex = 0.8) 
abline(h=seq(0.02, 0.085, by=0.001), col="gray", lty=3)
fan(t(APCsim$rates["80",,]),start=2017,probs=c(2.5,10,25,75,90,97.5),n.fan=4,ln=NULL,fan.col = colorRampPalette(c("royalblue", "white")), med.ln = TRUE, med.col = "navy")

#Examine forecasts for ages to be used when simulating
death.rates.30 <- subset.data.frame(death.rates, Age == 30, select = c(Year, Total))
death.rates.50 <- subset.data.frame(death.rates, Age == 50, select = c(Year, Total))

par(mfrow = c(2,3))
plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["30", ], 
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 30", lwd = 2)
lines(APCfor1$years, APCfor1$rates["30", ], col = "red", lwd = 2)
abline(h=c(0.0004, 0.0009, 0.0014), col="gray", lty=3)
points(death.rates.30, pch = 18, col = "blue", cex = 0.8) 

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["40", ], 
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 40", lwd = 2)
lines(APCfor1$years, APCfor1$rates["40", ], col = "red", lwd = 2)
abline(h=seq(0.0005, 0.0025, by = 0.0005), col="gray", lty=3)
points(death.rates.40, pch = 18, col = "blue", cex = 0.8) 

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["50", ], 
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 50", lwd = 2)
lines(APCfor1$years, APCfor1$rates["50", ], col = "red", lwd = 2)
abline(h=seq(0.0015, 0.007, by = 0.0005), col="gray", lty=3)
points(death.rates.50, pch = 18, col = "blue", cex = 0.8)

#Include approximate PIs
plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["30", ],
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 30", 
     xlim = c(1980, 2036), ylim = c(0.0002, 0.0012), lwd = 2)
abline(h=c(0.0002, 0.0012), col="gray", lty=3)
points(death.rates.30, pch = 18, col = "blue", cex = 0.8) 
fan(t(APCsim$rates["30",,]),start=2017,probs=c(2.5,10,25,75,90,97.5),n.fan=4,ln=NULL,fan.col = colorRampPalette(c("royalblue", "white")), med.ln = TRUE, med.col = "navy")

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["40", ],
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 40", 
     xlim = c(1980, 2036), ylim = c(0.0004, 0.0018), lwd = 2)
abline(h=c(0.0005, 0.0015), col="gray", lty=3)
points(death.rates.40, pch = 18, col = "blue", cex = 0.8) 
fan(t(APCsim$rates["40",,]),start=2017,probs=c(2.5,10,25,75,90,97.5),n.fan=4,ln=NULL,fan.col = colorRampPalette(c("royalblue", "white")), med.ln = TRUE, med.col = "navy")

plot(c(APCfit1$years, APCfor1$years), cbind(APCfor1$fitted, APCfor1$rates)["50", ],
     type = "l", xlab = "Year", ylab = "Mortality rate", main = "Age 50", 
     xlim = c(1980, 2036), ylim = c(0.001, 0.0055), lwd = 2)
abline(h=seq(0.001, 0.005, by = 0.001), col="gray", lty=3)
points(death.rates.50, pch = 18, col = "blue", cex = 0.8) 
fan(t(APCsim$rates["50",,]),start=2017,probs=c(2.5,10,25,75,90,97.5),n.fan=4,ln=NULL,fan.col = colorRampPalette(c("royalblue", "white")), med.ln = TRUE, med.col = "navy")

#Simulation of policies in force
AUmale<-StMoMoData(AU,series="male")
AUfemale<-StMoMoData(AU,series="female")

APCfit.m<-fit(APC1, data = AUmale, ages.fit=20:89, years.fit=1960:2016)
APCfit.f<-fit(APC1, data = AUfemale, ages.fit=20:89, years.fit=1960:2016)

## Males
APCsim.m <- simulate(APCfit.m, nsim = 1000, seed = 1, h = 11)

mSim30 <- extractCohort(APCsim.m$rates, age = 30, period = 2017)
mSim40 <- extractCohort(APCsim.m$rates, age = 40, period = 2017)
mSim50 <- extractCohort(APCsim.m$rates, age = 50, period = 2017)

mc.30.m<-rowMeans(mSim30)
mc.40.m<-rowMeans(mSim40)
mc.50.m<-rowMeans(mSim50)

pc30.m <- exp(-mc.30.m)
pc40.m <- exp(-mc.40.m)
pc50.m <- exp(-mc.50.m)

tpx30m = vector(length=10)
tpx40m = vector(length=10)
tpx50m = vector(length=10)

for (i in 1:10) {
  tpx30m[i] = prod(pc30.m[1:i])
  tpx40m[i] = prod(pc40.m[1:i])
  tpx50m[i] = prod(pc50.m[1:i])
}

M30<-round(3000*tpx30m)
M40<-round(2000*tpx40m)
M50<-round(1000*tpx50m)

rbind(M30,M40,M50)

## Females
APCsim.f <- simulate(APCfit.f, nsim=1000, seed=1)

fSim30 <- extractCohort(APCsim.f$rates, age = 30, period = 2017)
fSim40 <- extractCohort(APCsim.f$rates, age = 40, period = 2017)
fSim50 <- extractCohort(APCsim.f$rates, age = 50, period = 2017)

mc.30.f <- rowMeans(fSim30)
mc.40.f <- rowMeans(fSim40)
mc.50.f <- rowMeans(fSim50)

pc30.f <- exp(-fSim30)
pc40.f <- exp(-fSim40)
pc50.f <- exp(-fSim50)

tpx30f = vector(length = 10)
tpx40f = vector(length = 10)
tpx50f = vector(length = 10)

for (i in 1:10) {
  tpx30f[i] = prod(pc30.f[1:i])
  tpx40f[i] = prod(pc40.f[1:i])
  tpx50f[i] = prod(pc50.f[1:i])
}

F30<-round(2000*tpx30f)
F40<-round(1000*tpx40f)
F50<-round(500*tpx50f)

rbind(F30,F40,F50)
