
### Fomite-mediated transmission model - including viral decay on fomites

#### Setup ####

####################################################################

rm(list=ls()) #clear all variables

# Open Packages
library(mvtnorm)
library(mc2d) 
library(dplyr)
set.seed(12346)

#10000 simulations - variability
ndvar(10000)


# Source functions --------------------------------------------------------

# Make sure it is in the same working directory as current code
# Or add file path before file name
source("C:\\Users\\Fomite Model QMRAFunctions Transport Decay 12.19.21 Final.R")


# Set controls/interventions -------------------------------------------------------

#Master controls for both aerosol and close contact modules
m.Event <-"cough"
m.room.exchange <- mcstoc(runif, type="V",min=2.0, max =2.0)
m.Humidity<-"high"

#Intervention controls for both aerosol and close contact modules
i.Clean1<-0
i.Clean2<-0
i.Clean3<-0
i.Clean4<-0
i.Clean5<-0
i.Clean6<-0
i.Clean7<-0
i.Clean8<-0
i.index.mask<-"none"
i.susceptible.mask<-"none"
i.HW<-"no"
i.Glove<-"no"



# Aerosol function call ---------------------------------------------------

aero.dose <- Aerofunc(Event = m.Event, room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)

aero.dose.clean <- select(aero.dose,f1h:f8h, f1hconcentration:f8hconcentration)


# Close contact function call ---------------------------------------------

dose50601m <- Dosefunc(Event = m.Event, Volume.Fraction="50-60", Distance=1.0, Vol.Frac.Dist.Name = "50601m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)
dose601001m <- Dosefunc(Event = m.Event, Volume.Fraction="60-100", Distance=1.0, Vol.Frac.Dist.Name = "601001m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)
dose1007501m <- Dosefunc(Event = m.Event, Volume.Fraction="100+", Distance=1.0, Vol.Frac.Dist.Name = "1007501m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)

dose5060_1m <- select(dose50601m,starts_with("f"))

dose5060_1m <- select(dose5060_1m, -fomite1c, -fomite1p, -fomite2c, -fomite2p, -fomite3c, -fomite3p, -fomite4c, -fomite4p, -fomite5c, -fomite5p,
                      -fomite6c, -fomite6p, -fomite7c, -fomite7p, -fomite8c, -fomite8p)

dose60100_1m <- select(dose601001m,starts_with("f"))
dose60100_1m <- select(dose60100_1m, -fomite1c, -fomite1p, -fomite2c, -fomite2p, -fomite3c, -fomite3p, -fomite4c, -fomite4p, -fomite5c, -fomite5p,
                       -fomite6c, -fomite6p, -fomite7c, -fomite7p, -fomite8c, -fomite8p)

dose100750_1m <- select(dose1007501m,starts_with("f"))
dose100750_1m <- select(dose100750_1m, -fomite1c, -fomite1p, -fomite2c, -fomite2p, -fomite3c, -fomite3p, -fomite4c, -fomite4p, -fomite5c, -fomite5p,
                        -fomite6c, -fomite6p, -fomite7c, -fomite7p, -fomite8c, -fomite8p)

# Combine aerosol and close contact doses ---------------------------------

#adding aerosol and aerosol fomite to close contact doses
dose1m <-cbind (dose5060_1m, dose60100_1m, dose100750_1m, aero.dose.clean)

dose1m <- mutate( dose1m, 
                  f1m1h = f50601m1h + f601001m1h + f1007501m1h + f1h,
                  f1m2h = f50601m2h + f601001m2h + f1007501m2h + f2h,
                  f1m3h = f50601m3h + f601001m3h + f1007501m3h + f3h,
                  f1m4h = f50601m4h + f601001m4h + f1007501m4h + f4h,
                  f1m5h = f50601m5h + f601001m5h + f1007501m5h + f5h,
                  f1m6h = f50601m6h + f601001m6h + f1007501m6h + f6h,
                  f1m7h = f50601m7h + f601001m7h + f1007501m7h + f7h,
                  f1m8h = f50601m8h + f601001m8h + f1007501m8h + f8h,
                  fconc1m1h = fconc50601m1h + fconc601001m1h + fconc1007501m1h + f1hconcentration,
                  fconc1m2h = fconc50601m2h + fconc601001m2h + fconc1007501m2h + f2hconcentration,
                  fconc1m3h = fconc50601m3h + fconc601001m3h + fconc1007501m3h + f3hconcentration,
                  fconc1m4h = fconc50601m4h + fconc601001m4h + fconc1007501m4h + f4hconcentration,
                  fconc1m5h = fconc50601m5h + fconc601001m5h + fconc1007501m5h + f5hconcentration,
                  fconc1m6h = fconc50601m6h + fconc601001m6h + fconc1007501m6h + f6hconcentration,
                  fconc1m7h = fconc50601m7h + fconc601001m7h + fconc1007501m7h + f7hconcentration,
                  fconc1m8h = fconc50601m8h + fconc601001m8h + fconc1007501m8h + f8hconcentration)
dose1m_risk <- select(dose1m, f1m1h:f1m8h)
fconc1m <- select(dose1m, fconc1m1h:fconc1m8h)

# Calculate risk ----------------------------------------------------------

#dose-response parameter updated to Julian et al 2020 preprint
krisk=0.00680

#vaccination for susceptible worker
mRNAvaccine <- mcstoc(runif, type="V",min=0.01, max =0.14)
reducedVEvaccine <- mcstoc(runif, type="V",min=0.20, max =0.36)
deltavaccine <- mcstoc(runif, type="V",min=0.12, max =0.28)
betavaccine<- mcstoc(runif, type="V",min=0.04, max =0.25)
veat<- mcstoc(rtriang, type="V",min=0.052, mode= 0.115, max =0.177)
mRNAvaccine <- unmc (mRNAvaccine,drop = TRUE)
mRNAvaccine <- as.data.frame(mRNAvaccine)
reducedVEvaccine<-unmc (reducedVEvaccine,drop = TRUE)
reducedVEvaccine<-as.data.frame(reducedVEvaccine)
deltavaccine<-unmc (deltavaccine,drop = TRUE)
deltavaccine<-as.data.frame(deltavaccine)
betavaccine<-unmc (betavaccine,drop = TRUE)
betavaccine<-as.data.frame(betavaccine)
veat<-unmc (veat,drop = TRUE)
veat<-as.data.frame(veat)

# Pull combined doses through dose response for data frame output of risk

#aerosol module risk (>3m)
riskaero.df = 1-exp(-krisk*aero.dose.clean)
aeroriskvaxx <-cbind(riskaero.df, mRNAvaccine, reducedVEvaccine, deltavaccine, betavaccine, veat)
aeroriskvaxxfull <-mutate(aeroriskvaxx, 
                             
                             mRNAvaxxf1h = mRNAvaccine * f1h,
                             mRNAvaxxf2h = mRNAvaccine * f2h,
                             mRNAvaxxf3h = mRNAvaccine * f3h,
                             mRNAvaxxf4h = mRNAvaccine * f4h,
                             reducedVEvaxxf1h = reducedVEvaccine * f1h,
                             reducedVEvaxxf2h = reducedVEvaccine * f2h,
                             reducedVEvaxxf3h = reducedVEvaccine * f3h,
                             reducedVEvaxxf4h = reducedVEvaccine * f4h,
                             deltavaxxf1h = deltavaccine * f1h,
                             deltavaxxf2h = deltavaccine * f2h,
                             deltavaxxf3h = deltavaccine * f3h,
                             deltavaxxf4h = deltavaccine * f4h,
                             betavaxxf1h = betavaccine * f1h,
                             betavaxxf2h = betavaccine * f2h,
                             betavaxxf3h = betavaccine * f3h,
                             veatvaxxf1h = veat*f1h)

riskaero.quant<-as.data.frame(t(apply(aeroriskvaxxfull, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
riskaero.mean<- as.data.frame(colMeans(aeroriskvaxxfull))
riskaero.comb.stats<- cbind(riskaero.quant, riskaero.mean)
write.csv(riskaero.comb.stats, "C:\\Users\\Desktop\\aeroriskfomitedecay90.csv", row.names=TRUE)

#close contact 1m risk
risk1m.df = 1-exp(-krisk*dose1m_risk)
risk1mvaxx <-cbind(risk1m.df, mRNAvaccine, reducedVEvaccine, deltavaccine, betavaccine, veat)
risk1mvaxxfull <-mutate(risk1mvaxx,
                          
                          mRNAvaxxf1m1h = mRNAvaccine * f1m1h,
                          mRNAvaxxf1m2h = mRNAvaccine * f1m2h,
                          mRNAvaxxf1m3h = mRNAvaccine * f1m3h,
                          reducedVEvaxxf1m1h = reducedVEvaccine * f1m1h,
                          reducedVEvaxxf1m2h = reducedVEvaccine * f1m2h,
                          reducedVEvaxxf1m3h = reducedVEvaccine * f1m3h,
                          deltavaxxf1m1h = deltavaccine * f1m1h,
                          deltavaxxf1m2h = deltavaccine * f1m2h,
                          deltavaxxf1m3h = deltavaccine * f1m3h,
                          betavaxxf1m1h = betavaccine * f1m1h,
                          betavaxxf1m2h = betavaccine * f1m2h,
                          betavaxxf1m3h = betavaccine * f1m3h,
                          veatvaxxf1m1h = veat * f1m1h)

risk1m.quant<-as.data.frame(t(apply(risk1mvaxxfull, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
risk1m.mean<- as.data.frame(colMeans(risk1mvaxxfull))
risk1m.comb.stats<- cbind(risk1m.quant, risk1m.mean)
write.csv(risk1m.comb.stats, "C:\\Users\\Desktop\\risk1mfomitedecay90.csv", row.names=TRUE)

#close contact 1m fomite concentration
fconc1mcomb <-as.data.frame(t(apply(fconc1m, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
write.csv(fconc1mcomb, "C:\\Users\\Desktop\\fomiteconcentration1mdecay90.csv", row.names=TRUE)
      







