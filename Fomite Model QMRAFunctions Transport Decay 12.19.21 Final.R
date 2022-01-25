#Functions for fomite QMRA model 
#Incorporating viral decay on fomite

################################################################
#Defining function for aerosol module

# Aerosol Function --------------------------------------------------------

Aerofunc <- function(Event, room.exchange, Clean1, Clean2, Clean3, Clean4, Clean5, Clean6, Clean7, Clean8, Humidity, Infected.Mask, HW, Glove, Susceptible.Mask){
  
  #Hydraulic Diameter - in (meters)
  facility.length <- 10
  facility.width <- 10
  facility.height <- 10
  
  #Inner volume of facility units is m3 - looks good
  facility.volume <- (facility.height * facility.length * facility.width)
  #Inner facility area m^2
  facility.area = facility.length*facility.width
  
  #Used to calculate the Reynolds number later on (essentially this is a wind tunnel)
  dh <- ( 4 * (facility.length * facility.width) / (2 * (facility.length + facility.width)))
  
  #### Virus Calculations ###
  
  ###viral concentration in saliva, pfu/ml
  log.C.virus <- mcstoc(rtriang, type="V", min=6.1, mode=6.8, max=7.4)
  
  #2 times more viral particles
  #variant2.log.C.virus<- mcstoc(rtriang, type="V", min=6.38, mode=7.1, max=7.67)
  ###variant 10 times more viral particles, viral concentration in saliva, pfu/ml
  #variant10.log.C.virus <- mcstoc(rtriang, type="V", min=7.1, mode=7.8, max=8.4)
  ###variant 100 times more viral particles, viral concentration in saliva, pfu/ml
  #variant100.log.C.virus <- mcstoc(rtriang, type="V", min=8.1, mode=8.8, max=9.4)
  ###variant 1000 times (delta) more viral particles, viral concentration in saliva, pfu/ml
  #variant1000.log.C.virus <- mcstoc(rtriang, type="V", min=9.1, mode=9.8, max=10.4)
 
  
  ###volume of saliva expelled with a cough ml
  saliva.vol <- mcstoc(runif, type="V", min=0.0396, max=0.0484)
  
  ###viral concentration in saliva, pfu/ml reduced for vaccination among infected worker 2.8 fold reduction
  #VR.virus <- mcstoc(rtriang, type="V", min=5.7, mode=6.4, max=7.0)
 
  #viral concentration in saliva, pfu/ml reduced for vaccination among infected worker 4.5 fold reduction
  #VM.virus <- mcstoc(rtriang, type="V", min=5.4, mode=6.1, max=6.7)
  
  C.virus <- log.C.virus
  #C.virus <- variant10.log.C.virus
  #C.virus <- variant100.log.C.virus
  #C.virus <- variant1000.log.C.virus
  #C.virus <- variant2.log.C.virus
  #C.virus <- VR.virus
  #C.virus <- VM.virus
  
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7126899/
  ### Note: dpa is in cm rather than um for units to work later on; this range is representative of 2-49um (airborne) COUGH
  ### Note: dpa is in cm rather than um for units to work later on; this range is representative of 0.3um-20um (airborne) BREATH
  dpacough <- mcstoc(rtriang, type = "V", min = 0.0002, mode = 0.0006, max = 0.0049)
  dpabreath <- mcstoc(rtriang, type = "V", min = 0.00003, mode = 0.00008, max = 0.0020)
  
  ###Switch for droplet size by event
  ifelse(Event =="cough", dpa <- dpacough, dpa <- dpabreath)
  
  ##number of people infected can change, right now fixed to 1
  people=1
  
  ##total fraction of volume composed of droplet sizes 2-45 (ml/cough) --from paper reference above COUGH event
  v2.50 <- mcstoc(rtriang, type = "V", min =0.0000013938 , mode =0.000002349, max =0.0000026395)
  
  ##total fraction of volume composed of droplet sizes <0.6 to 2.2 (ml/breathing) Papineni paper -- BREATHING event
  vbreath <- mcstoc(runif, type = "V", min =0.000000000105743, max =0.000000000289581)
  
  ###Switch for volume fraction by event
  ifelse(Event =="cough", VF <- v2.50, VF <- vbreath)
  
  #breathing rate: #breaths per minute per person range 16-20 breaths/min, Fabian et al., Origin of Exhaled Breath Particles from Healthy and Human Rhinovirus-Infected 
  # converted to breaths per hour per person; this range represents both high and low producers
  breathrate <- mcstoc(runif, type="V", min=960, max=1200)
  
  #cough frequency: #coughs per hour per person
  coughfreq <- mcstoc(rempiricalD,type="V",values=c(0,1,2,3,4,5,6,7,8,9,10),prob=c(1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11 ))
  
  ###Switch for event frequency
  ifelse(Event =="cough", eventfreq <- coughfreq, eventfreq <- breathrate)

  #### Mask intervention switches for infected individual ####
  
  pwsm <- as.numeric(ifelse(Infected.Mask == "surgical", pwsm <- 1, pwsm <- 0))
  pwcm <- as.numeric(ifelse(Infected.Mask == "cloth", pwcm <- 1, pwcm <- 0))
  pwoN95 <- as.numeric(ifelse(Infected.Mask == "oN95", pwoN95 <- 1, pwoN95 <- 0))
  #pwdouble <- as.numeric(ifelse(Infected.Mask == "double", pwdouble <- 1, pwdouble <- 0))
  
  ### These are all log-reductions for each mask type
  
  ### Mask efficacies reported as log reductions
  
  ## (Source protection) surgical mask efficacy log reduction
  s.mask <- mcstoc(runif,type="V",min=0.387216,max=0.568636)
  
  ## (Source protection) cloth mask efficacy log reduction 
  c.mask <- mcstoc(runif,type="V",min=0.309804,max=0.619789)
  
  ## (Source protection) fit N95 efficacy log reduction
  oN95.mask <- mcstoc(runif, type="V", min=1.886057, max=2.154902)
  
  ## (Source protection) double masking (surgical mask then cloth) log reduction
  #double.mask <- mcstoc(runif, type="V", min=0.274, max=1.89)
  
  ### PPE mask efficacies reported as percent reduction for mask module
  ### Note: these are different parameters from the source protection
  
  s.mask.p <- mcstoc(runif, type="V", min=0.37, max=0.998)
  c.mask.p <- mcstoc(runif, type="V", min=0.17, max=0.887)
  oN95.mask.p <- 0.79
  #double.mask.p <- mcstoc(runif, type="V", min=0.40, max=0.9685)
  #### Shedding ####
  
  #distance from hand to face (m)
  disHF <- mcstoc(runif, type="V", min=0.05, max=0.1)
 
  #right angle of cone (degrees)
  angle <- mcstoc(runif, type="V", min=27.5, max=35)
  
  #converting degrees to radians
  angle_rad <- (angle*pi)/180
  
  radius <- disHF * tan(angle_rad)
  
  Ainoc <- pi * (radius^2)
  
  #Area of two hands (palm only) in m^2 that touches fomite surface, 245 cm2 for one palm
  hand.sa <- 0.0490
  
   ### E.virus is total viral shedding (PFU/hr) with mask switches added:
  ##Unit check: v2.50 (Volume.Fraction, is ml/cough); C.virus is PFU/ml; cough.freq is coughs per hour per person
  E.virus <- (VF * eventfreq * people)*(10^(C.virus-(pwsm*s.mask)-(pwcm*c.mask)-(pwoN95*oN95.mask)))
  
  #can do either hand.sa or Ainoc; PFU/m2
  C_hand_inf <- (saliva.vol*(10^(C.virus-(pwsm*s.mask)-(pwcm*c.mask)-(pwoN95*oN95.mask))))/Ainoc
  
  ### converting from PFU/hr to PFU/s):
  sigma <- E.virus * (1/3600) 
  
  #### Surface Area Calculations ####
  
  #number of individual cartons contaminated (note each of these is only one side of the carton)
  contamcartons <- 0
  contamcartons<- as.numeric(ifelse(coughfreq == 0, contamcartons <- 0, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 1, contamcartons <- 1, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 2, contamcartons <- 2, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 3, contamcartons <- 3, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 4, contamcartons <- 4, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 5, contamcartons <- 5, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 6, contamcartons <- 6, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 7, contamcartons <- 7, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 8, contamcartons <- 8, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 9, contamcartons <- 9, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 10, contamcartons <- 10, contamcartons <- contamcartons))
  
  #S.us surface area of facility (this stays the same from before)
  S.us<-facility.area
  
  #SA of top of carton to be exposed by worker 1 and 2;range from small to larger cartons (m2)
  cartontop.SA<-mcstoc(runif, type = "V", min=0.106, max =0.116)
  
  #composite surface area of tops of cartons (m2) of contaminated individual cartons (event 1 and 2)
  contamcompcartSA<-cartontop.SA*contamcartons
  
  #SA for single carton (m2)
  cartonSA<-mcstoc(runif, type = "V", min=0.41, max =0.54)
  
  #total number of individual cartons processed in one hour (36 to 54 per pallet with 4 pallets per hour)
  carton.number<-mcstoc(runif, type = "V", min=144, max =216)
  
  #composite surface area of all cartons (event 1 and 2) processed in one hour
  compcartonSA<-cartonSA*carton.number
  
   #number of sides of pallets contaminated by the second worker (pallet has 6 sides in total, only 4 are the available sides; 4 sides with 4 pallets = 16 side options)
  contampallet <- 0
  contampallet<- as.numeric(ifelse(coughfreq == 0, contampallet <- 0, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 1, contampallet <- 1, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 2, contampallet <- 2, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 3, contampallet <- 3, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 4, contampallet <- 4, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 5, contampallet <- 5, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 6, contampallet <- 6, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 7, contampallet <- 7, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 8, contampallet <- 8, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 9, contampallet <- 9, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 10, contampallet <- 10, contampallet <- contampallet))
 
   #number of pallets processed per hour
  pallets.processed<- 4
  
   #SA of the one side of plastic wrap (m2) this is representative of 10x10 or 9X4 or 9x6 cartons stacked
  singlesideplasticwrap.SA<-mcstoc(runif, type = "V", min=4.2, max =6.97)
  
  #SA of contaminated composite plastic wrap sides (factoring in the number of pallets sides contaminated by the second worker)
  contamcompplastic.SA<-singlesideplasticwrap.SA*contampallet
  
  #SA of the total plastic wrap (m2) this is representative of 10x10 or 9X4 or 9x6 cartons stacked
  plasticwrap.SA<-mcstoc(runif, type = "V", min=25.2, max =41.82)
 
  #SA of the total plastic wrap (m2) across the 4 pallets processed in 1h
  compplasticwrap.SA<-plasticwrap.SA*pallets.processed
  
  
  #### Air Dynamics and Temperature####
  
  #Temperature Distribution - This can be the user input variable, with selections from 60-100F by increments of 5F
  Temp.F <- 38
  
  #Absolute Temperature in K - This automatically converts the temperature input above to units of Kelvin
  Temp <- (Temp.F + 459.67)*(5/9)
  
  ### Updating the values to be in m2 / s 
  #Kinematic Viscosity of Air at 70F (m2/s) -- Corrected; based off of below, may update
  kv <- 0
  kv <- as.numeric(ifelse(Temp.F == 38, kv <- 1.357E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 40, kv <- 1.367E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 55, kv <- 1.441E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 60, kv <- 1.466E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 65, kv <- 1.491E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 70, kv <- 1.516E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 75, kv <- 1.5415E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 80, kv <- 1.567E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 85, kv <- 1.593E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 90, kv <- 1.619E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 95, kv <- 1.645E-5, kv <- kv))
  kv <- as.numeric(ifelse(Temp.F == 100, kv <- 1.671E-5, kv <- kv))
  
  #Dynamic Viscosity of Air at 70F ([μPa s], [N s/m2 *10-6]) --- Corrected, based off of below, may update
  dv <- 0
  dv <- as.numeric(ifelse(Temp.F == 38, dv <- 1.732E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 40, dv <- 1.737E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 55, dv <- 1.7781E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 60, dv <- 1.791E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 65, dv <- 1.8045E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 70, dv <- 1.818E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 75, dv <- 1.8315E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 80, dv <- 1.845E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 85, dv <- 1.858E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 90, dv <- 1.871E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 95, dv <- 1.884E-5, dv <- dv))
  dv <- as.numeric(ifelse(Temp.F == 100, dv <- 1.897E-5, dv <- dv))
  
  #Air density at Sea Level (kg/m^3) -- Correct. Consider Salinas valley
  air.den <- 1.225
  
  #Boltzmann Constant (J/K)
  k <- 1.38 * (10^(-23))
  
  #Mean Free Path of Air at Sea Level (cm)      
  lambda <- 3.4 * (10^(-6))
  
  #### Viral decay, settling velocity, and air exchange calculations ####
  
  ###confirmed this should be in seconds; converted to s later, here it is currently in hours
  #v.decay.lh <- mcstoc(runif, type = "V", max = 0.091083729, min = 0.046645167)
  #v.decay.hh <- mcstoc(runif, type = "V", max = 0.15267559, min = 0.067624115)
  ### viral decay per hour, updated to van Doremalen 2020 paper for aerosols (converted from TCID50 to PFU then
  ### calculated the half life and then converted to viral decay rate from the half life)
  ### also consistent with temp of 70 and RH of 40%, dhs.gov/science-and-technology/sars-airborne-calculator
  ### v.decay.hh updated to be 55F with RH 40% to reflect the colder temps
  v.decay.lh <- 0.614
  #v.decay.hh <- 0.614
  #v.decay.hh <- 0.1447
  v.decay.hh <- 0
  
  ifelse(Humidity =="high", v.decay <- v.decay.hh, v.decay <- v.decay.lh)
  
  #Settling Velocity (m/hr) per MERS Paper converted to m/s
  ###Changing to m instead of cm for dpa here
  vs <- (((0.108 * (dpa*10**-2)^2 * (1 + ((0.166) / (dpa*10**-2)) ))) / 3600)
  
  #Room ventalation (m3/s); converted from room exchanges per hour to total room volume replaced per hour 
  room.exchange <- room.exchange
  Q <- (room.exchange*facility.volume)/3600
  
  #Ac is the cross sectional area of the building m2
  Ac=facility.length * facility.width
  
  #U.avg is the mean air flow speed m/s 
  U.avg <- Q/Ac
  
  #### Air transport model: Friction velocity calculations in the market ####
  #Reynold's Number
  Re <- (U.avg)*(dh / kv)
  
  #Skin Friction Coefficient
  Cf <- 0.027 / ((Re)^(1/7))
  
  #Shear Stress
  Tw <- Cf * (0.5) * air.den * ((U.avg)^2)
  
  #Friction Velocity
  ###Had to update to being a vector for r.plus math to work
  u.star <- U.avg * sqrt((Cf / 2))
  
  ############################Deposition calculations in the market
  #R plus
  ###dpa needs to be in m here instead of cm -> X*E-2
  r.plus <- (((dpa*(10^-2)) * u.star) / (2 * kv))
  
  #Cunningham Correction Factor
  Cc <- 1 + ((2 * lambda) / dpa) * (1.257 + 0.4*exp((-1.1 * dpa) / (2 * lambda)))
  
  #Brownian Diffusivity
  ### dpa needs to be in m again here
  D <- (k * Temp * Cc) / (3 * pi * dv * (dpa * 10**-2))
  
  #Schmidt Number
  Sc <- kv / D
  
  #A for Integral Equation
  A <- (0.5) * log( ((10.92 * (Sc^(-1/3)) + 4.3)^(3)) / ((Sc^(-1)) + 0.0609) ) + 
    sqrt(3) * atan((8.6 - 10.92 * Sc^(-1/3)) / (sqrt(3) * 10.92 * Sc**(-1/3)))
  
  #B for Integral Equation
  B <- (0.5) * log( ((10.92 * Sc^(-1/3) + r.plus)^(3)) / (Sc^(-1) + (7.669 * 10^(-4))) * (r.plus)^(3)) + 
    sqrt(3) * atan( (2*r.plus - 10.92 * Sc^(-1/3)) / (sqrt(3) * 10.92 * Sc^(-1/3)))
  
  #Integral Equation
  I <- 3.64 * (Sc^(2/3)) * (A - B) + 39
  
  #### Deposition velocities for each plane ####
  
  #Deposition velocity on the upward-facing horizontal surface (m/s)
  v.us <- (vs / (1 - exp(-((vs * I) / u.star))))
  
  #### V.loss and standard, no intervention air concentrations + calculations ####
  ### Confirmed Ct values are in PFU/m^3
  ### Note: while the market paper says the timesteps are in hours, the units for delta t need to be s
  
  ### V loss is in m^3/s (removed settling velocity from v.loss term)
  V.loss <- Q + ((v.decay/3600) * facility.volume)
  
  ### In EQ 4, c(0) = 0, so we can assume that the concentration at time step 0 is 0 PFU / m^3
  Ct0 <- 0
  
  ### Assuming that all of our time steps will have a time difference of one hour (converted to seconds for the units)
  delta.t <- 3600
  ### So the Concentration inside the Marker at hour t = 1 is (PFU/m^3):
  Ct1 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  
  #units for Ctcarry over PFU/m3
  Ctcarry <- (Ct0 * exp( (-V.loss * delta.t) / facility.volume ))
  #units for total Ct1 PFU/m3
  totalCt1<- Ct1+Ctcarry
  fallout1<- totalCt1*v.us*S.us*delta.t
  ##fallout is in PFU units
  ### inhaleCt is in PFU/m3
  inhaleCt1 <- totalCt1-(fallout1/facility.volume)
  
  Ct2 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  
  Ctcarry2 <- ((totalCt1-(fallout1/facility.volume)) * exp( (-V.loss * delta.t) / facility.volume ))
  totalCt2<- Ct2+Ctcarry2
  
  ##fallout is in PFU units
  fallout2<- totalCt2*v.us*S.us*delta.t
  
  ### inhaleCt is in PFU/m3
  inhaleCt2 <- totalCt2-(fallout2/facility.volume)
  
  Ct3 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  
  Ctcarry3 <- ((totalCt2-(fallout2/facility.volume)) * exp( (-V.loss * delta.t) / facility.volume ))
  totalCt3<- Ct3+Ctcarry3
  
  ##fallout is in PFU units
  fallout3<- totalCt3*v.us*S.us*delta.t
  
  ### inhaleCt is in PFU/m3
  inhaleCt3 <- totalCt3-(fallout3/facility.volume)
  
  Ct4 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  
  Ctcarry4 <- ((totalCt3-(fallout3/facility.volume)) * exp( (-V.loss * delta.t) / facility.volume ))
  totalCt4<- Ct4+Ctcarry4
  
  ##fallout is in PFU units
  fallout4<- totalCt4*v.us*S.us*delta.t
  
  ### inhaleCt is in PFU/m3
  inhaleCt4 <- totalCt4-(fallout4/facility.volume)
  
  Ct5 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  
  Ctcarry5 <- ((totalCt4-(fallout4/facility.volume)) * exp( (-V.loss * delta.t) / facility.volume ))
  totalCt5<- Ct5+Ctcarry5
  
  ##fallout is in PFU units
  fallout5<- totalCt5*v.us*S.us*delta.t
  
  ### inhaleCt is in PFU/m3
  inhaleCt5 <- totalCt5-(fallout5/facility.volume)
  
  Ct6 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  
  Ctcarry6 <- ((totalCt5-(fallout5/facility.volume)) * exp( (-V.loss * delta.t) / facility.volume ))
  totalCt6<- Ct6+Ctcarry6
  
  ##fallout is in PFU units
  fallout6<- totalCt6*v.us*S.us*delta.t
  
  ### inhaleCt is in PFU/m3
  inhaleCt6 <- totalCt6-(fallout6/facility.volume)
  
  Ct7 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  
  Ctcarry7 <- ((totalCt6-(fallout6/facility.volume)) * exp( (-V.loss * delta.t) / facility.volume ))
  totalCt7<- Ct7+Ctcarry7
  
  ##fallout is in PFU units
  fallout7<- totalCt7*v.us*S.us*delta.t
  
  ### inhaleCt is in PFU/m3
  inhaleCt7 <- totalCt7-(fallout7/facility.volume)
  
  Ct8 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  
  Ctcarry8 <- ((totalCt7-(fallout7/facility.volume)) * exp( (-V.loss * delta.t) / facility.volume ))
  totalCt8<- Ct8+Ctcarry8
  
  ##fallout is in PFU units
  fallout8<- totalCt8*v.us*S.us*delta.t
  
  ### inhaleCt is in PFU/m3
  inhaleCt8 <- totalCt8-(fallout8/facility.volume)
  
  ############# Adjusting the Fomite Calculations & Adding Surface Cleaning Switches ####
  
  # Percent Reduction for cleaners being used - ranging from 99% - 99.9999% reduction
  sc.eff.p <- mcstoc(runif, type = "V", min = 0.999, max = 0.999)
  
  ###surface areas of stainless steel tables m^2; 0.5m by 0.5m
  #this represents the reach of an individual's arm
  sstable.sa <- 0.25
  
  #Area of 3 finger tips touching the surface. 2.1 cm2 for 3 fingers from one hand; * 2 hands = 4.2 cm2
  #converted to m2; Bouwknegt et al 2015
  fingers.sa <- 0.00042
  
  #Area of two hands (palm only) in m^2 that touches fomite surface, 245 cm2 for one palm
  hand.sa <- 0.0490
  
  #fomite stainless steel viral decay (from NEJM van Doremalen et al 2020); note: converted from TCID50 to PFU to viral decay rate
  #updated to be for plastic
  #fomite.decay <- 0.099139
  
  # New fomite calculations, final unit here is PFU
  # new fomite decay calculated based off of plastic pitol et al; converted to days
  fd <- mcstoc(runif, type = "V", min = 0.14, max = 0.22)
  # new decay time in days
  decaytime <- mcstoc(runif, type = "V", min = 90, max = 90)
  
  ### Fallout in PFU
  ### Generating composite contaminated cartons
  #Fomite 1 and fomite 2 units should be PFU/m2
  #Hour 1
  fomite1c <- (((fallout1/contamcompcartSA) * (hand.sa/compcartonSA)) * (1 - sc.eff.p * Clean1))*(exp(-fd*decaytime))
  ### Generating composite contaminated plastic wrap; inhaleCt1 is in PFU/m3
  fomite1p <- ((((inhaleCt1 *facility.volume)/contamcompplastic.SA) * (hand.sa/compplasticwrap.SA)) * (1 - sc.eff.p * Clean1))*(exp(-fd*decaytime)) 
  #Hour 2
  fomite2c <- (((fallout2/contamcompcartSA) * (hand.sa/compcartonSA)) * (1 - sc.eff.p * Clean2))
  ### Generating composite contaminated plastic wrap; inhaleCt1 is in PFU/m3
  fomite2p <- ((((inhaleCt2 *facility.volume)/contamcompplastic.SA) * (hand.sa/compplasticwrap.SA)) * (1 - sc.eff.p * Clean2))
  #Hour 3
  fomite3c <- (((fallout3/contamcompcartSA) * (hand.sa/compcartonSA)) * (1 - sc.eff.p * Clean3))
  ### Generating composite contaminated plastic wrap; inhaleCt1 is in PFU/m3
  fomite3p <- ((((inhaleCt3 *facility.volume)/contamcompplastic.SA) * (hand.sa/compplasticwrap.SA)) * (1 - sc.eff.p * Clean3))
  #Hour 4
  fomite4c <- (((fallout4/contamcompcartSA) * (hand.sa/compcartonSA)) * (1 - sc.eff.p * Clean4))
  ### Generating composite contaminated plastic wrap; inhaleCt1 is in PFU/m3
  fomite4p <- ((((inhaleCt4 *facility.volume)/contamcompplastic.SA) * (hand.sa/compplasticwrap.SA)) * (1 - sc.eff.p * Clean4))
  #Hour 5
  fomite5c <- (((fallout5/contamcompcartSA) * (hand.sa/compcartonSA)) * (1 - sc.eff.p * Clean5))
  ### Generating composite contaminated plastic wrap; inhaleCt1 is in PFU/m3
  fomite5p <- ((((inhaleCt5 *facility.volume)/contamcompplastic.SA) * (hand.sa/compplasticwrap.SA)) * (1 - sc.eff.p * Clean5))
  #Hour 6
  fomite6c <- (((fallout6/contamcompcartSA) * (hand.sa/compcartonSA)) * (1 - sc.eff.p * Clean6))
  ### Generating composite contaminated plastic wrap; inhaleCt1 is in PFU/m3
  fomite6p <- ((((inhaleCt6 *facility.volume)/contamcompplastic.SA) * (hand.sa/compplasticwrap.SA)) * (1 - sc.eff.p * Clean6))
  #Hour 7
  fomite7c <- (((fallout7/contamcompcartSA) * (hand.sa/compcartonSA)) * (1 - sc.eff.p * Clean7))
  ### Generating composite contaminated plastic wrap; inhaleCt1 is in PFU/m3
  fomite7p <- ((((inhaleCt7 *facility.volume)/contamcompplastic.SA) * (hand.sa/compplasticwrap.SA)) * (1 - sc.eff.p * Clean7))
  #Hour 8
  fomite8c <- (((fallout8/contamcompcartSA) * (hand.sa/compcartonSA)) * (1 - sc.eff.p * Clean8))
  ### Generating composite contaminated plastic wrap; inhaleCt1 is in PFU/m3
  fomite8p <- ((((inhaleCt8 *facility.volume)/contamcompplastic.SA) * (hand.sa/compplasticwrap.SA)) * (1 - sc.eff.p * Clean8))
  
  
  #fomite3 <- (((fomite2 + ((fallout3) * (sstable.sa/facility.area) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean3))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite4 <- (((fomite3 + ((fallout4) * (sstable.sa/facility.area) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean4))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite5 <- (((fomite4 + ((fallout5) * (sstable.sa/facility.area) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean5))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite6 <- (((fomite5 + ((fallout6) * (sstable.sa/facility.area) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean6))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite7 <- (((fomite6 + ((fallout7) * (sstable.sa/facility.area) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean7))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite8 <- (((fomite7 + ((fallout8) * (sstable.sa/facility.area) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean8))/fomite.decay)*(1-exp(-fomite.decay*1))
  
  #Hand Hygiene Module linking in with fomite/surface cleaning code
  
  ## hand wash efficacy log removal (2 log10 PFU removal converted to %reduction)
  hw.eff <- mcstoc(rtriang, type="V", min=0.99, mode=0.99, max=.99)
  
  # Probability of harvester handwashing
  prob.hw <- ifelse(HW == "yes", prob.hw <- 1, prob.hw <- 0)
  
  # Frequency of washing hands per hour (converted to minutes to work with math)
  Freq.hw <- mcstoc(rempiricalD,type="V",values=c(1,1,1,1,1),prob=c(0.2,0.2,0.2,0.2,0.2))
  
  # Probability of using gloves
  prob.glove <- ifelse(Glove == "yes", prob.glove <- 1, prob.glove <- 0)
  
  # Frequency of glove changes per hour (converted to minutes to work with math)
  freq.glove <- mcstoc(rempiricalD,type="V",values=c(1,1,1,1,1),prob=c(0.2,0.2,0.2,0.2,0.2))
  
  #fraction of pathogen transfer from hand to glove (4) during gloving process
  #f24<-mcstoc(runif,type="V",min=0,max=0.444)
  f24<-mcstoc(runif, type="V", min=1, max=1)
  
  #Area of fomite in m^2 that touches hand surface
  #Assuming to be equal to the area of hand that touches the fomite
  #Per the influenza paper
  #table.sa<- 0.0245
  
  #Fraction of pathogens transferred from fomite (1) to hand (2) 
  #NOTE: f12.hh is in a high humidity environment (40% - 65%)
  #NOTE: f12.lh is in a low humidity environment (15% - 32%)
  #Per https://dx.doi.org/10.1128%2FAEM.01030-13 Using MS2
  #updated to be for plastic at low humidity
  f12.hh <- mcstoc(rtriang, type = "V", min =0.067, mode = 0.217, max =0.367)
  #f12.hh <- mcstoc(rnorm, type = "V", mean = 0.374, sd = 0.16)
  f12.lh <- mcstoc(rtriang, type = "V", min = 0, mode=0.069, max= 0.158)
  
  #Humidity switch
  f12 <-ifelse(Humidity =="high", f12 <- f12.hh, f12<- f12.lh)
  
  #Fraction of pathogens transferred from hand (2) to fomite (1)
  #updated parameter based off of Alicia's EID paper
  f21 <- 0.025
  
  #Fraction of pathogens transferred from hand (2) to face (3) 
  #Assumed to be 35% per touch in the influenza paper
  #Making distribution based on table in influenza paper for lips
  #Updated based on Julian et al MS2 paper (same as above for transfer efficiencies)
  #NOTE: This was not RH dependent
  f23 <- mcstoc(rtriang, type = "V", min =0.137, mode = 0.20, max =0.263)
  
  #Frequency of contacts from hand to individual cartons per hour converted to (per minute)
  freq.hs <- (carton.number/60)
  
  #Frequency of contacts from hand to plastic wrap; distribution in times per hour converted to (per minute)
  freq.hs.plastic <- mcstoc(runif, type = "V", min=4/60, max =20/60)
  
  #Frequency of contacts from hand to face (per min)
  freq.hf <- 0.8
  
  #updated viral inactivation rate on hands per Nicas et al (see below): based off of influenza A. Note: this is already in units min^-1
  #https://onlinelibrary-wiley-com.proxy.library.emory.edu/doi/full/10.1111/j.1539-6924.2009.01253.x
  v.decay.skin<- mcstoc(runif, type = "V", min=0.92, max =1.47)
  v.decay.min <- v.decay.skin
  
  #converted from cm2 to m2 (all from https://aem.asm.org/content/aem/84/18/e00709-18.full.pdf)
  eyes.sa<-mcstoc(runif, type = "V", min =1E-5 , max =2E-4 )
  nose.sa<-mcstoc(runif, type = "V", min = 1E-5, max =1E-3 )
  mouth.sa <-mcstoc(runif, type = "V", min = 1E-4, max =4.1E-3 )
  
  face.sa<-eyes.sa+nose.sa+mouth.sa
  
   #Viral decay in minutes based on hand surface (collapsed this into one line)
  v.decay.min.hand <- v.decay.min + freq.hs * f21 + freq.hf * f23
 
  #Viral decay in minutes based on hand surface (collapsed this into one line)
  v.decay.infec <- v.decay.min + freq.hs * f12
  
  
  #Concentration of pathogens on hand per m^2
  #Assume at time 0, C(t) = 0
  time.min = 60
  T = mcstoc(runif, type = "V", min = 60, max = 60)
  
  C_surf_inf <- (((freq.hs * (Ainoc/hand.sa) * C_hand_inf * f21) / (v.decay.infec*time.min)) * 
                    (time.min + (exp(-v.decay.infec * time.min)-1)/v.decay.infec))*((1-hw.eff * prob.hw)^Freq.hw)
  
  ### new fomite here based off of hand contamination (PFU/m2)
  Surf_infect_cartons <- C_surf_inf * (hand.sa/compcartonSA) * (1 - sc.eff.p * Clean1)
  
  #Building in handwashing (freq efficacy prob of hw) and glove use/changes into the same step at hour 1.
  # concentration on hand from events 1-2
  C.hand.1c <- ((freq.hs * fomite1c * f12) / (v.decay.min.hand)) * 
    (1-exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
         
  #fraction of dose infectious, plan to run 0.01, 0.1, and 1
  infectious<- 1
  
  DT.hh1c <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.1c * f23 * T * infectious
  
  C.hand.1p <- ((freq.hs.plastic * fomite1p * f12) / (v.decay.min.hand)) * 
                                        (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
        
  DT.hh1p <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.1p * f23 * T * infectious
  
  ###NEW cough to hand fomite calculations
  #ifelse(prob.glove>0, C.hand.S <- (((freq.hs * Surf_infect_cartons * f12) / (v.decay.min.hand*time.min) * 
   #                                     (time.min + (exp(-v.decay.min.hand * time.min)-1)/v.decay.min.hand))*((1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove)), 
   #      C.hand.S <- (((freq.hs * Surf_infect_cartons * f12) / (v.decay.min.hand*time.min) * (time.min + (exp(-v.decay.min.hand * time.min)-1)/v.decay.min.hand))*((1-hw.eff * prob.hw)^ Freq.hw)))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hhS <- freq.hf * hand.sa* (fingers.sa/face.sa) * C.hand.S * f23 * T* infectious
  
   ### hour 2
  C.hand.2c <- ((freq.hs * fomite2c * f12) / (v.decay.min.hand)) * 
    (1 - exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
 
   DT.hh2c <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.2c * f23 * T * infectious
 
  C.hand.2p <- ((freq.hs.plastic * fomite2p * f12) / (v.decay.min.hand)) * 
    (1-exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh2p <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.2p * f23 * T * infectious
  
  ### hour 3
  C.hand.3c <- ((freq.hs * fomite3c * f12) / (v.decay.min.hand)) * 
    (1-exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh3c <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.3c * f23 * T * infectious
  
  C.hand.3p <- ((freq.hs.plastic * fomite3p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh3p <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.3p * f23 * T * infectious
 
  ### hour 4
  C.hand.4c <- ((freq.hs * fomite4c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh4c <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.4c * f23 * T * infectious
  
  C.hand.4p <- ((freq.hs.plastic * fomite4p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh4p <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.4p * f23 * T * infectious
  
  ### hour 5
  C.hand.5c <- ((freq.hs * fomite5c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh5c <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.5c * f23 * T * infectious
  
  C.hand.5p <- ((freq.hs.plastic * fomite5p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh5p <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.5p * f23 * T * infectious
 
  ### hour 6
  C.hand.6c <- ((freq.hs * fomite6c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh6c <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.6c * f23 * T * infectious
  
  C.hand.6p <- ((freq.hs.plastic * fomite6p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh6p <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.6p * f23 * T * infectious
  
  ### hour 7
  C.hand.7c <- ((freq.hs * fomite7c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh7c <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.7c * f23 * T * infectious
  
  C.hand.7p <- ((freq.hs.plastic * fomite7p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh7p <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.7p * f23 * T * infectious
  
  ### hour 8
  C.hand.8c <- ((freq.hs * fomite8c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh8c <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.8c * f23 * T * infectious
  
  C.hand.8p <- ((freq.hs.plastic * fomite8p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh8p <- freq.hf * hand.sa * (fingers.sa/face.sa) * C.hand.8p * f23 * T * infectious
  
  ###next timestep with carryover - hour 4
  #ifelse(prob.glove>0, C.hand.4 <- (C.hand.3-DT.hh3 + ((((freq.hs * fomite4 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove))), 
    #     C.hand.4 <- (C.hand.3 - DT.hh3 + (((freq.hs * fomite4 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^ Freq.hw)))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hh4 <- freq.hf * (fingers.sa/hand.sa) * C.hand.4 * f23 * T
  
  ###next timestep with carryover - hour 5 
  #ifelse(prob.glove>0, C.hand.5 <- (C.hand.4-DT.hh4 + ((((freq.hs * fomite5 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove))), 
    #     C.hand.5 <- (C.hand.4 - DT.hh4 + (((freq.hs * fomite5 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^ Freq.hw)))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hh5 <- freq.hf * (fingers.sa/hand.sa) * C.hand.5 * f23 * T
  
  ###next timestep with carryover - hour 6
  #ifelse(prob.glove>0, C.hand.6 <- (C.hand.5-DT.hh5 + ((((freq.hs * fomite6 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove))), 
    #     C.hand.6 <- (C.hand.5 - DT.hh5 + (((freq.hs * fomite6 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^ Freq.hw)))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hh6 <- freq.hf * (fingers.sa/hand.sa) * C.hand.6 * f23 * T
  
  ###next timestep with carryover - hour 7
  #ifelse(prob.glove>0, C.hand.7 <- (C.hand.6-DT.hh6 + ((((freq.hs * fomite7 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove))), 
    #     C.hand.7 <- (C.hand.6 - DT.hh6 + (((freq.hs * fomite7 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^ Freq.hw)))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hh7 <- freq.hf * (fingers.sa/hand.sa) * C.hand.7 * f23 * T
  
  ###next timestep with carryover - hour 8
  #ifelse(prob.glove>0, C.hand.8 <- (C.hand.7-DT.hh7 + ((((freq.hs * fomite8 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove))), 
    #     C.hand.8 <- (C.hand.7 - DT.hh7 + (((freq.hs * fomite8 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*((1-hw.eff * prob.hw)^ Freq.hw)))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hh8 <- freq.hf * (fingers.sa/hand.sa) * C.hand.8 * f23 * T
  
  #### DOSE and RISK CALCULATIONS ####
  
  #Beta is the deposition fraction of infectious virus into the URT and lungs
  lungdep=mcstoc(runif, type = "V", min = 1, max = 1)
  #p is the inhalation rate (m3/h); from exposure handbook for moderate activity 2.7E-02 m3/min converted to per hr
  inhalerate=mcstoc(runif, type="V", min=1.62, max=3.18)
  
  ####### Exposure time (hr)
  exposuretime=1
  
  #Susceptible individual mas usage
  pwsms <- as.numeric(ifelse(Susceptible.Mask == "surgical", pwsms <- 1, pwsms <- 0))
  pwcms <- as.numeric(ifelse(Susceptible.Mask == "cloth", pwcms <- 1, pwcms <- 0))
  pwoN95s <- as.numeric(ifelse(Susceptible.Mask == "oN95", pwoN95s <- 1, pwoN95s <- 0))
  #pwdoubles <- as.numeric(ifelse(Susceptible.Mask == "double", pwdoubles <- 1, pwdoubles <- 0))
  
  #dosetime1 <- (inhaleCt1*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p + pwoN95s*oN95.mask.p ))
  #dosetime2 <- (inhaleCt2*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p + pwoN95s*oN95.mask.p ))
  #dosetime3 <- (inhaleCt3*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p + pwoN95s*oN95.mask.p ))
  #dosetime4 <- (inhaleCt4*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p + pwoN95s*oN95.mask.p ))
  #dosetime5 <- (inhaleCt5*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p + pwoN95s*oN95.mask.p ))
  #dosetime6 <- (inhaleCt6*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p + pwoN95s*oN95.mask.p ))
  #dosetime7 <- (inhaleCt7*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p + pwoN95s*oN95.mask.p ))
  #dosetime8 <- (inhaleCt8*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p + pwoN95s*oN95.mask.p ))
  
  aero.dose.df <- as.data.frame(cbind(DT.hh1c,DT.hh1p, DT.hh2c,DT.hh2p, DT.hh3c,DT.hh3p,
                        DT.hh4c,DT.hh4p, DT.hh5c,DT.hh5p, DT.hh6c,DT.hh6p, DT.hh7c,DT.hh7p, DT.hh8c,DT.hh8p, fomite1c,fomite1p, 
                        fomite2c, fomite2p,fomite3c, fomite3p,fomite4c, fomite4p,fomite5c, fomite5p,
                        fomite6c, fomite6p,fomite7c, fomite7p,fomite8c, fomite8p))
  
  
  aero.dose.df <- mutate( aero.dose.df, 
                       f1h = DT.hh1c + DT.hh1p,
                       f2h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p,
                       f3h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p,
                       f4h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p,
                       f5h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p + DT.hh5c + DT.hh5p,
                       f6h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p + DT.hh5c + DT.hh5p + DT.hh6c + DT.hh6p,
                       f7h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p + DT.hh5c + DT.hh5p + DT.hh6c + DT.hh6p + DT.hh7c + DT.hh7p,
                       f8h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p + DT.hh5c + DT.hh5p + DT.hh6c + DT.hh6p + DT.hh7c + DT.hh7p + DT.hh8c + DT.hh8p,
                       f1hconcentration = fomite1c + fomite1p,
                       f2hconcentration = fomite2c + fomite2p,
                       f3hconcentration = fomite3c + fomite3p,
                       f4hconcentration = fomite4c + fomite4p,
                       f5hconcentration = fomite5c + fomite5p,
                       f6hconcentration = fomite6c + fomite6p,
                       f7hconcentration = fomite7c + fomite7p,
                       f8hconcentration = fomite8c + fomite8p)
  
  return(aero.dose.df)
  
}



# Close Contact Function -----------------------------------------------------------

#Function for close contact module
Dosefunc <- function(Event, Volume.Fraction, Distance, Vol.Frac.Dist.Name, room.exchange, Clean1, Clean2, Clean3, Clean4, Clean5, Clean6, Clean7, Clean8, Humidity, Infected.Mask, HW, Glove, Susceptible.Mask){
  
  #Hydraulic Diameter - in (meters) using measurements from the Market paper (we may want to put these to a distribution or change to be representative of our frozen industry)
  facility.length <- 10
  facility.width <- 10
  facility.height <- 10
  
  #Inner volume of facility units is m3 - looks good
  facility.volume <- (facility.height * facility.length * facility.width)
  #Inner facility area m2
  facility.area = facility.length*facility.width
  
  ###viral concentration in saliva (log PFU/ml)
  ### viral concentration following calibration for super spreader event is mcstoc(runif, type="V",min=8.3, max=8.8)
  
  log.C.virus <- mcstoc(rtriang, type="V", min=6.1, mode=6.8, max=7.4)
  
  ###variant 2 times more viral particles, viral concentration in saliva, pfu/ml
  #variant2.log.C.virus <- mcstoc(rtriang, type="V", min=6.38, mode=7.1, max=7.67)
  ###variant 10 times more viral particles, viral concentration in saliva, pfu/ml
  #variant10.log.C.virus <- mcstoc(rtriang, type="V", min=7.1, mode=7.8, max=8.4)
  ###variant 100 times more viral particles, viral concentration in saliva, pfu/ml
  #variant100.log.C.virus <- mcstoc(rtriang, type="V", min=8.1, mode=8.8, max=9.4)
  ###variant 1000 times more viral particles, viral concentration in saliva, pfu/ml
  #variant1000.log.C.virus <- mcstoc(rtriang, type="V", min=9.1, mode=9.8, max=10.4)
  
  
  ###viral concentration in saliva, pfu/ml reduced for vaccination among infected worker 2.8 fold reduction
  #VR.virus <- mcstoc(rtriang, type="V", min=5.7, mode=6.4, max=7.0)
  
  #viral concentration in saliva, pfu/ml reduced for vaccination among infected worker 4.5 fold reduction
  #VM.virus <- mcstoc(rtriang, type="V", min=5.4, mode=6.1, max=6.7)
  
  C.virus <- log.C.virus
  #C.virus <- variant10.log.C.virus
  #C.virus <- variant100.log.C.virus
  #C.virus <- variant1000.log.C.virus
  #C.virus <- variant2.log.C.virus
  #C.virus <- VR.virus
  #C.virus <- VM.virus
  
   ##number of infected individuals present in facility
  people=1
  
  ### Droplet volume fractions
  ### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7126899/
  
  ##total fraction of volume of droplet sizes 50-60 ml/cough--for use in the droplet contact transmission COUGH event
  v50.60 <- mcstoc(rtriang, type = "V", min = 0.00000351079, mode = 0.00000596834, max = 0.0000066705)
  
  ##total fraction of volume of droplet sizes 60-100 ml/cough--for use in the droplet contact transmission COUGH event
  v60.100 <- mcstoc(rtriang, type = "V", min = 0.00000106218, mode = 0.00000491064, max = 0.00000841824)
  
  ##total fraction of volume of droplet sizes >100-750 ml/cough--for use in the droplet contact transmission COUGH event
  v100.750 <- mcstoc(rtriang, type = "V", min = 0.003976453, mode = 0.006784165, max = 0.007646274)
  
  ##total fraction of volume composed of droplet sizes <0.6 to 2.2 (ml/breathing) Papineni paper -- BREATHING event
  vbreath <- mcstoc(runif, type = "V", min =0.000000000105743, max =0.000000000289581)
  
  #Swtich to automatically select Volume Fraction
  VF <- mcstoc(rtriang, type = "V", min = 0, max = 0, mode = 0)
  ifelse(Volume.Fraction == "50-60", VF <- v50.60, 
         ifelse(Volume.Fraction == "60-100", VF <- v60.100,
                ifelse(Volume.Fraction == "100+", VF <- v100.750, 
                       ifelse(Volume.Fraction == "breath", VF <- vbreath, VF == 0))))
  
  #cough frequency: #coughs per hour per person
  coughfreq <- mcstoc(rempiricalD,type="V",values=c(0,1,2,3,4,5,6,7,8,9,10),prob=c(1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11 ))
  
  #breathing rate: #breaths per minute per person range 16-20 breaths/min, Fabian et al., Origin of Exhaled Breath Particles from Healthy and Human Rhinovirus-Infected 
  # converted to breaths per hour per person; this range represents both high and low producers
  breathrate <- mcstoc(runif, type="V", min=960, max=1200)
  
  ###Switch for event frequency
  ifelse(Event =="cough", eventfreq <- coughfreq, eventfreq <- breathrate)
  
  #### Mask intervention switches####
  
  pwsm <- as.numeric(ifelse(Infected.Mask == "surgical", pwsm <- 1, pwsm <- 0))
  pwcm <- as.numeric(ifelse(Infected.Mask == "cloth", pwcm <- 1, pwcm <- 0))
  pwoN95 <- as.numeric(ifelse(Infected.Mask == "oN95", pwoN95 <- 1, pwoN95 <- 0))
  #pwdouble <- as.numeric(ifelse(Infected.Mask == "double", pwdouble <- 1, pwdouble <- 0))
  ### Mask efficacies reported as log reductions
  
  ## (Source protection) surgical mask efficacy log reduction
  s.mask <- mcstoc(runif,type="V",min=0.387216,max=0.568636)
  
  ## (Source protection) cloth mask efficacy log reduction 
  c.mask <- mcstoc(runif,type="V",min=0.309804,max=0.619789)
  
  ## (Source protection) fit N95 efficacy log reduction
  oN95.mask <- mcstoc(runif, type="V", min=1.886057, max=2.154902)
  
  ## (Source protection) double masking (surgical mask then cloth) log reduction
  #double.mask <- mcstoc(runif, type="V", min=0.274, max=1.89)
  
  ### PPE mask efficacies reported as percent reduction for mask module
  ### Note: these are different parameters from the source protection
  
  s.mask.p <- mcstoc(runif, type="V", min=0.37, max=0.998)
  c.mask.p <- mcstoc(runif, type="V", min=0.17, max=0.887)
  oN95.mask.p <- 0.79
  #double.mask.p <- mcstoc(runif, type="V", min=0.40, max=0.9685)
  
  #### Surface Area Calculations ####
  
  #number of individual cartons contaminated (note each of these is only one side of the carton)
  contamcartons <- 0
  contamcartons<- as.numeric(ifelse(coughfreq == 0, contamcartons <- 0, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 1, contamcartons <- 1, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 2, contamcartons <- 2, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 3, contamcartons <- 3, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 4, contamcartons <- 4, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 5, contamcartons <- 5, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 6, contamcartons <- 6, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 7, contamcartons <- 7, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 8, contamcartons <- 8, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 9, contamcartons <- 9, contamcartons <- contamcartons))
  contamcartons<- as.numeric(ifelse(coughfreq == 10, contamcartons <- 10, contamcartons <- contamcartons))
  
  #SA of top of carton to be exposed by worker 1 and 2;range from small to larger cartons (m2)
  cartontop.SA<-mcstoc(runif, type = "V", min=0.106, max =0.116)
  
  #composite surface area of tops of cartons (m2) of contaminated individual cartons (event 1 and 2)
  contamcompcartSA<-cartontop.SA*contamcartons
  
  #SA for single carton (m2)
  cartonSA<-mcstoc(runif, type = "V", min=0.41, max =0.54)
  
  #total number of individual cartons processed in one hour (36 to 54 per pallet with 4 pallets per hour)
  carton.number<-mcstoc(runif, type = "V", min=144, max =216)
  
  #composite surface area of all cartons (event 1 and 2) processed in one hour
  compcartonSA<-cartonSA*carton.number
  
   #number of pallets contaminated by the second worker (note each of these is only one side of the pallet)
  contampallet <- 0
  contampallet<- as.numeric(ifelse(coughfreq == 0, contampallet <- 0, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 1, contampallet <- 1, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 2, contampallet <- 2, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 3, contampallet <- 3, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 4, contampallet <- 4, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 5, contampallet <- 5, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 6, contampallet <- 6, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 7, contampallet <- 7, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 8, contampallet <- 8, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 9, contampallet <- 9, contampallet <- contampallet))
  contampallet<- as.numeric(ifelse(coughfreq == 10, contampallet <- 10, contampallet <- contampallet))
  
  #number of pallets processed per hour
  pallets.processed<- 4
  
  #SA of the one side of plastic wrap (m2) this is representative of 10x10 or 9X4 or 9x6 cartons stacked
  singlesideplasticwrap.SA<-mcstoc(runif, type = "V", min=4.2, max =6.97)
  
  #SA of contaminated composite plastic wrap sides (factoring in the number of pallet sides contaminated by the second worker)
  contamcompplastic.SA<-singlesideplasticwrap.SA*contampallet
  
  #SA of the total plastic wrap (m2) this is representative of 10x10 or 9X4 or 9x6 cartons stacked
  plasticwrap.SA<-mcstoc(runif, type = "V", min=25.2, max =41.82)
  
  #SA of the total plastic wrap (m2) across the 4 pallets processed in 1h
  compplasticwrap.SA<-plasticwrap.SA*pallets.processed
  
  
  ### Total viral shedding (PFU/hr) with mask switches:
  ## Unit check: volume.fraction is ml/cough; C.virus is PFU/ml; cough.freq is coughs per hour per person
  E.virus <- (VF* eventfreq * people)*(10^(C.virus-(pwsm*s.mask)-(pwcm*c.mask)-(pwoN95*oN95.mask)))
  
  ### Sigma (PFU / hr):
  sigma <- E.virus
  
  #### Viral decay, settling velocity, and air exchange calculations ####
  
  ### viral decay per hour, updated to van Doremalen 2020 paper for aerosols (converted from TCID50 to PFU then
  ### calculated the half life and then converted to viral decay rate from the half life)
  ### also consistent with temp of 70 and RH of 40%, dhs.gov/science-and-technology/sars-airborne-calculator
  v.decay.lh <- 0.614
  #v.decay.hh <- 0.614
  #v.decay.hh <- 0.1447
  v.decay.hh <- 0
  
  ifelse(Humidity =="high", v.decay <- v.decay.hh, v.decay <- v.decay.lh)
  
  ###Ventilation in the facility
  
  #Converted from room exchanges per hour to get full room volume per hour (m3/hr) 
  room.exchange <- room.exchange
  Q <- (room.exchange*facility.volume)
  
  #Ac is the cross sectional area of the building (m2)
  Ac=facility.length * facility.width
  
  #Viral air concentration calculations
  #v_loss is in (m^3/hr)
  V.loss <- Q  + (v.decay * facility.volume)
  
  ### c(0) = 0, so we can assume that the concentration at time step 0 is 0 PFU / m^3
  Ct0 <- 0
  
  ### Assuming that each time step is one hour, therefore delta time is 1 hour
  delta.t <- 1
  
  #adding particle distance distribution
  ppdis<-mcstoc(runif,type="V", min=0.01, max=0.22)
 
  #particle probability distance
  pp <- 0.00
  
  ifelse(Volume.Fraction == "breath" & Distance == 0, pp <- 1.0, pp <- pp)
  ifelse(Volume.Fraction == "breath" & Distance == 1.0, pp <- 1.0, pp <- pp)
  ifelse(Volume.Fraction == "breath" & Distance == 2.0, pp <- 1.0, pp <- pp)
  ifelse(Volume.Fraction == "breath" & Distance == 3.0, pp <- 1.0, pp <- pp)
  
  ifelse(Volume.Fraction == "50-60" & Distance == 0, pp <- 1.0, pp <- pp)
  ifelse(Volume.Fraction == "50-60" & Distance == 1.0, pp <- 0.82, pp <- pp)
  ifelse(Volume.Fraction == "50-60" & Distance == 2.0, pp <- 0.43, pp <- pp)
  ifelse(Volume.Fraction == "50-60" & Distance == 3.0, pp <- 0.19, pp <- pp)
  
  ifelse(Volume.Fraction == "60-100" & Distance == 0, pp <- 1.0, pp <- pp)
  ifelse(Volume.Fraction == "60-100" & Distance == 1.0, pp <- 0.44, pp <- pp)
  ifelse(Volume.Fraction == "60-100" & Distance == 2.0, pp <- 0.01, pp <- pp)
  
  ifelse(Volume.Fraction == "100+" & Distance == 0, pp <- 1.0, pp <- pp)
  ifelse(Volume.Fraction == "100+" & Distance == 1.0, pp <- ppdis, pp <- pp)
  
  #fallout calculation based on pp for the air compartment
  #adding distribution for fallout
  ppfalldis<-mcstoc(runif,type="V",min= 0.07, max=0.78)
  
  #particle probability distance
  ppfall <- 0.00
  
  ifelse(Volume.Fraction == "breath" & Distance == 0, ppfall <- 0, ppfall <- ppfall)
  ifelse(Volume.Fraction == "breath" & Distance == 1.0, ppfall <- 0, ppfall <- ppfall)
  ifelse(Volume.Fraction == "breath" & Distance == 2.0, ppfall <- 0, ppfall <- ppfall)
  ifelse(Volume.Fraction == "breath" & Distance == 3.0, ppfall <- 0, ppfall <- ppfall)
  
  ifelse(Volume.Fraction == "50-60" & Distance == 0, ppfall <- 0, ppfall <- ppfall)
  ifelse(Volume.Fraction == "50-60" & Distance == 1.0, ppfall <- 0.18, ppfall <- ppfall)
  ifelse(Volume.Fraction == "50-60" & Distance == 2.0, ppfall <- 0.18, ppfall <- ppfall)
  ifelse(Volume.Fraction == "50-60" & Distance == 3.0, ppfall <- 0.09, ppfall <- ppfall)
  
  ifelse(Volume.Fraction == "60-100" & Distance == 0, ppfall <- 0, ppfall <- ppfall)
  ifelse(Volume.Fraction == "60-100" & Distance == 1.0, ppfall <- 0.55, ppfall <- ppfall)
  ifelse(Volume.Fraction == "60-100" & Distance == 2.0, ppfall <- 0.04, ppfall <- ppfall)
  
  ifelse(Volume.Fraction == "100+" & Distance == 0, ppfall <- 0, ppfall <- ppfall)
  ifelse(Volume.Fraction == "100+" & Distance == 1.0, ppfall <- ppfalldis, ppfall <- ppfall)
  
  ### Air viral concentration in the Market at hour t = 1 is (PFU/m^3). No carryover for close contact.
  ### At hour t = X:
  Ct1 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  Ct1pp<-Ct1*pp
  fallout1=Ct1*ppfall
  
  Ct2 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  Ct2pp<-Ct2*pp
  fallout2=Ct2*ppfall
  
  Ct3 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  Ct3pp<-Ct3*pp
  fallout3=Ct3*ppfall
  
  Ct4 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma 
  Ct4pp<-Ct4*pp
  fallout4=Ct4*ppfall
  
  Ct5 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma 
  Ct5pp<-Ct5*pp
  fallout5=Ct5*ppfall
  
  Ct6 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma  
  Ct6pp<-Ct6*pp
  fallout6=Ct6*ppfall
  
  Ct7 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma 
  Ct7pp<-Ct7*pp
  fallout7=Ct7*ppfall
  
  Ct8 <- (1 / V.loss) * (1 - exp( (-V.loss * delta.t) / facility.volume)) * sigma 
  Ct8pp<-Ct8*pp
  fallout8=Ct8*ppfall
  
  ############# Fomite Calculations and Surface Cleaning Switches ####
  
  # Percent Reduction for cleaners being used - ranging from 99.9% - 99.99% reduction (3 log to 4 log)
  sc.eff.p <- mcstoc(runif, type = "V", min = 0.999, max = 0.999)
  
  ###surface areas of stainless steel tables m^2; 0.5 by 0.5 (representative of arm reach)
  sstable.sa <- 0.25
  
  #Area of 3 finger tips touching the surface. 2.1 cm2 for 3 fingers from one hand; * 2 hands = 4.2 cm2
  #converted to m2; Bouwknegt et al 2015
  fingers.sa <- 0.00042
  
  #Area of two hands (palm only) in m^2 that touches fomite surface, 245 cm2 for one palm
  hand.sa <- 0.0490
  
  #designating floor around fomite m^2; 
  floor<-1
  
  #fomite stainless steel viral decay (from NEJM van Doremalen et al 2020); note: converted from TCID50 to PFU to viral decay rate
  fomite.decay <- 0.099139
  
  # new fomite decay calculated based off of plastic pitol et al; converted to days
  fd <- mcstoc(runif, type = "V", min = 0.14, max = 0.22)
  # new decay time in days
  decaytime <- mcstoc(runif, type = "V", min = 90, max = 90)
  
  # Fomite calculations; amount of virus (PFU) on a table at each time step. Final unit is PFU
    # (PFU/m3) * (m3) *(1/m2)*(m2/m2) * (m2/m2)* (unitless) * (1 - %Reduction * Unit-less)
  #final unit is PFU/m2
  #fomite from fallout associated with composite individual boxes
  fomite1c <- (((fallout1 * facility.volume)/contamcompcartSA) * (hand.sa/compcartonSA)* (1 - sc.eff.p * Clean1))*(exp(-fd*decaytime))
  fomite1p <- (((Ct1pp * facility.volume)/contamcompplastic.SA) * (hand.sa/ compplasticwrap.SA) * (1 - sc.eff.p * Clean1))*(exp(-fd*decaytime))
  fomite2c <- (((fallout2 * facility.volume)/contamcompcartSA) * (hand.sa/compcartonSA)* (1 - sc.eff.p * Clean2))
  fomite2p <- (((Ct2pp * facility.volume)/contamcompplastic.SA) * (hand.sa/ compplasticwrap.SA) * (1 - sc.eff.p * Clean2))
  fomite3c <- (((fallout3 * facility.volume)/contamcompcartSA) * (hand.sa/compcartonSA)* (1 - sc.eff.p * Clean3))
  fomite3p <- (((Ct3pp * facility.volume)/contamcompplastic.SA) * (hand.sa/ compplasticwrap.SA) * (1 - sc.eff.p * Clean3))
  fomite4c <- (((fallout4 * facility.volume)/contamcompcartSA) * (hand.sa/compcartonSA)* (1 - sc.eff.p * Clean4))
  fomite4p <- (((Ct4pp * facility.volume)/contamcompplastic.SA) * (hand.sa/ compplasticwrap.SA) * (1 - sc.eff.p * Clean4))
  fomite5c <- (((fallout5 * facility.volume)/contamcompcartSA) * (hand.sa/compcartonSA)* (1 - sc.eff.p * Clean5))
  fomite5p <- (((Ct5pp * facility.volume)/contamcompplastic.SA) * (hand.sa/ compplasticwrap.SA) * (1 - sc.eff.p * Clean5))
  fomite6c <- (((fallout6 * facility.volume)/contamcompcartSA) * (hand.sa/compcartonSA)* (1 - sc.eff.p * Clean6))
  fomite6p <- (((Ct6pp * facility.volume)/contamcompplastic.SA) * (hand.sa/ compplasticwrap.SA) * (1 - sc.eff.p * Clean6))
  fomite7c <- (((fallout7 * facility.volume)/contamcompcartSA) * (hand.sa/compcartonSA)* (1 - sc.eff.p * Clean7))
  fomite7p <- (((Ct7pp * facility.volume)/contamcompplastic.SA) * (hand.sa/ compplasticwrap.SA) * (1 - sc.eff.p * Clean7))
  fomite8c <- (((fallout8 * facility.volume)/contamcompcartSA) * (hand.sa/compcartonSA)* (1 - sc.eff.p * Clean8))
  fomite8p <- (((Ct8pp * facility.volume)/contamcompplastic.SA) * (hand.sa/ compplasticwrap.SA) * (1 - sc.eff.p * Clean8))
  
  #fomite3 <- ((((fomite2 + ((fallout3) * (facility.volume) * (sstable.sa/floor) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean3)))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite4 <- ((((fomite3 + ((fallout4) * (facility.volume) * (sstable.sa/floor) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean4)))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite5 <- ((((fomite4 + ((fallout5) * (facility.volume) * (sstable.sa/floor) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean5)))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite6 <- ((((fomite5 + ((fallout6) * (facility.volume) * (sstable.sa/floor) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean6)))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite7 <- ((((fomite6 + ((fallout7) * (facility.volume) * (sstable.sa/floor) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean7)))/fomite.decay)*(1-exp(-fomite.decay*1))
  #fomite8 <- ((((fomite7 + ((fallout8) * (facility.volume) * (sstable.sa/floor) * (hand.sa/sstable.sa))) * (1 - sc.eff.p * Clean8)))/fomite.decay)*(1-exp(-fomite.decay*1))
  
  #Hand Hygiene Module linking in with fomite/surface cleaning code
  
  ## hand wash efficacy log removal (log10 PFU removal 2 log) --converted to %reduction
  hw.eff <- mcstoc(rtriang, type="V", min=0.99, mode=0.99, max=.99)
  
  # Probability of harvester handwashing
  prob.hw <- ifelse(HW == "yes", prob.hw <- 1, prob.hw <- 0)
  
  # Frequency of washing hands per hour (converted to minutes to work with math)
  Freq.hw <- mcstoc(rempiricalD,type="V",values=c(1,1,1,1,1),prob=c(0.2,0.2,0.2,0.2,0.2))
  
  # Probability of harvester using gloves (Jaykus et al 2009)
  # Probability of using gloves
  prob.glove <- ifelse(Glove == "yes", prob.glove <- 1, prob.glove <- 0)
  
  # Frequency of glove changes per hour (converted to minutes to work with math)
  freq.glove <- mcstoc(rempiricalD,type="V",values=c(1,1,1,1,1),prob=c(0.2,0.2,0.2,0.2,0.2))
  
  #fraction of pathogen transfer from hand to glove (4) during gloving process
  #https://www-ncbi-nlm-nih-gov.proxy.library.emory.edu/pmc/articles/PMC4136105/
  f24<-mcstoc(runif,type="V",min=1,max=1)
  
  #Area of fomite in m^2 that touches hand surface
  #Assuming to be equal to the area of hand that touches the fomite
  
  #Fraction of pathogens transferred from fomite (1) to hand (2) 
  #NOTE: f12.hh is in a high humidity environment (40% - 65%)
  #NOTE: f12.lh is in a low humidity environment (15% - 32%)
  #updated f12.hh to be for acrylic, same paper, low humidity
  #Per https://dx.doi.org/10.1128%2FAEM.01030-13 Using MS2
  f12.hh <- mcstoc(rtriang, type = "V", min =0.067, mode = 0.217, max =0.367)
  #f12.hh <- mcstoc(rnorm, type = "V", mean = 0.374, sd = 0.16)
  f12.lh <- mcstoc(rtriang, type = "V", min = 0, mode=0.069, max= 0.158)
  
  #Humidity switch
  f12 <-ifelse(Humidity =="high", f12 <- f12.hh, f12<- f12.lh)
  
  #Fraction of pathogens transferred from hand (2) to fomite (1)
  #Using this parameter from Alicia's EID paper
  f21 <- 0.025
  
  #Fraction of pathogens transferred from hand (2) to face (3) 
  #Updated based on Julian et al MS2 paper (same as above for transfer efficiencies)
  #NOTE: This was not RH dependent
  f23 <- mcstoc(rtriang, type = "V", min =0.137, mode = 0.20, max =0.263)
  
  #Frequency of contacts from hand to individual cartons per hour converted to (per minute)
  freq.hs <- (carton.number/60)
  
  #Frequency of contacts from hand to plastic wrap; distribution in times per hour converted to (per minute); equivalent to 1 touch per pallet up to 5 per pallet
  freq.hs.plastic <- mcstoc(runif, type = "V", min=4/60, max =20/60)
  
  #Frequency of contacts from hand to face (per min) 
  freq.hf <- 0.8
  
  #updated viral inactivation rate on hands per Nicas et al (see below): based off of influenza A. Note: this is already in units min^-1
  #https://onlinelibrary-wiley-com.proxy.library.emory.edu/doi/full/10.1111/j.1539-6924.2009.01253.x
  v.decay.skin<- mcstoc(runif, type = "V", min=0.92, max =1.47)
  v.decay.min <- v.decay.skin
  
  #converted from cm2 to m2 (all from https://aem.asm.org/content/aem/84/18/e00709-18.full.pdf)
  eyes.sa<-mcstoc(runif, type = "V", min =1E-5 , max =2E-4 )
  nose.sa<-mcstoc(runif, type = "V", min = 1E-5, max =1E-3 )
  mouth.sa <-mcstoc(runif, type = "V", min = 1E-4, max =4.1E-3 )
  
  face.sa<-eyes.sa+nose.sa+mouth.sa
  
    #Viral decay in minutes based on hand surface (collapsed this into one line)
  v.decay.min.hand <- v.decay.min + freq.hs * f21 + freq.hf * f23
  
   #Concentration of pathogens on hand (total PFU)
  #Assume at time 0, C(t) = 0
  time.min = 60
  T = mcstoc(runif, type = "V", min = 60, max = 60)
  
  #Building in handwashing (freq efficacy prob of hw) and glove use/changes into the same step at hour 1. Units are in PFU
 #concentation on hand and dose from events 1 and 2 composite
  #units on hand are PFU/m2
   C.hand.1c <- ((freq.hs * fomite1c * f12) / (v.decay.min.hand)) * 
     (1-exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
         
  #fraction of dose infectious, plan to run 0.01, 0.1, and 1
   infectious<- 1
  
  #dose has units of PFU 
  DT.hh1c <- freq.hf * hand.sa * (fingers.sa/face.sa)* C.hand.1c * f23 * T * infectious
  
  C.hand.1p <- ((freq.hs.plastic * fomite1p * f12) / (v.decay.min.hand)) * 
                                                          (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
        
  DT.hh1p <- freq.hf * hand.sa* (fingers.sa/face.sa)*C.hand.1p * f23 * T * infectious
  
 #hour 2
  C.hand.2c <- ((freq.hs * fomite2c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh2c <- freq.hf * hand.sa * (fingers.sa/face.sa)* C.hand.2c * f23 * T * infectious
  
  C.hand.2p <- ((freq.hs.plastic * fomite2p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh2p <- freq.hf * hand.sa* (fingers.sa/face.sa)*C.hand.2p * f23 * T * infectious
  
  #hour 3
  C.hand.3c <- ((freq.hs * fomite3c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh3c <- freq.hf * hand.sa * (fingers.sa/face.sa)* C.hand.3c * f23 * T * infectious
  
  C.hand.3p <- ((freq.hs.plastic * fomite3p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh3p <- freq.hf * hand.sa* (fingers.sa/face.sa)*C.hand.3p * f23 * T * infectious
  
  #hour 4
  C.hand.4c <- ((freq.hs * fomite4c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh4c <- freq.hf * hand.sa * (fingers.sa/face.sa)* C.hand.4c * f23 * T * infectious
  
  C.hand.4p <- ((freq.hs.plastic * fomite4p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh4p <- freq.hf * hand.sa* (fingers.sa/face.sa)*C.hand.4p * f23 * T * infectious
  
  #hour 5
  C.hand.5c <- ((freq.hs * fomite5c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh5c <- freq.hf * hand.sa * (fingers.sa/face.sa)* C.hand.5c * f23 * T * infectious
  
  C.hand.5p <- ((freq.hs.plastic * fomite5p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh5p <- freq.hf * hand.sa* (fingers.sa/face.sa)*C.hand.5p * f23 * T * infectious
  
  #hour 6
  C.hand.6c <- ((freq.hs * fomite6c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh6c <- freq.hf * hand.sa * (fingers.sa/face.sa)* C.hand.6c * f23 * T * infectious
  
  C.hand.6p <- ((freq.hs.plastic * fomite6p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh6p <- freq.hf * hand.sa* (fingers.sa/face.sa)*C.hand.6p * f23 * T * infectious
  
  #hour 7
  C.hand.7c <- ((freq.hs * fomite7c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh7c <- freq.hf * hand.sa * (fingers.sa/face.sa)* C.hand.7c * f23 * T * infectious
  
  C.hand.7p <- ((freq.hs.plastic * fomite7p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh7p <- freq.hf * hand.sa* (fingers.sa/face.sa)*C.hand.7p * f23 * T * infectious
  
  #hour 8
  C.hand.8c <- ((freq.hs * fomite8c * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh8c <- freq.hf * hand.sa * (fingers.sa/face.sa)* C.hand.8c * f23 * T * infectious
  
  C.hand.8p <- ((freq.hs.plastic * fomite8p * f12) / (v.decay.min.hand)) * 
    (1- exp(-v.decay.min.hand * time.min))*((1-hw.eff * prob.hw)^Freq.hw) 
  
  DT.hh8p <- freq.hf * hand.sa* (fingers.sa/face.sa)*C.hand.8p * f23 * T * infectious
   
  ###next timestep with carryover - hour 3
 # ifelse(prob.glove>0, C.hand.3 <- (C.hand.2-DT.hh2 + ((((freq.hs * fomite3 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove)), 
  #       C.hand.3 <- (C.hand.2 - DT.hh2 + (((freq.hs * fomite3 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^ Freq.hw))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hh3 <- freq.hf * (fingers.sa/hand.sa)*  C.hand.3 * f23 * T
  
  ###next timestep with carryover - hour 4
  #ifelse(prob.glove>0, C.hand.4 <- (C.hand.3-DT.hh3 + ((((freq.hs * fomite4 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove)), 
    #     C.hand.4 <- (C.hand.3 - DT.hh3 + (((freq.hs * fomite4 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^ Freq.hw))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hh4 <- freq.hf * (fingers.sa/hand.sa)*C.hand.4 * f23 * T
  
  ###next timestep with carryover - hour 5
  #ifelse(prob.glove>0, C.hand.5 <- (C.hand.4-DT.hh4 + ((((freq.hs * fomite5 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove)), 
    #     C.hand.5 <- (C.hand.4 - DT.hh4 + (((freq.hs * fomite5 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^ Freq.hw))
  
  #Dose Transferred to facial membrane over interval [0,T]
#  DT.hh5 <- freq.hf * (fingers.sa/hand.sa)* C.hand.5 * f23 * T
  
  ###next timestep with carryover - hour 6
 # ifelse(prob.glove>0, C.hand.6 <- (C.hand.5-DT.hh5 + ((((freq.hs * fomite6 * f12) / v.decay.min.hand) * 
  #                                                        (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove)), 
   #      C.hand.6 <- (C.hand.5 - DT.hh5 + (((freq.hs * fomite6 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^ Freq.hw))
  
  #Dose Transferred to facial membrane over interval [0,T]
  #DT.hh6 <- freq.hf * (fingers.sa/hand.sa)* C.hand.6 * f23 * T
  
  ###next timestep with carryover - hour 7
  #ifelse(prob.glove>0, C.hand.7 <- (C.hand.6-DT.hh6 + ((((freq.hs * fomite7 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove)), 
    #     C.hand.7 <- (C.hand.6 - DT.hh6 + (((freq.hs * fomite7 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^ Freq.hw))
  
  #Dose Transferred to facial membrane over interval [0,T]
 # DT.hh7 <- freq.hf * (fingers.sa/hand.sa)* C.hand.7 * f23 * T
  
  ###next timestep with carryover - hour 8
  #ifelse(prob.glove>0, C.hand.8 <- (C.hand.7-DT.hh7 + ((((freq.hs * fomite8 * f12) / v.decay.min.hand) * 
   #                                                       (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^Freq.hw)*((1-f24)^freq.glove)), 
    #     C.hand.8 <- (C.hand.7 - DT.hh7 + (((freq.hs * fomite8 * f12) / v.decay.min.hand) * (1 - exp(-v.decay.min.hand * time.min)))*(1-hw.eff * prob.hw)^ Freq.hw))
  
  #Dose Transferred to facial membrane over interval [0,T]
#  DT.hh8 <- freq.hf * (fingers.sa/hand.sa)* C.hand.8 * f23 * T
  
  ##############Susceptible Mask and Dose Calculations######################################
  #Beta is the deposition fraction of infectious virus into the lungs
  lungdep=mcstoc(runif, type = "V", min = 1, max = 1)
  
  #p is the inhalation rate (m3/h); from exposure handbook for moderate activity 2.7E-02 m3/min converted to per hr
  inhalerate=mcstoc(runif, type="V", min=1.62, max=3.18)
  
  ####### unit in hours
  exposuretime=1
  
  pwsms <- as.numeric(ifelse(Susceptible.Mask == "surgical", pwsms <- 1, pwsms <- 0))
  pwcms <- as.numeric(ifelse(Susceptible.Mask == "cloth", pwcms <- 1, pwcms <- 0))
  pwoN95s <- as.numeric(ifelse(Susceptible.Mask == "oN95", pwoN95s <- 1, pwoN95s <- 0))
  #pwdoubles <- as.numeric(ifelse(Susceptible.Mask == "double", pwdoubles <- 1, pwdoubles <- 0))
  
  #Dose Calculations with Mask interventions
  #dose in PFU = (PFU/m3 * unitless * m3/hr * hr)*(unitless*%Reduction)
  
  #dosetime1 <- (Ct1pp*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p +  pwoN95s*oN95.mask.p ))
  #dosetime2 <- (Ct2pp*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p +  pwoN95s*oN95.mask.p ))
  #dosetime3 <- (Ct3pp*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p +  pwoN95s*oN95.mask.p ))
  #dosetime4 <- (Ct4pp*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p +  pwoN95s*oN95.mask.p ))
  #dosetime5 <- (Ct5pp*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p +  pwoN95s*oN95.mask.p ))
  #dosetime6 <- (Ct6pp*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p +  pwoN95s*oN95.mask.p ))
  #dosetime7 <- (Ct7pp*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p +  pwoN95s*oN95.mask.p ))
  #dosetime8 <- (Ct8pp*lungdep*inhalerate*exposuretime)*(1 - (pwsms*s.mask.p + pwcms*c.mask.p +  pwoN95s*oN95.mask.p ))
  
  dose.df <- as.data.frame(cbind(DT.hh1c,DT.hh1p, DT.hh2c,DT.hh2p, DT.hh3c,DT.hh3p,
                                 DT.hh4c,DT.hh4p, DT.hh5c,DT.hh5p, DT.hh6c,DT.hh6p, DT.hh7c,DT.hh7p, DT.hh8c,DT.hh8p, fomite1c,fomite1p, 
                                 fomite2c, fomite2p,fomite3c, fomite3p,fomite4c, fomite4p,fomite5c, fomite5p,
                                 fomite6c, fomite6p,fomite7c, fomite7p,fomite8c, fomite8p))
  
  
  dose.df <- mutate(   dose.df, 
                        fVolFracDist1h = DT.hh1c + DT.hh1p,
                        fVolFracDist2h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p,
                        fVolFracDist3h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p,
                        fVolFracDist4h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p,
                        fVolFracDist5h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p + DT.hh5c + DT.hh5p,
                        fVolFracDist6h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p + DT.hh5c + DT.hh5p + DT.hh6c + DT.hh6p,
                        fVolFracDist7h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p + DT.hh5c + DT.hh5p + DT.hh6c + DT.hh6p + DT.hh7c + DT.hh7p,
                        fVolFracDist8h = DT.hh1c + DT.hh1p + DT.hh2c + DT.hh2p + DT.hh3c + DT.hh3p + DT.hh4c + DT.hh4p + DT.hh5c + DT.hh5p + DT.hh6c + DT.hh6p + DT.hh7c + DT.hh7p + DT.hh8c + DT.hh8p,
                        fconcVolFracDist1h = fomite1c + fomite1p,
                        fconcVolFracDist2h = fomite2c + fomite2p,
                        fconcVolFracDist3h = fomite3c + fomite3p,
                        fconcVolFracDist4h = fomite4c + fomite4p,
                        fconcVolFracDist5h = fomite5c + fomite5p,
                        fconcVolFracDist6h = fomite6c + fomite6p,
                        fconcVolFracDist7h = fomite7c + fomite7p,
                        fconcVolFracDist8h = fomite8c + fomite8p)
  
  names(dose.df) <- gsub("VolFracDist", Vol.Frac.Dist.Name,names(dose.df))
  
  return(dose.df)
  
}
