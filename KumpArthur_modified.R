# ---------------------------------------------------------- #
# ************ KUMP & ARTHUR STYLE BOX MODEL *************** #
# This script builds on the equations of Kump and Arthur,    #
# 1999 (Chemical Geology) [herein referred to as KA99]       #
# to develop a box model of the C cycle that solves          #
# transient perturbations to the system.                     #
# ---------                                                  #
# T. Kukla (Stanford Univ. 2019)                             #
# ---------------------------------------------------------- #
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)

# *************************************************************************************************** #
## -- CONFIGURATION NOTES -- ##
## [1] Forward model solving carbon fluxes and their isotopes
## [2] Fborg is a function of Phosphorus availability and burial
## [3] DB (fractionation difference b/t org and inorg) does not respond
#      to changes in the phosphorus reservoir but does respond to 
#      changes in pCO2
## [4] U is a factor increase in mountain erosion (set to 1 as default)
## [5] We do NOT solve carbonate chemistry
## [6] Three possible perturbations: 
# -------- (1) a C-cycle perturbation of some factor with a set isotopic composition of -6 per mille
# -------- (2) a phosphorus recycling perturbation 
# -------- (3) a mountain erosion perturbation 
# *************************************************************************************************** #

rm(list=ls())     # clear global environment
setwd('/Users/tylerkukla/Documents/Stanford/Teaching_TA/EarthSysEvol_winter2019/Lecture_3_KumpArthur')


# -------- SET INITIAL CONDITIONS --------
#... *****Time domain (x1,000 yrs)*****
t.start <- 1         # year to start model run (in kyr)
t.end <- 2500        # year to end model run (in kyr)
dt <- 1              # time step (kyr)
t_kyr <- seq(t.start, t.end, by=dt)  # sequence of time steps
dur <- length(t_kyr) # number of time steps


#... *****Fluxes (x10^12mol C/kyr)*****
# CARBON FLUXES
Fw.carb <- vector(length = dur) ; Fw.carb[] <- 34000    # [molC kyr-1] carbonate weathering flux
Fw.org <- vector(length = dur) ; Fw.org[] <- 10000      # [molC kyr-1] organic C weathering flux
Fw <- vector(length = dur) ; Fw[] <- Fw.carb[1] + Fw.org[1]    # [molC kyr-1] total weathering flux
Fw.sil <- vector(length = dur) ; Fw.sil[] <- 6000       # [molC kyr-1] silicate weathering flux
Fvolc <- vector(length = dur) ; Fvolc[] <- 6000         # [molC kyr-1] volcanic degassing (+ metamorphism) C flux
Fb.org <- vector(length = dur) ; Fb.org[] <- 10000      # [molC kyr-1] burial of org C
Fb.carb <- vector(length = dur) ; Fb.carb[] <- 40000    # [molC kyr-1] burial of carbonate C
F.in <- vector(length = dur) ; F.in[] <- Fw[1] + Fvolc[1]  # [molC kyr-1] total input of all weathered and volc C (approximately mantle isotopic comp) (Fw.prime in KA99)

# PHOSPHORUS FLUXES
Fw.phos <- vector(length = dur) ; Fw.phos[] <- 47       # [molP04 kyr-1] weathered PO4
Fb.phos <- vector(length = dur) ; Fb.phos[] <- 47       # [molP04 kyr-1] buried PO4
sil.P <- 0.58     # Fraction of P delivered by silicate weathering (Shields and Mills 2017--PNAS)
carb.P <- 0.21    # Fraction of P delivered by carbonate weathering (Shields and Mills 2017--PNAS)
org.P <- 0.21     # Fraction of P delivered by organic weathering (Shields and Mills 2017--PNAS)


#... *****Reservoirs*****
# CARBON AND PHOSPHORUS
pCO2 <- rep(560, dur)          # [ppmv] atmospheric CO2 
RCO2 <- pCO2 / pCO2            # factor change in pCO2 (for SilWx feedback)
pCO2.res <- vector(length = dur) ; pCO2.res[] <- 100000          # [molC 10^12] atmospheric pCO2 reservoir
HCO3.res <- vector(length = dur) ; HCO3.res[] <- 3.8e6           # [molC 10^12] oceanic bicarbonate reservoir 
M0 <- vector(length = dur) ; M0[1] <- pCO2.res[1] + HCO3.res[1]  # [molC 10^12] mass of inorganic carbonate box
PO4 <- rep(3100, dur)          # [molP 10^12] mass of phosphate


#... *****Scaling exponents and factors*****
# ALPHAS - weathering / erosion exponents
sil.exp <- 0.55   # Uplift exponent for silicate weathering
carb.exp <- 0.9   # Uplift exponent for carbonate weathering
worg.exp <- 0.7   # Uplift exponent for organic weathering
borg.exp <- 0.15  # Uplift exponent for organic carbon burial

# SCALING PARAMETERS
carb.clim.exp <- 0.05      # Exponent carbonate weathering Arrhenius reaction (not implemented as a changeable parameter)
R.n <- rep(0.8, dur)       # [0-to-1] Exponent for silicate reactivity of land surface (Caves et al., 2016--EPSL)
U <- rep(1, dur)           # Mountain erosion factor increase (uplift param)
P.exp <- rep(1, dur)       # Exponent for phosphorus recycling 


#... *****Isotopes*****
DB <- vector(length = dur) ; DB[] <- (((159.5*0.25) + 38.39) / (0.034*pCO2[1])) - 33     # [per mille] inorg C to org C fractionation (Delta-b in KA99)
d.carb <- vector(length = dur) ; d.carb[] <- 0            # [per mille] C isotopic composition of carbonate
d.in <- vector(length = dur) ; d.in[] <- -6               # [per mille] C isotopic composition of weathered / volc inputs (approximately mantle value) (delta prime in paper)
d.org <- vector(length = dur) ; d.org[] <- DB[1] + d.carb[1]     # [per mille] C isotopic composition of organic matter


# -------- DEFINE A PERTURBATION TO THE SYSTEM --------
## pert_fun -- This function helps calculate a perturbation vector and can be used on any input parameter
#  INPUTS: (1) InitialValue: the value of the perturbed parameter in its un-perturbed (initial) state
#              [!!!-- for a flux perturbation InitialValue = 0 --!!!]
#          (2) Pstart: the location of the year in "t_kyr" when the perturbation starts (if dt=1 then this is just the year the perturbation starts)
#          (3) Pend: the same but for the end of the perturbation (Pend must be greater than Pstart)
#          (4) FactorChange: the factor that we multiply the InitialValue by during the years of the perturbation
#          (5) t_kyr: NO NEED TO CHANGE - set as default to the time-series previously defined
pert_fun <- function(InitialValue=0, Pstart, Pend, NewValue, timeLoc=t_kyr){
  if(Pstart > Pend){stop("Pstart must be less than Pend!")}
  perturbVec <- vector(length=length(t_kyr))
  perturbVec <- ifelse(timeLoc >= Pstart & timeLoc <= Pend, NewValue, InitialValue)
}

# EXAMPLE PERTURBATIONS
## [VOLCANISM]
F.extra <- pert_fun(InitialValue <- 0, Pstart <- 1e3, Pend <- 1.5e3, NewValue <- Fvolc[1]*2)   # [molC kyr] change in a flux of carbon (often the volcanic flux)
d13C.extra <- -6      # [per mille] isotopic composition of the additional carbon
P.exp <- pert_fun(InitialValue <- P.exp, Pstart <- 1e3, Pend <- 1.3e3, NewValue <- P.exp[1]*1.25)  # [] increase in phosphorus recycling

## [METHANE]
F.extra <- pert_fun(InitialValue <- 0, Pstart <- 1e3, Pend <- 1.5e3, NewValue <- Fvolc[1]/2)   # [molC kyr] change in a flux of carbon (often the volcanic flux)
d13C.extra <- -55      # [per mille] isotopic composition of the additional carbon
P.exp <- pert_fun(InitialValue <- P.exp, Pstart <- 1e3, Pend <- 1.3e3, NewValue <- P.exp[1]*1)  # [] increase in phosphorus recycling

## [MTN UPLIFT]
# F.extra <- pert_fun(InitialValue <- 0, Pstart <- 1e3, Pend <- 1.5e3, NewValue <- 0)   # [molC kyr] change in a flux of carbon (often the volcanic flux)
# d13C.extra <- -6      # [per mille] isotopic composition of the additional carbon
# U <- pert_fun(InitialValue <- U, Pstart <- 1e3, Pend <- 3e3, NewValue <- 2.5)  # [] increase in phosphorus recycling



# **************************************************************************** #
# ----------------------------- RUN THE MODEL -------------------------------- #
# **************************************************************************** #
# ------- C-CYCLE transient box model -------
#... loop from time 2 (time 1 is already defined) to the last timestep (which is equal to "dur")
for(tstep in 2:dur){
  #... FIRST calculate the fluxes at this time step
  # Volcanism
  Fvolc[tstep] <- Fvolc[tstep]
  
  # Silicate Weathering
  Fw.sil[tstep] <- Fw.sil[1] * (U[tstep]**sil.exp) * (pCO2[tstep-1] / pCO2[1])**R.n[tstep]   # simplified functional form of the feedback
  
  # Carbonate Fluxes--equations are from Shields and Mills 2017 (PNAS) [assuming constant temp; set first 288 to "temp" (then define temp from pCO2) if you want to change temp]
  Fw.carb[tstep] <- Fw.carb[1]* (U[tstep]**carb.exp) *exp(carb.clim.exp*(288-288))
  Fb.carb[tstep] <- Fw.carb[tstep] + Fw.sil[tstep]    # quasi steady state assumption of KA99
  
  # Organic and P weathering fluxes
  Fw.org[tstep] <- Fw.org[1] * (U[tstep]**worg.exp)
  Fw.phos[tstep] <- Fw.phos[1] * (((sil.P)*(Fw.sil[tstep]/Fw.sil[1])) + ((carb.P)*(Fw.carb[tstep]/Fw.carb[1])) + ((org.P)*(Fw.org[tstep]/Fw.org[1])))
  
  # Organic and P burial fluxes
  Fb.org[tstep] <- Fb.org[1] * (U[tstep]**borg.exp) * ((PO4[tstep-1]/PO4[1])**P.exp[tstep])
  Fb.phos[tstep] <- Fb.phos[1] * (Fb.org[tstep]/Fb.org[1])**(1/P.exp[tstep])
  
  #... SECOND calculate the reservoir sizes using finite difference
  M0[tstep] <- M0[tstep-1] + (Fw.org[tstep] + Fvolc[tstep] + F.extra[tstep] - Fb.org[tstep] - Fw.sil[tstep])*dt
  PO4[tstep] <- PO4[tstep-1] + (Fw.phos[tstep] - Fb.phos[tstep])*dt
  
  #... THIRD scale the new mass of the inorganic C reservoir to pCO2
  pCO2[tstep] <- ((M0[tstep]/M0[1])**2) * pCO2[1]
  
  #... FOURTH calculate isotopic values 
  Fw.prime <- Fw.carb[tstep] + Fw.org[tstep] + Fvolc[tstep]             # the Fw_prime flux from KA99
  DB[tstep] <- (((159.5*0.25) + 38.39) / (0.034*pCO2[tstep])) - 33  # the org-inorg fractionation (Delta_B in KA99)
  #... and solve d13C carb
  d.carb[tstep] <- d.carb[tstep-1] + ((Fw.prime*(d.in[tstep-1] - d.carb[tstep-1]) - Fb.org[tstep]*DB[tstep] + F.extra[tstep]*(d13C.extra - d.carb[tstep-1])) / M0[tstep])*dt
  #... and now d13C org
  d.org[tstep] <- DB[tstep] + d.carb[tstep]
  
  # -------- and loop through again! 

}

#... when we're all done, we will bring the output together in a nice, single data frame
df <- as_tibble(cbind(t_kyr, M0, pCO2, d.carb, d.org, DB, PO4, Fb.phos, Fb.org, Fw.phos, 
                      Fw.org, Fb.carb, Fw.carb, Fw.sil, Fvolc, F.extra, P.exp, U, R.n))



# --------------------- DATA VISUALIZATION ------------------------
# PERTURBATIONS
#... F-extra
p_Fext <- ggplot(df) + 
  geom_line(aes(x=t_kyr, y=F.extra), color='#2A2C2B', size=1) +
  labs(x="Time (x1000 yrs)", y="Carbon Perturbation [mol C kyr-1]") +
  theme_linedraw()

#... phosphorus recycling
p_Pexp <- ggplot(df) + 
  geom_line(aes(x=t_kyr, y=P.exp), color='#2A2C2B', size=1) +
  labs(x="Time (x1000 yrs)", y="Phosphorus recycling exponent") +
  theme_linedraw()

#... pCO2 
p_pCO2 <- ggplot(df) + 
  geom_line(aes(x=t_kyr, y=pCO2), color='#8E2800', size=1) +
  labs(x="Time (x1000 yrs)", y=expression("pCO"['2'])) +
  theme_linedraw()

#... OUTPUT
#... d.carb
p_dcarb <- ggplot(df) + 
  geom_line(aes(x=t_kyr, y=d.carb), color='#3E606F', size=1) +
  labs(x="Time (x1000 yrs)", y=expression('δ'^'13'*'C'['Carbonate'])) +
  theme_linedraw()

#... d.org
p_dorg <- ggplot(df) + 
  geom_line(aes(x=t_kyr, y=d.org), color='#468966', size=1) +
  labs(x="Time (x1000 yrs)", y=expression('δ'^'13'*'C'['Organic'])) +
  theme_linedraw()

#... Fw.sil
p_fwsil <- ggplot(df) + 
  geom_line(aes(x=t_kyr, y=Fw.sil), color='#EB7F00', size=1) +
  labs(x="Time (x1000 yrs)", y="Silicate Weathering Flux") +
  theme_linedraw()


# ----- arrange it all together
ggarrange(p_Fext, p_Pexp, p_pCO2, p_fwsil, p_dcarb, p_dorg, ncol=2, nrow=3, align="hv")

