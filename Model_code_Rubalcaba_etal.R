######### A model to simulate steady-steady oxygen consumption of oxygen in aquatic ectotherms ###########
## Reference:
# Rubalcaba, J. G., Verberk, W. C. E. P., Hendriks, A. J., Saris, B., Woods, H. A. 
# "Oxygen limitation may affect the temperature- and size-dependence of metabolism in aquatic ectotherms". PNAS.
######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######

OxBalance <- function(Tw,    # Water temperature (ºC)
                      S=35,  # Salinity (parts per thousand)
                      v=0.01,# water velocity (m/s)
                      Co=1,  # Oxygen concentration (0-1, where 1=equilibrium water-atmosphere)
                      M,     # Body mass (g)
                      b0=1e5,# Normalization constant MTE
                      E=0.6, # Activation energy (eV)
                      a=1e-5,# Scaling constant (a) /exponent (b) (Surface area ~ a*Mass^b)
                      b=0.7,
                      delta=1,# Increase in ox demand MMR = delta*RMR
                      Km=500, # Michelis-Menten constant (mg/m3)
                      h=1     # Hill coeficient (1=Michaelis-Menten; h>1:cooperative binding; h<1:competitive binding)
){
  ## ANIMAL MORPHOLOGY
  M = M * 1000 # mass in mg
  dens = 1e9   # Density (mg/m3)
  V = M/dens   # Body volume
  A = a*(M/1000)^b # Surface area (m2)
  L = sqrt(A)  # Characteristic length
  
  ## WATER PROPRIETIES
  rho = 1e3 # Density of water (units for the Re number formula: kg/m3)
  Ts <- log((298.15 - Tw)/(273.15 + Tw))
  Fs = exp(S*(-0.00624523 - 0.00737614*Ts - 0.0103410*Ts^2 - 0.00817083*Ts^3) - 4.88682e-7 * S^2)
  rho_ox = 1.42905 * exp(2.00907 + 3.22014*Ts + 4.05010*Ts^2 + 4.94457*Ts^3 + 0.256847*Ts^4 + 3.88767*Ts^5) * Fs * 1e3  # Dissolved oxygen (mg/m3) Garcia and Gordon (1992)
  
  mu = (0.00000000000277388442*Tw^6-0.00000000124359703683*Tw^5+0.00000022981389243372*Tw^4-0.00002310372106867350*Tw^3+ 0.00143393546700877000*Tw^2-0.06064140920049450000*Tw+1.79157254681817000000)/1000 ## Water dyn visc (Pa s)
  D = 7.4e-8 * (Tw + 273) * (2.26 * 18)^0.5 / (mu*1000 * 25.6^0.6) * 1e-4  # diffusion coefficient (m2/s)
  hc_funct <- function(D, L, rho, v, mu) D / L * (0.6 * (L * rho * v / mu)^(1/2) * (mu / (rho * D))^(1/3)) # forced convection (Garner and Grafton, 1954)
  time_transform = 60 * 60 # sec to hours
  hc <- hc_funct(D, L, rho, v, mu) * time_transform # Convection mass transfer coefficient (m/h)
  
  ox_out <- Co * rho_ox # Oxygen concentration in water (mg/m3)
 
  ## METABOLIC RATE according to MTE
  k = 8.61e-5 # Boltzmann constant (eV/K)
  met_rate <- delta * b0 * M^(3/4) * exp(-E/(k*(273+Tw))) * time_transform # (mg/h)
  
  # Steady state
  ox_in <- -(A*hc*Km - A*hc*ox_out + met_rate - sqrt(A^2*hc^2*(ox_out^2 + Km^2) + 2*A^2*hc^2*ox_out*Km - 2*A*hc*ox_out*met_rate + 2*A*hc*met_rate*Km + met_rate^2)) / (2*A*hc)
  
  ## The actual metabolic rate is the MTE rate corrected by Michaelis-Menten kynetics
  actual_met_rate <- met_rate * ox_in^h / (ox_in^h + Km^h)
  
  return(list("CI"=ox_in, "CO"=ox_out, "hc"=hc, "met_rate"=met_rate, "actual_met_rate"=actual_met_rate))
}

########  EXAMPLE DATA

Tw <- seq(-10, 25, length.out = 100) # Water temperature range (ºC)
S <- 0 # Salinity (ppt) -freshwater
M <- 100 # Body mass (g)
b0 <- 7600 # Normalization constant
E <- 0.6 # Activation energy (eV)
delta <- 5 # Max FAS
Km <- 600 # Michaelis constant (mg m-3)
h <- 1 # Hill's coefficient
v <- 0.1 # water velocity
a <- 0.00118 # Gill surface ~ mass allometry (Pauly & Cheung, 2018)
b <- 0.79

# Resting metabolic rate (RMR) is computed for delta=1, while maximum metabolic rate (MMR) is delta > 1
out_RMR <- OxBalance(Tw, S=S, v=v, M=M, b0=b0, a=a, b=b, delta=1, Km=Km, h=h)
out_MMR <- OxBalance(Tw, S=S, v=v, M=M, b0=b0, a=a, b=b, delta=delta, Km=Km, h=h)

# Maximum metabolic rate
plot(out_MMR$actual_met_rate/M ~ Tw, ylim=c(0, 0.1), type="l", lwd=2, col="grey80", 
     ylab=expression(paste("Mass-specific metabolic rate ( ", mg[O2], g^-1, h^-1, ")")),
     xlab="Water temperature (ºC)")
# Resting metabolic rate
lines(out_RMR$actual_met_rate/M ~ Tw, lwd=2)

# Aerobic scope (MMR - RMR)
AS <- out_MMR$actual_met_rate - out_RMR$actual_met_rate

plot(AS ~ Tw,  type="l", lwd=2, 
     ylab=expression(paste("Aerobic scope ( ", mg[O2], g^-1, h^-1, ")")),
     xlab="Water temperature (ºC)")

