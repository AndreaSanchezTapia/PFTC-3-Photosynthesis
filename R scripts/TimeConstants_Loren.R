######################################################################################
# Leaf thermal time constants
# Sean Michaletz sean.michaletz@gmail.com
# March 2015
######################################################################################

# Loren, make sure to run all code in section 2 (contains needed libraries)


#################################################################
### code chunk number 1: paths, libraries, load data 
#################################################################

# Clear memory
rm(list=ls(all=T))

# Load packages
library(ggplot2)
library(plyr)

# Set working directory
setwd("path")

# Load trait data
df <- read.csv("filename", header=T)

# Load climate data
climate <- read.csv("filename", header=T)

# Add climate data to trait data file
df$Elev2_m <- climate$Elev_m[ match(df$Dataset, climate$Dataset)
df$MAT <- climate$MAT[ match(df$Dataset, climate$Dataset)]
df$MAP <- climate$MAP[ match(df$Dataset, climate$Dataset)]
df$invBT_MA <- climate$invBT_MA[ match(df$Dataset, climate$Dataset)]
df$reh_MA <- climate$reh_MA[ match(df$Dataset, climate$Dataset)]
df$wnd_MA <- climate$wnd_MA[ match(df$Dataset, climate$Dataset)]
df$Lgs <- climate$Lgs[ match(df$Dataset, climate$Dataset)]
df$GST <- climate$GST[ match(df$Dataset, climate$Dataset)]
df$GSP <- climate$GSP[ match(df$Dataset, climate$Dataset)]
df$invBT_GS <- climate$invBT_GS[ match(df$Dataset, climate$Dataset)]
df$reh_GS <- climate$reh_GS[ match(df$Dataset, climate$Dataset)]
df$wnd_GS <- climate$wnd_GS[ match(df$Dataset, climate$Dataset)]
df$SOLR_GS <- climate$SOLR[ match(df$Dataset, climate$Dataset)]
df$V120cm_GS <- climate$V120cm[ match(df$Dataset, climate$Dataset)]
df$wnd120cm_MA <- climate$wnd120cm_MA[ match(df$Dataset, climate$Dataset)]
df$wnd120cm_GS <- climate$wnd120cm_GS[ match(df$Dataset, climate$Dataset)]
df$tSun_h_geosphere <- climate$GStsun[ match(df$Dataset, climate$Dataset)]
df$tSun_h_Kearney <- climate$tSUN_GS[ match(df$Dataset, climate$Dataset)]

# combine Glopnet and BIEN stomatal conductance data
df$Gs_molm2s1 <- ifelse(is.na(df$Gs_molm2s1), df$GsGlopnet_molm2s1, df$Gs_molm2s1)
df$Gs_molm2s1 <- ifelse(!is.na(df$GsGlopnet_molm2s1) & df$GsGlopnet_molm2s1 > df$Gs_molm2s1, df$GsGlopnet_molm2s1, df$Gs_molm2s1)
# Remove rows with "NA" for stomatal conductance
df <- df[!is.na(df$Gs_molm2s1),]

# remove rows without LDMC
df <- df[-which(is.na(df$LDMC_g.g)), ]


#################################################################
### code chunk number 2: libraries
### modified from code by Colin Osbourne
###
###  Library of physical constants and properties
###  (modified from c++ code for polar forest model)
###
#################################################################

# Atmospheric Pressure (Pa) at elevation (m)
# atmospheric_pressure = 101325 # standard atmosphere 
# http://en.wikipedia.org/wiki/Atmospheric_pressure
atmospheric_pressure <- function(elevation)
{
  pressure <- 101325*((1-((0.0065*elevation)/288.15))^((gravity*0.0289644)/(gas_constant*0.0065)))
  return(pressure)
}

# Air density (kg m-3), linear fit to combined data from Appendix 4 in 
# Jones (2014) and Table Table A4 from Incropera (2002)
AirDensity <- function(Tair)
{
  density = 1.29 - 0.00402 * (Tair)
  return(density)
}

## BOUNDARY LAYER CONDUCTANCE
## Given windspeed (m/s) and characteristic dimension (m), this
## returns boundary layer conductance in mm s-1 (divide by 1000 to get m s-1).
## Note that according to Jones, this includes both sides of leaf.
## Note also that conductances multiplied by 1.5 to account for outdoor 
## conditions (turbulence, etc.).  See p. 59 of Jones (2014), 
## p. 224 of Campbell and Norman (1998), and p. 231 of Bonan (2008)
## See Eqn. 3.31, p. 59 in Jones (2014)
# Checked OK, Michaletz Oct 2014
boundary <- function(geometry, windspeed, dimension)
{
  # Use flat plate for broadleaves
  if(geometry %in% "B") {
    bound_cond = 1.5 * 6.62 * (windspeed / dimension)^0.5
  } # Use cylinder for needles
  if(geometry %in% "N") {
    bound_cond = 1.5 * 4.03 * (windspeed^0.6 / dimension^0.4)
  }
  return(bound_cond)
}

# Leaf emissivity for longwave radiation (dimensionless)
#						- value from p. 27 Gates (1980)
emissivity = 0.96

# Universal Gas Constant, R	(J mol-1 K-1)
gas_constant = 8.31447

# Gravitational Acceleration	(m s-2)
gravity = 9.80665

# Latent Heat of Vaporisation for Water	(MJ g-1)
#						- temp dependency (Eq. 44, Friend 1995)
#						- temp in degrees C
#						- tested 19/05/99
# Incropera fit in J/mol: Hvap =(-0.0438*Tair+45.055)*1000
latent_heat <- function(air_temperature)
{
  latent = 3.1512 * 1E-3 - 2.38 * 1E-6 * (air_temperature + 273.15)
  if (latent > 0.002513) latent = 0.002513	# value at -5 degrees C
  
  return(latent)
}

# Psychrometer Constant	(Pa K-1)
#				- temperature dependency following Eq. 43 in Friend (1995) and P. 105 Jones (2014)
#				- temp in degrees C
#				- tested 19/05/99
psychrometer <- function(air_temperature, atmospheric_pressure, specific_heat_capacity)
{
  psych = (atmospheric_pressure * specific_heat_capacity) / (0.622 * 1E9 * latent_heat(air_temperature))
  if (psych > 68.0) psych = 68.0	# value for very high temps
  if (psych < 64.6) psych = 64.6	# value at -5 degrees C
  
  return(psych)
}

# Specific heat capacity of air (J kg-1 K-1) from Table A.4 in Incropera
# polynomial fit shows this is 1007 from 0 to around 50 C
specific_heat_capacity = 1007

# Saturation vapor pressure (Pa)
#  					- Tetens formula in Shuttleworth (1993) and cited p. 837 in New et al. (1999)
#    				- temperature in degrees C
sat_VP <- function(air_temp)
{
  saturation_VP = 610.8 * exp((17.27 * air_temp) / (237.3 + air_temp))
  
  # trap very low values of sat_VP
  # (this would need looking at if you're interested in low temps)
  if (saturation_VP < 402) saturation_VP = 402 # sat_VP over ice at -5 degrees C
  
  return(saturation_VP)
}

# Saturation Vapour Pressure  - Pa
#					- temp dependence from Buck (1981)
#					- temperature in degrees C
#					- tested 19/05/99
#sat_VP <- function(air_temp)
#{
#  sat_VP = 613.75 * exp((17.502 * air_temp) / (240.97 + air_temp))
#  if (sat_VP < 402.) sat_VP = 402. # sat_VP over ice at -5 degrees C
#  
#  return(sat_VP)
#}

# Slope of the Relationship between Saturation VP and Temperature (Pa K-1)
# This is exponential fit to data from Appendix 4 in Jones (2014) (actually s=46.433*e^(0.0546*Ta))
# Friend (1997): the slope of the saturation vapour concentration versus
# temperature is calculated between air temperature and air temperature 
# plus 0.1 K, rather than between air and foliage temperature (Eq. (53) of Friend, 1995). 
# This last approximation is necessary because no iteration to solve the energy balance 
# is used, yet this slope needs to be known to calculate foliage temperature.
slope_sat_VP <- function(air_temperature)
{
  slope = 48.7 * exp(0.0532 * air_temperature)
  if (slope < 32) slope = 32	# slope at -5 degrees C
  
  return(slope)
}

# Stefan-Boltzman Constant (W m-2 K-4)
#				- used to calculate longwave radiation emission
StefBoltz = 5.670373 * 1E-8  


#################################################################
### code chunk number 3: Leaf trait calculations
#################################################################

# Absorptivity
# Gates (1980) uses 0.7 for conifers and 0.6 for deciduous
df$Absorptance <- ifelse(df$LeafType %in% "N", 0.7, 0.6) # this requires LeafType column, "N" for needle, "B" for braodleaf

# Total leaf area (m2)
# 3.14 * projected area for needle leaves, 2 * projected area for broadleaves
# (Chapin et al. Principles of Terrestrial Ecosystem Ecology)
#df$Arealeaf_m2 <- ifelse(df$LeafType %in% "N", pi * df$ProjAreaLeaf_m2, 2 * df$ProjAreaLeaf_m2)

# Characteristic dimension (m)
# Use width when known. When unknown, we can assume shape is 
# square and take L = 0.7*sqrt(area) based on Table 7.5
# in Campbell and Norman (1998), although this will result 
# in artificially large L for species with compound leaves
# such as "code" column 2469 (Didymopanax morototoni).
df$Dimension_m <- 0.7 * df$ProjAreaLeaf_m2^0.5 # This avoids using TRY and doesn't change results

# Add phi conversion factor for broadleaves and needle-leaves
df$phi <- ifelse(df$LeafType %in% "N", 1/pi, 1/2) # this requires LeafType column, "N" for needle, "B" for braodleaf


#################################################################
### code chunk number 4: Climate calculations
#################################################################

# Calculate VPD (Pa), P_atm (Pa), and air density (kg m-3) using functions defined in physical library
for(i in 1:nrow(df))
{
  # calculate saturation vapour pressure based on air temperature and add to data frame
  df$SatVP[i] <- sat_VP(df$GST[i])
  
  # calculate VPD and add to the data frame using Eq. 2.15 from 
  # Landsberg and Sands (2010) VPD = sat_VP * (100-RH)/100
  df$VPD[i] <- sat_VP(df$GST[i]) * (1-(df$reh_GS[i]/100))
  
  # calculate atmospheric pressure
  df$Patm_Pa[i] <- atmospheric_pressure(df$Elev2_m[i])
  
  # calculate air density
  df$AirDens_kgm3[i] <- AirDensity(df$GST[i])
  
  # Psychrometer Constant  (Pa K-1)
  df$psych_PaK[i] <- psychrometer(df$GST[i], df$Patm_Pa[i], specific_heat_capacity)
  
  #Slope of the Relationship between Saturation VP and Temperature (Pa K-1)
  df$slope[i] <- slope_sat_VP(df$GST[i])
}

# Set all wind velocities to 0.3 m/s
df$U_ms <- 0.3 # this can be treated as variable, contant, or mean value from a site


#################################################################
### code chunk number 5: Calculate conductances
#################################################################
# define a function to convert gs from molar to velocity units (m s-1)
# Eqn. 3.23a. See section 3.2.6 on p. 52 in Jones (2014) for explanation
convert_gs <- function(molar_gs, airTemp, AtmPressure)
{
  if (molar_gs <= 0.001) gs = 0
  else gs = molar_gs * gas_constant * (airTemp + 273.15) / AtmPressure
  return(gs)
}

# calculate combined boundary layer and stomatal conductances to water and heat
# add this information to every line in the data frame
for(i in 1:nrow(df))
{
  # convert gs from molar to velocity units (m s-1)
  df$aero_gs[i] <- convert_gs(df$Gs_molm2s1[i], df$GST[i], df$Patm_Pa[i])
  
  # calculate leaf boundary layer conductance (mm s-1)
  df$ga[i] <- boundary(df$LeafType[i], df$U_ms[i], df$Dimension_m[i])
  
  # convert boundary layer conductance from mm s-1 to m s-1
  df$ga[i] <- df$ga[i] / 1000
  
  # calculate boundary layer conductance to heat (m s-1)
  df$gH[i] <- df$ga[i]
  
  # calculate radiation conductance (m s-1); p. 10 Incropera, p. 101 Jones.
  df$gR[i] <- (4*emissivity*StefBoltz*(df$GST[i] + 273.15)^3)/(df$AirDens_kgm3[i] * specific_heat_capacity)
  
  # calculate total conductance to water (m s-1) = boundary layer + stomatal in series 
  # This implicitly means ga from both sides. Measured aero_gs must also
  # implicitly include amphi/hypostomatous differences.
  # Eq. in paragraph on top right of p. 105 in Jones (2014)
  # Eqs. 7.13 and 7.14, p. 146 Monson and Baldocchi 
  df$gW[i] <- df$ga[i] * df$aero_gs[i] / (df$ga[i] + df$aero_gs[i])
}
#write.csv(df, "CodeCheck.csv")

#################################################################
### code chunk number 6: Calculate site-specific h and t
#################################################################
# 1) Convection only, h = rho*c_p*g_H; Eq A.2 from Campbell (1977)
# Heat transfer coefficient (W m-2 C-2)
df$h1_Wm2 <- df$AirDens_kgm3 * specific_heat_capacity * df$gH
# Thermal tims constant (s)
df$t1_s <- df$phi * df$LMA_kgm2 * ((4180 / (df$LDMC_g.g * df$h1_Wm2)) + ((1902.6 - 4180) / df$h1_Wm2))

# 2) Convection and radiation, h = rho*c_p*(g_H + g_R); Eq 15.9 from Monteith and Unsworth (2013)
# Heat transfer coefficient (W m-2 C-2)
df$h2_Wm2 <- df$AirDens_kgm3 * specific_heat_capacity * (df$gH + df$gR)
# Thermal tims constant (s)
df$t2_s <- df$phi * df$LMA_kgm2 * ((4180 / (df$LDMC_g.g * df$h2_Wm2)) + ((1902.6 - 4180) / df$h2_Wm2))

# 3) Convection, radiation, and transpiration, h = rho*c_p*(g_H + g_R + g_W*s/gamma)
df$h3_Wm2 <- df$AirDens_kgm3 * specific_heat_capacity * (df$gH + df$gR + df$gW * df$slope / df$psych_PaK)
# Thermal tims constant (s)
df$t3_s <- df$phi * df$LMA_kgm2 * ((4180 / (df$LDMC_g.g * df$h3_Wm2)) + ((1902.6 - 4180) / df$h3_Wm2))

