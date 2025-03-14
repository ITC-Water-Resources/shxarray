# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#

a_earth=0.6378136460e+07 #mean radius in meter
rho_sea=1.025e3 #average density of sea water kg/m^3
rho_water=1.e3 # average density of water
rho_earth=5517.0 #average density of the Earth
rho_ice=931.0 #density of ice kg/m^3 taken from G. Spada and friends
g=9.80665e0 # mean gravity m/s^2



#flattening accroding to GRS80
grs80_C20=-4.84166854896120e-04

#Mean moments of interita from the earth ( kg m^2)
ixx_A=8.0102e+37  #~= (I_xx+Iyy)/2
izz_C=8.0365e+37

#average earth rotation rate in rad/s
ohm_earth=7.292115467064e-5

#Chandler frequency ( radians per second) from iers constants
ohm_chandler=1.679064144e-7

#some standard degree 2 elastic body and load Love numbers
k2loadprem=-0.3054020195e+00
k2bodyprem=0.303
l2bodyprem=0.0855
h2bodyprem=0.612
