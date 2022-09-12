# This script is my scratchpad for developing an REPT alogrithm #

#First let's just make a homogenization script

import numpy as np

fuel_type = 'un';

#Material Densities
buf_rho = 1; #g/cc
ipyc_rho = 1.9; #g/cc
sic_rho = 3.2;  #g/cc
opyc_rho = 1.87; #g/cc
matrix_rho = 1.75; #g/cc  - graphite at 1000K


# Enrichment
enr = 0.195; #U235 Enrichment
u_mol = (235*enr+238*(1-enr))

# Add in coatings #
f_or = 0.02135; # TRISO Fuel Kernel Outer Radius [cm]
buf_or = 0.03135; # Buffer Outer Radius [cm]
ipyc_or = 0.03485; #Inner PyC Radius [cm]
sic_or = 0.03835; #SiC Radius
opyc_or = 0.04325; #Outer PyC Radius
tr_or = 0.04235; # TRISO Particle Outer Radius [cm]

#Particle Volumes - From Ramey and Petrovic

kn_V = 4/3*np.pi*f_or**3; #Fuel Kernel Volume [cc]
buf_V = 4/3*np.pi*(buf_or**3-f_or**3); # Buffer Volume [cc]
ipyc_V = 4/3*np.pi*(ipyc_or**3-buf_or**3); #Inner PyC Volume [cc]
sic_V = 4/3*np.pi*(sic_or**3-ipyc_or**3); #SiC Volume [cc]
opyc_V = 4/3*np.pi*(opyc_or**3-sic_or**3); # Outer PyC Volume
tr_V =  4/3*np.pi*tr_or**3; #TRISO Particle Volume [cc]

# Various Fuel Properties
u_fraction = {
        "un" : (u_mol)/(u_mol+14.007),
        "uco": (u_mol)/(u_mol+15.999+14.007),
        "uo2": (u_mol)/(u_mol+15.999*2)}

density = {
        "un" : 11.3,
        "uco": 13.63,
        "uo2": 10.97} #Density g/cc

#Fuel Rod Attributes
pl_or = .6; #FCM Pellet Outer Radius [cm]
pl_len = 10; #FCM Pellet Length [cm]
pl_V = np.pi*pl_or**2*pl_len; #Pellet Volume [cc]

#Fuel System Design 
pf = 0.34; #Packing Fraction
num_P = round(pl_V*pf/tr_V) #Number of Particles in Defined Pellet Volume
tot_f_V = num_P*kn_V; #Total Fuel Kernel Volume in Defined Pellet [cc]

# Properties for Homogenization

tot_U_g =u_fraction[fuel_type]*tot_f_V*density[fuel_type]; #Total gU in fuel rod
tot_U235_g = enr*tot_U_g; #Total gU235 in the rod
tot_U238_g = (1-enr)*tot_U_g; #Total gU238 in the rod 

#Mass of other Particle Materials
tot_buf_g = buf_V*num_P*buf_rho;
tot_pyc_g = ipyc_V*num_P*ipyc_rho+opyc_V*num_P*opyc_rho;
tot_sic_g = sic_V*num_P*sic_rho;

#Mass of non-particle constituents
tot_mat_g = pl_V*(1-pf)*matrix_rho;

#Mass of Fuel Pellet
tot_pel_g = tot_mat_g+tot_U_g+tot_buf_g+tot_sic_g+tot_pyc_g;

#Pellet Fractions
pel_mass_fraction=np.divide([tot_U235_g,tot_U238_g,tot_buf_g,tot_pyc_g,tot_sic_g,tot_mat_g],tot_pel_g);



