#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 22:13:45 2022

@author: jonah
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 21:50:03 2022

@author: jonah
"""

# The previous script was a reusable homogenizaition script
## Now let's take a stab at writng this as a re-usable REPT Function

def REPT(mat,den,triso,pellet,rept):
    """ This script is used to homogenize a TRISO fuel pellet to 
    simplify Monte Carlo modeling. 
    
    Inputs:
        mat = [abbreviation for fuel type, enrichment, fuel pack frac]
        den = [Array of material densities]
        triso = [outer radii of each TRISO layer]
        pellet = [fuel pellet dimensions]
        rept = radius to limit fuel to
    
    Output:
        Weight fractions for single material definition
    """
    import numpy as np
    
    # Enrichment
    u_mol = (235*mat[1]+238*(1-mat[1]));
    
    # Various Fuel Properties
    u_fraction = {
            "un" : (u_mol)/(u_mol+14.007),
            "uco": (u_mol)/(u_mol+15.999+14.007),
            "uo2": (u_mol)/(u_mol+15.999*2)}
    
    density = {
            "un" : 11.3,
            "uco": 13.63,
            "uo2": 10.97} #Density g/cc
    
    #Material Densities
    [buf_rho, ipyc_rho, sic_rho, opyc_rho, matrix_rho] = den;
    
    #Triso Dimensions
    
    [f_or, buf_or, ipyc_or, sic_or, opyc_or, tr_or] = triso;
    
    
    #Particle Volumes
    
    kn_V = 4/3*np.pi*f_or**3; #Fuel Kernel Volume [cc]
    buf_V = 4/3*np.pi*(buf_or**3-f_or**3); # Buffer Volume [cc]
    ipyc_V = 4/3*np.pi*(ipyc_or**3-buf_or**3); #Inner PyC Volume [cc]
    sic_V = 4/3*np.pi*(sic_or**3-ipyc_or**3); #SiC Volume [cc]
    opyc_V = 4/3*np.pi*(opyc_or**3-sic_or**3); # Outer PyC Volume
    tr_V =  4/3*np.pi*tr_or**3; #TRISO Particle Volume [cc]
    
    #Fuel Rod Attributes
    pl_or, pl_len = pellet;
    pl_V = np.pi*pl_or**2*pl_len; #Pellet Volume [cc]
    
    #Fuel System Design 
    pf = mat[2]; #Packing Fraction
    num_P = round(pl_V*pf/tr_V) #Number of Particles in Defined Pellet Volume
    tot_f_V = num_P*kn_V; #Total Fuel Kernel Volume in Defined Pellet [cc]


    # Properties for Homogenization

    tot_U_g =u_fraction[mat[0]]*tot_f_V*density[mat[0]]; #Total gU in fuel rod
    tot_U235_g = mat[1]*tot_U_g; #Total gU235 in the rod
    tot_U238_g = (1-mat[1])*tot_U_g; #Total gU238 in the rod 
    
    #Mass of other Particle Materials
    tot_buf_g = buf_V*num_P*buf_rho;
    tot_pyc_g = ipyc_V*num_P*ipyc_rho+opyc_V*num_P*opyc_rho;
    tot_sic_g = sic_V*num_P*sic_rho;
    
    #Establish two matrix regions
    tot_mat_g = pl_V*(1-pf)*matrix_rho;
    reg1_frac = (rept**2/pellet[0]**2);
    reg1_mat_g = tot_mat_g*reg1_frac;
    
    #Mass of Fuel Pellet where fuel exists
    reg1_pel_g = reg1_mat_g+tot_U_g+tot_buf_g+tot_sic_g+tot_pyc_g;
    
    #Pellet Fractions
    reg1_mass_fraction=np.divide([tot_U235_g,tot_U238_g,tot_buf_g,tot_pyc_g,tot_sic_g,reg1_mat_g],reg1_pel_g);
    
    return reg1_mass_fraction


# Define Test Case
mat1 = ['un', 0.195,.34];
den1= [1, 1.9, 3.2, 1.8, 1.75]; 
triso1= [0.02135, 0.03135, 0.03485, 0.03835, 0.04325, 0.04235] 
pellet1=[0.6,10];
rept1=0.3

#Run Test Case
print(REPT(mat1,den1,triso1,pellet1,rept1))