#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:06:06 2023

@author: nandinisahu
"""

import pandas as pd
import numpy as np
from scipy.integrate import dblquad
from scipy.optimize import minimize_scalar
#import uncertainties as unc
#import uncertainties.unumpy as unp
from multiprocessing import Pool

# Set the directory and import the data
df1 = pd.read_csv("/Users/nandinisahu/OneDrive - UNSW/GL_codes/AGEL_pilot_sample_Modelling_Result.csv")
#df1 = df1[df1.Name != 'AGEL105100-055628']
n = 7
print(df1['Name'][n])
PHOTFLAM= df1['PHOTFLAM (ergs/cm2/A/e-)'][n] #Inverse sensitivity, ergs/cm2 A/e-
STmag_ZP=21.10 #STmag zeropoint magnitude to use here
PHOTPLAM=df1['PHOTPLAM (Angstroms)'][n]

sample_size=10 #500

#Sersic1
q1=np.random.normal(df1['q1_L'][n], df1['e_q1_L'][n], sample_size)
pa1 = np.random.normal(df1['PA1_L'][n], df1['e_PA1_L'][n], sample_size)
amp1pixel = np.random.normal(df1['amp_L1'][n], df1['e_amp_L1'][n], sample_size)  # Assuming no uncertainty for amp1pixel
amp1 = amp1pixel / (0.08 * 0.08)  # No uncertainty in conversion
re1 = np.random.normal(df1['R_sersic1'][n], df1['e_R_sersic1'][n], sample_size)
n1 = np.random.normal(df1['n_sersic1'][n], df1['e_n_sersic1'][n], sample_size)
k1 = 2 * n1 - (1 / 3) + 4 / (n1 * 405) + 46 / (25515 * (n1 ** 2))  

# Sersic2
q2 = np.random.normal(df1['q2_L'][n], df1['e_q2_L'][n], sample_size)
pa2 = np.random.normal(df1['PA2_L'][n], df1['e_PA2_L'][n], sample_size)
amp2pixel = np.random.normal(df1['amp_L2'][n], df1['e_amp_L2'][n], sample_size)  # Assuming no uncertainty for amp2pixel
amp2 = amp2pixel / (0.08 * 0.08)  # No uncertainty in conversion
re2 = np.random.normal(df1['R_sersic2'][n], df1['e_R_sersic2'][n], sample_size)
n2 = np.random.normal(df1['n_sersic2'][n], df1['e_n_sersic2'][n], sample_size)
k2 = 2 * n2 - (1 / 3) + 4 / (n2 * 405) + 46 / (25515 * (n2 ** 2))  

qavg = (q1 + q2) / 2

def f(a, n, k, re, q, x, y):
    return a * np.exp(-k * ((np.sqrt(x**2 + y**2 / q**2) / re)**(1 / n) - 1))

def total_flux_gen(i):
    amp11 = amp1[i]
    n11=n1[i]
    k11=k1[i]
    re11=re1[i]
    q11=q1[i]
    amp22=amp2[i]
    n22=n2[i]
    k22=k2[i]
    re22=re2[i]
    q22=q2[i]
    q_avg=qavg[i]
    flux, _ = dblquad(lambda x, y: f(amp11, n11, k11, re11, q11, x, y)+f(amp22, n22, k22, re22, q22, x, y), -1000, 1000, lambda x: -1000 / q_avg, lambda x: 1000 / q_avg)
    return flux


def main():
    i=range(0,sample_size,1)
    p=Pool(7)
    result=p.map(total_flux_gen, i)
    return result

if __name__== "__main__":
    count = main()
    gal_tot_flux=np.mean(count)
    err_gal_tot_flux=np.std(count)
    print('total flux (e-/s):', gal_tot_flux)
    print('uncertainty in total flux:',err_gal_tot_flux)
    Lens_mag_STmag=-2.5*np.log10(gal_tot_flux*PHOTFLAM)-STmag_ZP  #flux units e-/s        
    Lens_mag=Lens_mag_STmag -5*np.log10(PHOTPLAM)+18.692 #AB mag
    Lens_mag_plus=-2.5*np.log10((gal_tot_flux+err_gal_tot_flux)*PHOTFLAM)-STmag_ZP -5*np.log10(PHOTPLAM)+18.692
    Lens_mag_minus=-2.5*np.log10((gal_tot_flux-err_gal_tot_flux)*PHOTFLAM)-STmag_ZP -5*np.log10(PHOTPLAM)+18.692
    e_Lens_mag=(Lens_mag_minus-Lens_mag_plus)/2
    print('Mag_gal_AB:',Lens_mag)
    print('err_Mag_gal_AB:', e_Lens_mag)
