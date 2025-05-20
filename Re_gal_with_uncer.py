#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 11:48:23 2023

@author: nandinisahu
"""

import pandas as pd
import numpy as np
from scipy.integrate import dblquad
from scipy.optimize import minimize_scalar
#import uncertainties as unc
#import uncertainties.unumpy as unp
#from multiprocessing import Pool

# Set the directory and import the data
df1 = pd.read_csv("/Users/nandinisahu/OneDrive - UNSW/GL_codes/AGEL_pilot_sample_Modelling_Result.csv")
#df1 = df1[df1.Name != 'AGEL105100-055628']
n = 0
print(df1['Name'][n])

sample_size=500
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

gal_tot_flux=df1['Total_Flux_count_Sersic_fit'][n]
e_gal_tot_flux=df1['e_Total_Flux_count_Sersic_fit'][n]

def f(a, n, k, re, q, x, y):
    return a * np.exp(-k * ((np.sqrt(x**2 + y**2 / q**2) / re)**(1 / n) - 1))

def h(x, y):
    return f(np.mean(amp1), np.mean(n1), np.mean(k1), np.mean(re1), np.mean(q1), x, y)

def g(x, y):
    return f(np.mean(amp2), np.mean(n2), np.mean(k2), np.mean(re2), np.mean(q2), x, y)

# Define the flux function
def flux_func(x, y):
    return h(x, y) + g(x, y)

def flux_ratio(rr):
    fluxx, _ = dblquad(flux_func, -rr, rr, lambda x: -rr / np.mean(qavg), lambda x: rr / np.mean(qavg))
    return abs(flux_tot / fluxx - 2)


# Define the function to minimize considering uncertainties
flux_tot=gal_tot_flux
# Find the value of Re_gal that satisfies the condition
result1 = minimize_scalar(flux_ratio, bounds=(0.001, 10), method='bounded')
Re_gal1 = result1.x #minimizer  (result1.fun #minimum)

flux_tot=gal_tot_flux+e_gal_tot_flux
result2 = minimize_scalar(flux_ratio, bounds=(0.001, 10), method='bounded')
Re_gal2 = result2.x

flux_tot=gal_tot_flux-e_gal_tot_flux
result3 = minimize_scalar(flux_ratio, bounds=(0.001, 10), method='bounded')
Re_gal3 = result3.x

e_Re_gal=np.abs(Re_gal3-Re_gal2)/2

# Print the value of Re_gal and its uncertainty
print("Re_gal:")
print(Re_gal1)

print("err_Re_gal:")
print(e_Re_gal)
