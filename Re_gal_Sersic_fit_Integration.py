#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 14:39:42 2023

@author: nandinisahu
"""

import pandas as pd
import numpy as np
from scipy.integrate import dblquad
from scipy.optimize import minimize_scalar

# Set the directory and import the data
df1 = pd.read_csv("/Users/nandinisahu/OneDrive - UNSW/GL_codes/AGEL_pilot_sample_Modelling_Result.csv")
df1=df1[df1.Name != 'AGEL105100-055628']
n=4

# Extract relevant parameters
q1 = df1['q1_L'][n]
pa1 = df1['PA1_L'][n]
amp1pixel = df1['amp_L1'][n]
amp1 = df1['amp_L1'][n]/(0.08 * 0.08)
re1 = df1['R_sersic1'][n]
n1 = df1['n_sersic1'][n]
k1 = 2 * n1 - (1 / 3) + 4 / (n1 * 405) + 46 / (25515 * (n1 ** 2))

q2 = df1['q2_L'][n]
pa2 = df1['PA2_L'][n]
amp2pixel = df1['amp_L2'][n]
amp2 = df1['amp_L2'][n]/(0.08 * 0.08)
re2 = df1['R_sersic2'][n]
n2 = df1['n_sersic2'][n]
k2 = 2 * n2 - (1 / 3) + 4 / (n2 * 405) + 46 / (25515 * (n2 ** 2))

qavg = (q1 + q2) / 2

# in following change q1 to which ever q values is the dominent component

# Define the functions f, h, and g
def f(a, n, k, re, q, x, y):
    return a * np.exp(-k * ((np.sqrt(x**2 + y**2 / q**2) / re)**(1 / n) - 1))

def h(x, y):
    return f(amp1, n1, k1, re1, q1, x, y)

def g(x, y):
    return f(amp2, n2, k2, re2, q2, x, y)

# Define the flux function
def flux_func(x, y):
    return h(x, y) + g(x, y)

# Perform the numerical integration for complete flux
flux, _ = dblquad(flux_func, -1000, 1000, lambda x: -1000 / q1, lambda x: 1000 / q1) #qavg

print('flux count e-/s')
print(flux)


# Define the function to minimize
def flux_ratio(rr):
    fluxx, _ = dblquad(flux_func, -rr, rr, lambda x: -rr / q1, lambda x: rr / q1)
    return abs(flux / fluxx - 2)

# Find the value of Re_gal that satisfies the condition
result = minimize_scalar(flux_ratio, bounds=(0.001, 10), method='bounded')
reall = result.x

# Print the value of Re_gal
print("Re_gal:")
print(reall)

# Print the fluxx
#print("Fluxx:")
#print(fluxx)
