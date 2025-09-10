"""
Created on Wed Sep 27 15:28:38 2023

@author: chatgpt and @nandini
"""

import numpy as np
from uncertainties import ufloat
import pandas as pd

def f(q, r, g):
    # Calculate the function value
    return ((2 / (1 + q)) ** (0.5 / g)) * (q ** 0.5) * r

def propagate_uncertainty(q, r, g, u_q, u_r, u_g):
    # Calculate the partial derivatives of f with respect to q, r, and g
    df_dq = (0.5 * (2 / (1 + q)) ** (0.5 / g) * q ** (-0.5) * r) * u_q
    df_dr = ((2 / (1 + q)) ** (0.5 / g) * (0.5 * q ** 0.5)) * u_r
    df_dg = (-(2 / (1 + q)) ** (0.5 / g) * (0.25 * q ** 0.5 / g ** 2)) * u_g

    # Calculate the total uncertainty using error propagation formula
    total_uncertainty = np.sqrt(df_dq ** 2 + df_dr ** 2 + df_dg ** 2)
    return total_uncertainty

filename='AGEL_pilot_sample_Modelling_Result.csv'
df1=pd.read_csv(filename)

n=7 #system index, starts from 0


# Input values and uncertainties
q_value = df1['q_m']           #nominal_value
r_value = df1['R_Eins']        #nominal_value
g_value = df1['gam']           #nominal_value
u_q = df1['e_q_m']             #std_dev
u_r = df1['e_R_Eins']          #std_dev
u_g = df1['e_gam']             #std_dev

# Create uncertainty objects for q, r, and g
q = ufloat(q_value[n], u_q[n])
r = ufloat(r_value[n], u_r[n])
g = ufloat(g_value[n], u_g[n])

# Calculate the function value
result = f(q, r, g)

# Propagate the uncertainty
total_uncertainty = propagate_uncertainty(q_value[n], r_value[n], g_value[n], u_q[n], u_r[n], u_g[n])

#print("spherical equivalent Einstein radius: ", result)
#print("Total Uncertainty: ", total_uncertainty)
print("spherical equivalent Einstein radius:: {:.3f}".format(result.nominal_value))
print("Total Uncertainty: {:.3f}".format(total_uncertainty))
