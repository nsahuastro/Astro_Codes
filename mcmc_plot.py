"""
Created on Thu Oct 20 18:32:46 2022

@author: nandinisahu
"""

import corner
import matplotlib.pyplot as plt
import pandas as pd
import os

#Main_dir = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_Huang/DESI-329.6820+02.9584'
Main_dir = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/DESJ2335-5152/trial4'

input_data_filename='esr_SR_DESJ2335-5152_type16_iter1_best.002_iter1.mcmc'
input_file = os.path.join(Main_dir, input_data_filename)

## first see which column number to read and plot
#
df=pd.read_table(input_file,skiprows=(1),sep="\s+")
print(df.info()) 
#df1.keys()

#cmap = plt.cm.get_cmap('gist_rainbow') #testing varying colors

### comment above and use those numbers in 'usecols' in following dataframes

'''
#Plot lens mass profile para #7,8,9,10,11,12
df1=pd.read_table(input_file, skiprows=(1),sep="\s+", usecols=[7,8,9,10,11,12,13,14])#,usecols=[0,4]
col_name=df1.columns 
#figure = corner.corner(df1,color='green',labels=col_name) 
figure = corner.corner(df1,color='green',
                       labels=[r"$x$", r"$y$",r"$b/a$", r"$\theta$", r"$\theta_E$",r"$\gamma$'",r"$\gamma_{ext}$",r"$\theta_{ext}$"],
                       title_quantiles=[0.16, 0.5, 0.84],
                       title_fmt=".3f",
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True,
                       smooth=0.5                    
                       )
'''

# doubts:
##truths={4.447, 5.366, 0.8327, 2.524, 3.4600, 0.423} ??
##quantiles=[0.16, 0.5, 0.84] #2 sigma level for the 1D histogram
##what value should smooth have ?

'''
#Plot lens light profile para 
df2=pd.read_table(input_file, skiprows=(1000),sep="\s+", usecols=[7,8,9,10,11,12,13,14,15,16,17,18,19,20])#,usecols=[0,4]
col_name=df2.columns 
figure = corner.corner(df2,color='green',
                       labels=[r"$x1$", r"$y1$", r"$q1_l$", r"$PA1_l$", r"$Amp1$",r"$R1_{e}$",r"$n1$", r"$x2$", r"$y2$", r"$q2_l$",r"$PA2_l$", r"$Amp2$",r"$R2_{e}$",r"$n2$"],
                       title_quantiles=[0.16, 0.5, 0.84],
                       title_fmt=".3f",
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True,
                       smooth=0.5                       
                       )
'''


#Plot single Sersic lens light profile para 
df2=pd.read_table(input_file, skiprows=(1),sep="\s+", usecols=[7,8,9,10,11,12,13])#,usecols=[0,4]
col_name=df2.columns 
figure = corner.corner(df2,color='green',
                       labels=[r"$x$", r"$y$", r"$q_l$", r"$PA_l$", r"$Amp$",r"$R_{e}$",r"$n$"],
                       title_quantiles=[0.16, 0.5, 0.84],
                       title_fmt=".3f",
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True,
                       smooth=0.5                       
                       )


plt.savefig('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/DESJ2335-5152/trial4/esr_SR_DESJ2335-5152_type16_iter1_best.002_iter1.pdf')
