"""
Created on Tue Feb  7 11:48:36 2023

@author: nandinisahu
"""

from astropy.io import fits
import numpy as np

fits_file='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/PS1J2336/perturber_mask.fits'

infile=fits.open(fits_file)

img_data=infile[0].data
header=infile[0].header

#convert to npz
np.savez_compressed('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/PS1J2336/test', img_data=img_data, header=header)

#load npz
npz_file_load = np.load('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/PS1J2336/test.npz')

#checks
print(npz_file_load.files)
print(npz_file_load['img_data'])
print(npz_file_load['header'])
