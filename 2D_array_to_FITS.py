#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 14:42:42 2024

@author: chatgpt and nandinisahu
"""

from astropy.io import fits
import numpy as np

# Create a 2D array (replace this with your own data)
#data = np.random.random((200, 200))

#data=np.ones((240, 240))

data = np.zeros((240, 240), dtype=int)
data[:195, :]=1   #for first lens mask
#data[196:, :]=1  #FOR SECOND LENS Mask or just reverse first one



# Specify the file name for the FITS file
fits_filename = '/Users/nandinisahu/OneDrive - UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/GLEE_model/lensmask1.fits'

# Create a FITS header (you can customize this as needed)
header = fits.Header()
header['SIMPLE'] = True
header['BITPIX'] = -64  # 64-bit floating point format
header['NAXIS'] = 2
header['NAXIS1'] = data.shape[1]
header['NAXIS2'] = data.shape[0]

# Write the data and header to the FITS file
fits.writeto(fits_filename, data, header, overwrite=True)

print(f"FITS file '{fits_filename}' has been created.")
