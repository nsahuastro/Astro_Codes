#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 00:58:58 2022

@author: nandinisahu
"""

import numpy as np
from astropy.io import fits
#from astropy.wcs import WCS
fits_file1='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F140W_cutout_half_arcmin.fits'
fits_file2='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F200LP_cutout_half_arcmin_scaled.fits'

file1=fits.open(fits_file1)
file2=fits.open(fits_file2)

dataf1=file1[0].data
dataf2=file2[0].data

green_image= np.sqrt(dataf1**2 + dataf2**2)

hdu=fits.PrimaryHDU(green_image)
hdul=fits.HDUList([hdu])
hdul.writeto('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_green.fits')