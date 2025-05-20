#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 21:23:50 2022

@author: nandinisahu
"""
#import os
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from reproject import reproject_interp

#Main_dir ='../Vy_DCLS0854-0424/'
#input_file1 = os.path.join(Main_dir, 'DCLS0854-0424_F140W_drz_sci_cutout.fits')
#input_file2 = os.path.join(Main_dir, 'DCLS0854-0424_F200LP_drc_sci_cutout.fits')

fits_file1='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F140W_cutout_half_arcmin.fits'
fits_file2='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F200LP_cutout_half_arcmin.fits'

hdu1=fits.open(fits_file1)[0]
hdu2=fits.open(fits_file2)[0]
## note:  get_pkg_data_filename only works on local files so above PATH didnt work here
#hdu1 = fits.open(get_pkg_data_filename('DCLS0854-0424/DCLS0854-0424_F140W_cutout.fits'))[0]
#hdu2 = fits.open(get_pkg_data_filename('DCLS0854-0424/DCLS0854-0424_F200LP_cutout.fits'))[0]

'''
ax1 = plt.subplot(1,2,1, projection=WCS(hdu1.header))
ax1.imshow(hdu1.data, origin='lower',vmin=0., vmax=1)
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')
ax1.set_title('HST F140W')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu2.header))
ax2.imshow(hdu2.data, origin='lower', vmin=0, vmax=1)
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
ax2.coords['dec'].set_axislabel_position('r')
ax2.coords['dec'].set_ticklabel_position('r')
ax2.set_title('HST F200LP')
'''
array, footprint = reproject_interp(hdu2, hdu1.header)

ax1 = plt.subplot(1,2,1, projection=WCS(hdu1.header))
ax1.imshow(array, origin='lower', vmin=-0.1, vmax=0.1)
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')
ax1.set_title('Reprojected HST F200LP image')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu1.header))
ax2.imshow(footprint, origin='lower', vmin=0, vmax=1.5)
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
ax2.coords['dec'].set_axislabel_position('r')
ax2.coords['dec'].set_ticklabel_position('r')
ax2.set_title('HST F200LP image footprint')

fits.writeto('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F200LP_cutout_half_arcmin_scaled.fits', array, hdu1.header, overwrite=True)

