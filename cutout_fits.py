#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 00:38:47 2022

@author: nandinisahu
"""

import numpy as np
import os
import h5py
from astropy.io import fits
from astropy.nddata.utils import Cutout2D 
from astropy import units as u 
from astropy.coordinates import SkyCoord 
from astropy.wcs import WCS
import matplotlib.pyplot as plt

'''
#convert .h5 to .fits
input_file = 'DESJ0102+0158/image_DESJ0102+0158_F140W_galaxy.h5'
#input_file = 'DESJ0102+0158/image_DESJ0102+0158_F200LP_galaxy.h5'

f = h5py.File(input_file, 'r')
print(f.keys())
image_data = f['image_data'][()]

hdu = fits.PrimaryHDU(image_data)
hdul = fits.HDUList([hdu])

#change filename here
out_fits_file='DESJ0102+0158/image_DESJ0102+0158_F140W_galaxy.fits'
#fits_file='image_DESJ0102+0158_F200LP_galaxyfits'
hdul.writeto(out_fits_file, overwrite=True)

f.close()


'''
#read .fits file

#fits_file='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F200LP_drc_sci.fits'

fits_file='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F140W_drz_sci.fits'
#jf543j010_drc.fits

#fits_file='DCLS0854-0424/F200LP_WFC3/reduced_data/DCLS0854-0424_F200LP_drc_sci.fits'
#fits_file='DCLS0854-0424/DCLS0854-0424_F200LP_drc_sci.fits'

hdul=fits.open(fits_file)

#position1=SkyCoord(hdul[0].header['RA_APER']*u.deg,hdul[0].header['DEC_APER']*u.deg) #header does not have exact target coordinates

position=SkyCoord((226.93810653)*u.deg,(5.38230017)*u.deg)
#position=(3495.0573,3702.32)
#size=(200*u.pixel, 200*u.pixel) 
#size=(0.0083333*u.deg, 0.0083333*u.deg) 
size=(0.5*u.arcmin, 0.5*u.arcmin)

wcs1 = WCS(hdul[0].header)

cutout=Cutout2D(hdul[0].data, position, size, wcs=wcs1)
print(cutout.data)


### show
ax1 = plt.subplot(1,1,1) #projection=WCS(hdu1.header)
shw1=ax1.imshow(cutout.data, origin='lower',vmin=-0.3, vmax=0.1,cmap='RdBu_r')
plt.colorbar(shw1, shrink=1)
###

hdu2=fits.PrimaryHDU()
#hdul2=fits.HDUList([hdu2])
# Update image data from the cutout
hdu2.data = cutout.data

# Update the WCS from the cutout
hdu2.header.update(cutout.wcs.to_header())

#hdu2.header['Targ_RA']=(83.455527, '[deg] RA at centre of the lens')
#hdu2.header['Targ_DEC']=(-25.615115, '[deg] DEC at centre of the lens')
#hdu2.header['Targ_RA']=(83.4555704, '[deg] RA at centre of the lens')
#hdu2.header['Targ_DEC']=(-25.6151177, '[deg] DEC at centre of the lens')

#fits.writeto('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/DESJ0042-3718/DESJ0042-3718_F140W_cutout.fits',cutout.data,hdu2.header,overwrite=True)
#fits.writeto('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/DESJ0042-3718/DESJ0042-3718_F200LP_cutout.fits',cutout.data,hdu2.header,overwrite=True)

fits.writeto('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F140W_cutout_half_arcmin.fits',cutout.data,hdu2.header,overwrite=True)
#fits.writeto('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/DESJ0042-3718/image_DESJ0042-3718_F200LP.fits',cutout.data,hdu2.header,overwrite=True)

