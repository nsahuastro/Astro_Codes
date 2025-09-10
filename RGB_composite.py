"""
Created on Fri Jun  2 13:41:38 2023

@author: nandinisahu
"""

#import pyfits
from astropy.io import fits
import numpy as np
#import pylab as py
import matplotlib.pyplot as plt
#from astropy.visualization import make_lupton_rgb
#import img_scale

 
fits_red_img ='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F140W_cutout_half_arcmin.fits'
fits_green_img ='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_green.fits'
fits_blue_img = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_F200LP_cutout_half_arcmin_scaled.fits'

red_img=fits.open(fits_red_img)[0]
green_img=fits.open(fits_green_img)[0]
blue_img=fits.open(fits_blue_img)[0]

#
img = np.zeros((red_img.shape[0], red_img.shape[1], 3), dtype=float)
img[:,:,0] = np.sqrt(red_img.data)/2 #red_img.data/2 #np.sqrt(red_img.data/10)
img[:,:,1] = np.sqrt(green_img.data)/3 #green_img.data/2 #np.sqrt(green_img.data/10)
img[:,:,2] = 15*blue_img.data


#image = make_lupton_rgb(red_img.data, green_img.data, 10*blue_img.data, Q=10, stretch=0.5)
#plt.imshow(image)

plt.imshow(img,origin='lower',vmin=-1, vmax=1,aspect='equal') #,cmap='Spectral'
#plt.title('DCLS1507')
plt.xticks([])  # Remove x-axis ticks
plt.yticks([])  # Remove y-axis ticks

plt.savefig('/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/data/DCLS1507+0522_rgb.pdf',dpi=1200)
