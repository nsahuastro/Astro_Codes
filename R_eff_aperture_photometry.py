"""
Created on Tue Oct 18 10:19:22 2022

@author: nandinisahu
"""

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  
from matplotlib.patches import Ellipse
import numpy as np
import sep
#import cma
#from reproject import reproject_interp

#fits_file1='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_Huang/DESI-329.6820+02.9584/SR_DESI-329.6820+02.9584_typ16_iter1_best_es001_im.fits'
#fits_file_full='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_Huang/DESI-329.6820+02.9584/DESI-329.6820+02.9584_F140W_drz_sci.fits'
fits_file1='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/DESJ0537-4647/trial3/SR_DESJ0537-4647_t2_type16_best_es001_im.fits'
fits_file_full='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/DESJ0537-4647/DESJ0537-4647_F140W_drz_sci.fits'

hdu1=fits.open(fits_file1)[0]
#hdu2=fits.open(fits_file2)[0]
hdu_full=fits.open(fits_file_full)[0]

arc_grid=hdu1.data[0] # can choose either hdu2.data[0] or hdu2.data[1], difference can be seen from residual 
arc_grid_header=hdu1.header
#src_grid=hdu2.data
#src_grid_header=hdu2.header
full_file_header=hdu_full.header

Observed=hdu1.data[5]
Reconstructed=hdu1.data[6]
Lens_light=hdu1.data[7]
#Lens_light2=hdu1.data[5]-hdu1.data[1]
Normalized_Residual=hdu1.data[9]
#Source_light=hdu2.data


#----------------------------Aperture photometry Using sep package
Lens_light_sep=Lens_light.byteswap(inplace=True).newbyteorder()
#Lens_light_sep=Lens_light2
#Lens_light_sep=Lens_light.byteswap().newbyteorder()
##Input array with dtype '>f8' has non-native byte order. Only native byte order arrays are supported. To change the byte order of the array 'data', do 'data = data.byteswap().newbyteorder()'

bkg = sep.Background(Lens_light_sep)  #Background object

bkg_2D=bkg.back() #Background as 2D array of the size of 'Lens_light' grid.
bkg_rms_2D = bkg.rms() #evaluate the background noise as 2-d array, same size as original image
bkg_mean=bkg.globalback # get a "global" mean of the image background
bkg_rms=bkg.globalrms #get a "global" noise of the image background

#plt.imshow(bkg_rms_2D, cmap='gray', origin='lower') #, interpolation='nearest'
#plt.colorbar()

#--------------objects detection  (could be single galaxy or multiple)

# subtract the background IF REQUIRED
#Lens_light_sep_sub = Lens_light_sep - bkg 
###data_sub = bkg.subfrom(Lens_light.byteswap().newbyteorder())
f_bkg=0  #should be some non-zero value if there is background noise
objects = sep.extract(Lens_light_sep, f_bkg, err=bkg.globalrms) #object detection #bkg.globalrms
objects.dtype.names  ## available fields
lens_x=objects['x']
lens_y=objects['y']


#plot objected whether background subtracted or not
fig, ax = plt.subplots()
m, s = np.mean(Lens_light_sep), np.std(Lens_light_sep)
im=ax.imshow(Lens_light_sep, cmap='gray', origin='lower', vmin=m-s, vmax=m+s)

f=6 #factor deciding ellipse size
#f1=0.037976993223888804
#f2=0.04375941635638057

# plot an ellipse for each object
for i in range(len(objects)):
    e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                width=f*objects['a'][i],
                height=f*objects['b'][i],
                angle=objects['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)

#--------------Flux calculation

#aperture photometry with a 'R' pixel radius at the locations of the objects:

#circle
#R=6*objects['a'] #a multiple of objects semi major axis
#flux, fluxerr, flag = sep.sum_circle(Lens_light_sep, objects['x'], objects['y'],R, err=bkg.globalrms, gain=1.0)

#ellipse
#gain adds poisson noise to flux error. use 1 or None or some other well thought factor
flux, fluxerr, flag = sep.sum_ellipse(Lens_light_sep, objects['x'], objects['y'], f*objects['a'],f*objects['b'],objects['theta'], err=bkg.globalrms, gain=None)


# show the objects results:
for i in range(len(objects)):
    print("object {:d}: flux = {:f} +/- {:f}".format(i, flux[i], fluxerr[i]))

    
#Lens_mag=-2.5*np.log10(flux/(10**(ZP/-2.5)))
#ref: https://hst-docs.stsci.edu/wfc3dhb/chapter-9-wfc3-data-analysis/9-1-photometry#:~:text=Formally%2C%20the%20HST%20VEGAMAG%20system,(Fobject%2FFvega)
PHOTFLAM= full_file_header['PHOTFLAM'] #  #Inverse sensitivity, ergs/cm2/A/e-
#EXP_time= full_file_header['EXPTIME']  #exposure time seconds
STmag_ZP=21.10 #STmag zeropoint magnitude
Lens_mag_STmag=-2.5*np.log10(flux*PHOTFLAM)-STmag_ZP  #flux units e-/s
#AB
PHOTPLAM=full_file_header['PHOTPLAM']
Lens_mag=Lens_mag_STmag -5*np.log10(PHOTPLAM)+18.692
for i in range(len(objects)):
    print("object {:d}: Lens apparant mag = {:f} mag".format(i, Lens_mag[i]))

#--------------half light radius calculation
r, flag = sep.flux_radius(Lens_light_sep, objects['x'], objects['y'], f*objects['a'], 0.5,  subpix=5) #normflux=flux,

q_L=objects['b']/objects['a'] #axis ratio

for i in range(len(objects)):
    print("object {:d}: axis_ratio={:f}, PA={:f} rad, R_eff = {:f} pix".format(i, q_L[i], objects['theta'][i], r[i]))
