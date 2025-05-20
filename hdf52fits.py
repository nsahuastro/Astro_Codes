"""
Created on Wed Sep 28 10:30:20 2022

@author: Nandini Sahu
"""

#ßimport numpy as np
import os
import h5py
from astropy.io import fits

#Main_dir = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/'
#Main_dir = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_Huang/'
Main_dir='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/GLEE_model/'

object_name = 'DCLS1507+0522' #psf_model #objnameß
band = '200LP' 
#band = 'F200LP'
camera = 'WFC3'

#output_data_dir = os.path.join(Main_dir, '{}/'.format(object_name)) 
#cutout_dir = os.path.join(Main_dir, '{}/'.format(object_name))   

output_data_dir=Main_dir
cutout_dir=Main_dir
 
##output_data_dir = os.path.join(Main_dir, '{}/{}_{}/reduced_data/'.format(object_name,band,camera))    
##cutout_dir = os.path.join(Main_dir, 'cutout/{}/'.format(object_name,band,camera))   

#input_data_filename = 'image_{}_{}.h5'.format(object_name,band)
input_data_filename = 'psf_model_F200LP.h5'


input_file = os.path.join(cutout_dir, input_data_filename)
#input_file = os.path.join(cutout_dir,'errmap_DESJ2335-5152_F140W.h5')

f = h5py.File(input_file, 'r')
print(f.keys())
#image_data = f['image_data'][()]          #for image
psf_data = f['kernel_point_source'][()]    #for psf


#hdu = fits.PrimaryHDU(image_data)
hdu = fits.PrimaryHDU(psf_data)
hdul = fits.HDUList([hdu])
#out_fits_filename='image_{}_{}.fits'.format(object_name,band)
out_fits_filename='psf_model_F200LP.fits'
out_file=os.path.join(output_data_dir, out_fits_filename)
hdul.writeto(out_file)

#f.close()