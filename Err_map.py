"""
Created on Wed Sep 28 12:07:10 2022

@author: Nandini Sahu
"""
import os
import numpy as np
import h5py
from astropy.io import fits

Main_dir = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/'
#Main_dir = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_Huang/'

object_name = 'DESJ2219-4348'
band = 'F140W' 
#band = 'F200LP'
camera = 'WFC3'
output='errmap'

output_data_dir = os.path.join(Main_dir, '{}/'.format(object_name)) 
cutout_dir = os.path.join(Main_dir, '{}/'.format(object_name))   

input_data_filename = 'image_{}_{}_galaxy.h5'.format(object_name,band)
#input_data_filename = 'image_{}_{}.fits'.format(object_name,band)
input_file = os.path.join(cutout_dir, input_data_filename)
 

### for .h5 file 
f = h5py.File(input_file, 'r')    
#print(f.keys())
data_reduced = f['image_data'][()]
background_rms = f['background_rms'][()]
exposure_map = f['exposure_time'][()]
ra_at_xy_0 = f['ra_at_xy_0'][()]
dec_at_xy_0 = f['dec_at_xy_0'][()]
transform_pix2angle = f['transform_pix2angle'][()]
f.close()

#print(background_rms )

### for fits
#f=fits.open(input_file)
    
error_map = np.nan_to_num( np.sqrt(data_reduced / exposure_map + background_rms**2), nan=background_rms)

#error_map_times_2=error_map*2
#error_map_times_1_8=error_map*1.8
#error_map_times_1_5=error_map*1.5
#error_map_times_1_3=error_map*1.3
#error_map_times_1_2=error_map*1.2
#error_map_by_2=error_map/2
#error_map_by_1_8=error_map/1.8
#error_map_by_1_6=error_map/1.6
#error_map_by_1_5=error_map/1.5
error_map_by_1_4=error_map/1.4
#error_map_by_1_3=error_map/1.3
#error_map_by_1_2=error_map/1.2
#error_map_by_1_1=error_map/1.1

''' saving in h5 '''
'''
out_h5_filename='{}_{}_{}.h5'.format(output,object_name,band)
out_h5file=os.path.join(output_data_dir, out_h5_filename)
h5file= h5py.File(out_h5file, 'w')
h5file.create_dataset('image_data', data=error_map)
'''


''' saving in fits '''

#h5f=h5py.File(h5file,'r')
hdu = fits.PrimaryHDU(error_map_by_1_4)
hdul = fits.HDUList([hdu])
out_fits_filename='{}_{}_{}_by_1_4.fits'.format(output,object_name,band)
out_fitsfile=os.path.join(output_data_dir, out_fits_filename)
hdul.writeto(out_fitsfile, overwrite=True)
