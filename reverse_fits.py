"""
Created on Wed Dec 14 16:58:55 2022

@author: nandinisahu
"""

#import numpy as np
from astropy.io import fits
import numpy as np

def safe_division(a, b):
    try:
        result = a / b
        return result
    except ZeroDivisionError:
        # Replace division by zero with zero
        print("Warning: Division by zero. Replacing with zero.")
        return 0

#fits_file1='/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/HST_16773/DCLS1507+0522/perturber_mask.fits'
fits_file1='/Users/nandinisahu/OneDrive - UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/GLEE_model/lensmask1.fits'

#fits_file2='/Users/nandinisahu/OneDrive - UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/GLEE_model/source2_arc_mask.fits'

#read file
file1=fits.open(fits_file1)
#file2=fits.open(fits_file2)


#manipulate
#Divide
#new_data=np.sqrt(1/((0.05**4)*file1[0].data)) #std dev from weight map
#new_data=np.sqrt(safe_division(1,(file1[0].data))) #(0.05**4)*

#subtract or reverse mask
#new_data=file1[0].data-file2[0].data
new_data=1-file1[0].data

hdu=fits.PrimaryHDU(new_data)
hdul=fits.HDUList([hdu])
hdul.writeto('/Users/nandinisahu/OneDrive - UNSW/AGEL/compound_lens/compound_lens_dcls1507-NS-main/GLEE_model/lens_mask_s2.fits')


