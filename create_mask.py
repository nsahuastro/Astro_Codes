"""
Created on Tue Feb 27 14:25:40 2024

@author: nandinisahu
"""

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import CCDData
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.visualization import ZScaleInterval
import numpy as np
import pyregion

def create_mask_fits(ds9_regions_file, output_fits_file, image_data):
    # Load DS9 regions file using pyregion
    with open(ds9_regions_file, 'r') as f:
        region_string = f.read()
    ds9_regions = pyregion.parse(region_string)

    # Create an empty mask
    mask = np.zeros_like(image_data, dtype=bool)

    # Create mask based on DS9 regions
    for region in ds9_regions:
        if region.name == 'circle':
            # Get the circle parameters
            center = region.coord_list[:2]
            radius = region.coord_list[2]

            # Create a mask for the circle
            y, x = np.ogrid[:image_data.shape[0], :image_data.shape[1]]
            mask_circle = ((x - center[0]) ** 2 + (y - center[1]) ** 2) <= radius ** 2
            mask |= mask_circle

    # Save the mask to a FITS file
    header = fits.Header()
    header['SIMPLE'] = True
    header['BITPIX'] = 8
    header['NAXIS'] = 2
    header['NAXIS1'] = image_data.shape[1]
    header['NAXIS2'] = image_data.shape[0]

    # Assuming you have WCS information in your image_data, update the header accordingly
    wcs = WCS(header)

    # Save the mask FITS file
    mask_fits = fits.PrimaryHDU(mask.astype(np.uint8), header=header)
    mask_fits.writeto(output_fits_file, overwrite=True)

# Example usage
if __name__ == "__main__":
    # Specify the DS9 regions file and the output FITS file
    ds9_regions_file = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/AGEL_DR2/KCWI/Tucker_Keerthi/DCLS1507/counter_image2.reg'
    output_fits_file = '/Users/nandinisahu/Library/CloudStorage/OneDrive-UNSW/AGEL/AGEL_DR2/KCWI/Tucker_Keerthi/DCLS1507/counter_image2.fits'

    # Assuming you have an image already loaded
    # Replace this with your actual image data loading logic
    image_data = np.random.random((100,100))

    create_mask_fits(ds9_regions_file, output_fits_file, image_data)
