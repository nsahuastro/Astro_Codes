"""
Created on Thu Apr 25 20:02:56 2024

@author: nandinisahu
"""
from astropy.coordinates import SkyCoord
from astropy import units as u

def convertDegreesToHMS(ra_deg:float ,dec_deg:float)->str:
    '''
    returns ra and dec in hms from degrees using astropy
    '''
    c = SkyCoord(
        ra    = ra_deg*u.degree,
        dec   = dec_deg*u.degree
    )
    return(c.to_string('hmsdms').replace('h',':').replace('d',':').replace('m',':').replace('s',''))


cc=convertDegreesToHMS(139.89600000, 30.53230000)

print(cc)
