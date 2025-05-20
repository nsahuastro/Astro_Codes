from astropy.io import fits
import numpy as np
import tkinter, tkinter.filedialog
import re
import math
import os
import subprocess
from shutil import copyfile
from scipy import ndimage
from functools import reduce
import astropy.io.fits as pyfits
import time
import random
root = tkinter.Tk()
root.withdraw()
import re
from astropy.cosmology import LambdaCDM

cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

#Input redshifts
z_lens=[]
z_lens.append(float(input("zlens: ")))
#z_lens=0.856
np.array(z_lens)

z_source=[]
z_source.append(float(input("zsource: ")))
#z_source=2.392
np.array(z_source)



#Calculate angular diameter distances
Dd=cosmo.angular_diameter_distance(z_lens)
Dd=float(Dd.value)
Ds=cosmo.angular_diameter_distance(z_source)
Ds=float(Ds.value)
Dds=cosmo.angular_diameter_distance_z1z2(z_lens,z_source)
Dds=float(Dds.value)
Dt=round((1+z_lens[0])*(Dd*Ds/Dds),8)

#Calculate time-delay distance
print("Dd: " + str(Dd) + '\n')
print("Ds: " + str(Ds) + '\n')
print("Dds: " + str(Dds) + '\n')
print("Dds/Ds: " + str(Dds/Ds) + '\n')
print("Dt: " + str(Dt) + '\n')
