import numpy
from matplotlib.pyplot import *
%matplotlib inline
import scipy
from scipy import integrate
from astropy.cosmology import Planck15 as cosmo


#Overall purpose:
#Using Lyman alpha data from the MUSYC survery in order to create a graph of Luminosity function and therefore SFR density
#Later: Finding the number of galaxies in each sample
#In general this code aims to devlope a "Platinum sample":
#Galaxies that LSST can use to better measure the Baryon Acoustic Oscillation (BAO) feature, and thereby to better probe dark energy.
#Note that the broader the redshift bin, the more that projection washes out the BAO signal.
#Determine the number density of Extreme Emission Line galaxies (EELs) that will have significantly improved photometric redshifts from LSST due to their strong emission lines triggering one or more of the ugrizy photometric bands.
#We’ll be looking for very strong Lymanalpha, Halpha, [OIII], and [OII] emission lines.
#The ugrizy filter options with corresponding final coadded depths (5 sigma) are as follows:
#u: 26.1, g: 27.4, r: 27.5, i: 26.8, z: 26.1, y: 24.9



filt = input("Filter?  (uband,gband,rband,iband,zband,yband)")
print("filt = ", filt)
