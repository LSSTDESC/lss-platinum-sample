#figure out how many galaxies are represented in each sample
#using the current converts
#then lyman alpha
#graenta et al
#Using the old Lyman alpha
#total rate in python

#Overall purpose:
#Using Lyman alpha data from the MUSYC survery in order to create a graph of Luminosity function and therefore SFR density
#In general this code aims to devlope a "Platinum sample":
#Galaxies that LSST can use to better measure the Baryon Acoustic Oscillation (BAO) feature, and thereby to better probe dark energy.
#Note that the broader the redshift bin, the more that projection washes out the BAO signal.
#Determine the number density of Extreme Emission Line galaxies (EELs) that will have significantly improved photometric redshifts from LSST due to their strong emission lines triggering one or more of the ugrizy photometric bands.
#We’ll be looking for very strong Lymanalpha, Halpha, [OIII], and [OII] emission lines.


import numpy
from matplotlib.pyplot import *
%matplotlib inline
import scipy
from scipy import integrate
from astropy.cosmology import Planck15 as cosmo
