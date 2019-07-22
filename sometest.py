import numpy
from numpy import loadtxt
#should reorganize this somehow
from matplotlib.pyplot import *
import scipy
from scipy import integrate
from astropy.cosmology import Planck15 as cosmo


log10Lstar = (m1*z) + b1
log10phistar = (m2*z) + b2

zlow = 1
zhigh = 3

log10Lstarlow = 39
log10Lstarhigh = 46

log10phistarlow = -7
log10phistarhigh = 0


m1 = (log10Lstarhigh - log10Lstarlow) / (zhigh - zlow)
m2 = (log10phistarhigh - log10phistarlow) / (zhigh - zlow)

print("m1 = ",m1)
print("m2 = ",m2)


b1test1 = log10Lstarlow - (m1*zlow)
b1test2 = log10Lstarhigh - (m1*zhigh)

print("b1test1 = ",b1test1)
print("b1test2 = ",b1test2)

b2test1 = log10phistarlow - (m2*zlow)
b2test2 = log10phistarhigh - (m2*zhigh)

print("b2test1 = ",b2test1)
print("b2test2 = ",b2test2)