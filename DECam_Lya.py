print("This code uses the Schechter luminosity functions")
print("from linear parametrization from Ciardullo et al 2012 for Lymanalpha 1216 Angstroms")
print("to find the expected number of galaxies in possible new narrowband filters")
print("for the DECam on the CTIO telescope in Chile:")


import numpy
from numpy import loadtxt
#should reorganize this somehow
from matplotlib.pyplot import *
import scipy
from scipy import integrate
from astropy.cosmology import Planck15 as cosmo


#this code calculates some stuff first, then defines a few functions, and then has a bunch of if statements that run the functions 


#need to model DECam filter as a gaussian
#have it be 2000 points in each array

lambda_Lymanalpha = 121.6 #in nm
print("lambda_Lymanalpha",lambda_Lymanalpha)


# #testband
# #need the following to calculate the gaussian-like filter transmission
# z_array_testband = numpy.linspace #NO IDEA - NEED TO CHECK ?? Paper values or our values?  (2.59-0.08,2.59+0.08,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
# wav_array_testband = (z_array_testband+1)*lambda_Lymanalpha
# FWHM_val_testband = 50 #nm
# stand_dev_testband = FWHM_val_testband/(2*numpy.sqrt(2*numpy.log(2)))
# lambda_c_testband = 500.7 #nm
# mean_mu_testband = lambda_c_testband
# #calculating the gaussian-like filter transmission
# gaussian_dist_testband = (1/(stand_dev_testband*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_testband-mean_mu_testband)**2)/(2*(stand_dev_testband**2))))
# #the two values to plot
# DECam_filter_testband = gaussian_dist_testband #transmission
# DECam_wav_testband = wav_array_testband #wavelengths
# FWHMlow_wav_testband = wav_array_testband[0]+((wav_array_testband[-1]-wav_array_testband[0])/2)-(FWHM_val_testband/2) #nm
# FWHMhigh_wav_testband = wav_array_testband[-1]-((wav_array_testband[-1]-wav_array_testband[0])/2)+(FWHM_val_testband/2) #nm
# ABmag = 25.4

#there are:
#four possible gband narrowband filters
#three possible rband filters
#one possible zband filter

#extra filter
#rband4
#need the following to calculate the gaussian-like filter transmission
z_array_rband4 = numpy.linspace(4.54-0.04,4.54+0.04,num=2000)
wav_array_rband4 = (z_array_rband4+1)*lambda_Lymanalpha
FWHM_val_rband4 = 10 #nm 
stand_dev_rband4 = FWHM_val_rband4/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_rband4 = 673.5 #nm
mean_mu_rband4 = lambda_c_rband4
#calculating the gaussian-like filter transmission
gaussian_dist_rband4 = (1/(stand_dev_rband4*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_rband4-mean_mu_rband4)**2)/(2*(stand_dev_rband4**2))))
#the two values to plot
DECam_filter_rband4 = gaussian_dist_rband4 #transmission
DECam_wav_rband4 = wav_array_rband4 #wavelengths
FWHMlow_wav_rband4 = wav_array_rband4[0]+((wav_array_rband4[-1]-wav_array_rband4[0])/2)-(FWHM_val_rband4/2) #nm
FWHMhigh_wav_rband4 = wav_array_rband4[-1]-((wav_array_rband4[-1]-wav_array_rband4[0])/2)+(FWHM_val_rband4/2) #nm


#NEED TO FIX THE REST - PROPAGATE CORRECTED REDSHIFT RANGES INTO THE REST OF THE CODE AND SPREADSHEET

#four possible gband filters:

#gband1
#need the following to calculate the gaussian-like filter transmission
z_array_gband1 = numpy.linspace(2.59-0.04,2.59+0.04,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_gband1 = (z_array_gband1+1)*lambda_Lymanalpha
FWHM_val_gband1 = 10 #nm
stand_dev_gband1 = FWHM_val_gband1/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_gband1 = 436.3 #nm
mean_mu_gband1 = lambda_c_gband1
#calculating the gaussian-like filter transmission
gaussian_dist_gband1 = (1/(stand_dev_gband1*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_gband1-mean_mu_gband1)**2)/(2*(stand_dev_gband1**2))))
#the two values to plot
DECam_filter_gband1 = gaussian_dist_gband1 #transmission
DECam_wav_gband1 = wav_array_gband1 #wavelengths
FWHMlow_wav_gband1 = wav_array_gband1[0]+((wav_array_gband1[-1]-wav_array_gband1[0])/2)-(FWHM_val_gband1/2) #nm
FWHMhigh_wav_gband1 = wav_array_gband1[-1]-((wav_array_gband1[-1]-wav_array_gband1[0])/2)+(FWHM_val_gband1/2) #nm

#gband2
#need the following to calculate the gaussian-like filter transmission
z_array_gband2 = numpy.linspace(2.59-0.018,2.59+0.018,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_gband2 = (z_array_gband2+1)*lambda_Lymanalpha
FWHM_val_gband2 = 4.4 #nm
stand_dev_gband2 = FWHM_val_gband2/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_gband2 = 436.3 #nm
mean_mu_gband2 = lambda_c_gband2
#calculating the gaussian-like filter transmission
gaussian_dist_gband2 = (1/(stand_dev_gband2*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_gband2-mean_mu_gband2)**2)/(2*(stand_dev_gband2**2))))
#the two values to plot
DECam_filter_gband2 = gaussian_dist_gband2 #transmission
DECam_wav_gband2 = wav_array_gband2 #wavelengths
FWHMlow_wav_gband2 = wav_array_gband2[0]+((wav_array_gband2[-1]-wav_array_gband2[0])/2)-(FWHM_val_gband2/2) #nm
FWHMhigh_wav_gband2 = wav_array_gband2[-1]-((wav_array_gband2[-1]-wav_array_gband2[0])/2)+(FWHM_val_gband2/2) #nm

#gband3
#need the following to calculate the gaussian-like filter transmission
z_array_gband3 = numpy.linspace(3.12-0.04,3.12+0.04,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_gband3 = (z_array_gband3+1)*lambda_Lymanalpha
FWHM_val_gband3 = 10 #nm
stand_dev_gband3 = FWHM_val_gband3/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_gband3 = 500.7 #nm
mean_mu_gband3 = lambda_c_gband3
#calculating the gaussian-like filter transmission
gaussian_dist_gband3 = (1/(stand_dev_gband3*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_gband3-mean_mu_gband3)**2)/(2*(stand_dev_gband3**2))))
#the two values to plot
DECam_filter_gband3 = gaussian_dist_gband3 #transmission
DECam_wav_gband3 = wav_array_gband3 #wavelengths
FWHMlow_wav_gband3 = wav_array_gband3[0]+((wav_array_gband3[-1]-wav_array_gband3[0])/2)-(FWHM_val_gband3/2) #nm
FWHMhigh_wav_gband3 = wav_array_gband3[-1]-((wav_array_gband3[-1]-wav_array_gband3[0])/2)+(FWHM_val_gband3/2) #nm

#gband4
#need the following to calculate the gaussian-like filter transmission
z_array_gband4 = numpy.linspace(3.12-0.027,3.12+0.027,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_gband4 = (z_array_gband4+1)*lambda_Lymanalpha
FWHM_val_gband4 = 5 #nm
stand_dev_gband4 = FWHM_val_gband4/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_gband4 = 500.7 #nm
mean_mu_gband4 = lambda_c_gband4
#calculating the gaussian-like filter transmission
gaussian_dist_gband4 = (1/(stand_dev_gband4*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_gband4-mean_mu_gband4)**2)/(2*(stand_dev_gband4**2))))
#the two values to plot
DECam_filter_gband4 = gaussian_dist_gband4 #transmission
DECam_wav_gband4 = wav_array_gband4 #wavelengths
FWHMlow_wav_gband4 = wav_array_gband4[0]+((wav_array_gband4[-1]-wav_array_gband4[0])/2)-(FWHM_val_gband4/2) #nm
FWHMhigh_wav_gband4 = wav_array_gband4[-1]-((wav_array_gband4[-1]-wav_array_gband4[0])/2)+(FWHM_val_gband4/2) #nm


#three possible rband filters:

#rband1
#need the following to calculate the gaussian-like filter transmission
z_array_rband1 = numpy.linspace(4.40-0.04,4.40+0.04,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_rband1 = (z_array_rband1+1)*lambda_Lymanalpha
FWHM_val_rband1 = 10 #nm
stand_dev_rband1 = FWHM_val_rband1/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_rband1 = 656.3 #nm
mean_mu_rband1 = lambda_c_rband1
#calculating the gaussian-like filter transmission
gaussian_dist_rband1 = (1/(stand_dev_rband1*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_rband1-mean_mu_rband1)**2)/(2*(stand_dev_rband1**2))))
#the two values to plot
DECam_filter_rband1 = gaussian_dist_rband1 #transmission
DECam_wav_rband1 = wav_array_rband1 #wavelengths
FWHMlow_wav_rband1 = wav_array_rband1[0]+((wav_array_rband1[-1]-wav_array_rband1[0])/2)-(FWHM_val_rband1/2) #nm
FWHMhigh_wav_rband1 = wav_array_rband1[-1]-((wav_array_rband1[-1]-wav_array_rband1[0])/2)+(FWHM_val_rband1/2) #nm

#rband2
#need the following to calculate the gaussian-like filter transmission
z_array_rband2 = numpy.linspace(4.40-0.027,4.40+0.027,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_rband2 = (z_array_rband2+1)*lambda_Lymanalpha
FWHM_val_rband2 = 6.6 #nm
stand_dev_rband2 = FWHM_val_rband2/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_rband2 = 656.3 #nm
mean_mu_rband2 = lambda_c_rband2
#calculating the gaussian-like filter transmission
gaussian_dist_rband2 = (1/(stand_dev_rband2*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_rband2-mean_mu_rband2)**2)/(2*(stand_dev_rband2**2))))
#the two values to plot
DECam_filter_rband2 = gaussian_dist_rband2 #transmission
DECam_wav_rband2 = wav_array_rband2 #wavelengths
FWHMlow_wav_rband2 = wav_array_rband2[0]+((wav_array_rband2[-1]-wav_array_rband2[0])/2)-(FWHM_val_rband2/2) #nm
FWHMhigh_wav_rband2 = wav_array_rband2[-1]-((wav_array_rband2[-1]-wav_array_rband2[0])/2)+(FWHM_val_rband2/2) #nm

#rband3
#need the following to calculate the gaussian-like filter transmission
z_array_rband3 = numpy.linspace(4.45-0.07,4.45+0.07,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_rband3 = (z_array_rband3+1)*lambda_Lymanalpha
FWHM_val_rband3 = 17 #nm
stand_dev_rband3 = FWHM_val_rband3/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_rband3 = 662.2 #nm
mean_mu_rband3 = lambda_c_rband3
#calculating the gaussian-like filter transmission
gaussian_dist_rband3 = (1/(stand_dev_rband3*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_rband3-mean_mu_rband3)**2)/(2*(stand_dev_rband3**2))))
#the two values to plot
DECam_filter_rband3 = gaussian_dist_rband3 #transmission
DECam_wav_rband3 = wav_array_rband3 #wavelengths
FWHMlow_wav_rband3 = wav_array_rband3[0]+((wav_array_rband3[-1]-wav_array_rband3[0])/2)-(FWHM_val_rband3/2) #nm
FWHMhigh_wav_rband3 = wav_array_rband3[-1]-((wav_array_rband3[-1]-wav_array_rband3[0])/2)+(FWHM_val_rband3/2) #nm


#one possible zband:

#zband
#need the following to calculate the gaussian-like filter transmission
z_array_zband = numpy.linspace(6.93-0.04,6.93+0.04,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_zband = (z_array_zband+1)*lambda_Lymanalpha
FWHM_val_zband = 9.5 #nm
stand_dev_zband = FWHM_val_zband/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_zband = 964 #nm
mean_mu_zband = lambda_c_zband
#calculating the gaussian-like filter transmission
gaussian_dist_zband = (1/(stand_dev_zband*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_zband-mean_mu_zband)**2)/(2*(stand_dev_zband**2))))
#the two values to plot
DECam_filter_zband = gaussian_dist_zband #transmission
DECam_wav_zband = wav_array_zband #wavelengths
FWHMlow_wav_zband = wav_array_zband[0]+((wav_array_zband[-1]-wav_array_zband[0])/2)-(FWHM_val_zband/2) #nm
FWHMhigh_wav_zband = wav_array_zband[-1]-((wav_array_zband[-1]-wav_array_zband[0])/2)+(FWHM_val_zband/2) #nm


#zband_tenth
#need the following to calculate the gaussian-like filter transmission
z_array_zband_tenth = numpy.linspace(6.93-0.04,6.93+0.04,num=2000) #THIS IS IN REDSHIFT, IT NEEDS TO BE IN LAMBDA!!!!
wav_array_zband_tenth = (z_array_zband_tenth+1)*lambda_Lymanalpha
FWHM_val_zband_tenth = 9.5 #nm
stand_dev_zband_tenth = FWHM_val_zband_tenth/(2*numpy.sqrt(2*numpy.log(2)))
lambda_c_zband_tenth = 964 #nm
mean_mu_zband_tenth = lambda_c_zband_tenth
#calculating the gaussian-like filter transmission
gaussian_dist_zband_tenth = (1/(stand_dev_zband_tenth*numpy.sqrt(2*numpy.pi)))*(numpy.e**((-(wav_array_zband_tenth-mean_mu_zband_tenth)**2)/(2*(stand_dev_zband_tenth**2))))
#the two values to plot
DECam_filter_zband_tenth = gaussian_dist_zband_tenth #transmission
DECam_wav_zband_tenth = wav_array_zband_tenth #wavelengths
FWHMlow_wav_zband_tenth = wav_array_zband_tenth[0]+((wav_array_zband_tenth[-1]-wav_array_zband_tenth[0])/2)-(FWHM_val_zband_tenth/2) #nm
FWHMhigh_wav_zband_tenth = wav_array_zband_tenth[-1]-((wav_array_zband_tenth[-1]-wav_array_zband_tenth[0])/2)+(FWHM_val_zband_tenth/2) #nm


#now to plot the transmission functions of the narrowband filters
figure(4)

plot(DECam_wav_gband1,DECam_filter_gband1,"b",label="gband1")
plot(DECam_wav_gband2,DECam_filter_gband2,"b",label="gband2")
plot(DECam_wav_gband3,DECam_filter_gband3,"b",label="gband3")
plot(DECam_wav_gband4,DECam_filter_gband4,"b",label="gband4")

plot(DECam_wav_rband1,DECam_filter_rband1,"c",label="rband1")
plot(DECam_wav_rband2,DECam_filter_rband2,"c",label="rband2")
plot(DECam_wav_rband3,DECam_filter_rband3,"c",label="rband3")

plot(DECam_wav_zband,DECam_filter_zband,"y",label="zband")

xlim(400,1000)
ylim(0,0.25)
ylabel("Throughput")
xlabel("Wavelengths (nm)")
title("DECam Possible Narrowband Filter Transmission Functions")

legend(loc="upper right")

show()



print("(units are in Angstroms):")
print("Lymanalpha 1216")
emline = input("Plot Schechter Luminosity Function for which emission line?  (just type in: allends)")
print("emline = ", emline)



#this is the stuff that I need to do now: (actually updated)
#integrate the comovingphi array wrt the distance using astropy.cosmology (see end)


#FIX THESE NUMBERS
#the code will be able to use any filter
print("The possible narrowband filter options with corresponding final coadded depths (5 sigma) are as follows:")
print("gband1:  1 hour: 24.9  ;  10 hours: 26.2  ;  3 hours: 25.5  ;  6 hours: 25.9 ")
print("gband2:  1 hour: 24.5  ;  10 hours: 25.7  ;  3 hours: 25.1  ;  6 hours: 25.5 ")
print("gband3:  1 hour: 24.9  ;  10 hours: 26.2  ;  3 hours: 25.5  ;  6 hours: 25.9 ")
print("gband4:  1 hour: 24.6  ;  10 hours: 25.9  ;  3 hours: 25.2  ;  6 hours: 25.6 ")
print("rband1:  1 hour: 24.6  ;  10 hours: 25.8  ;  3 hours: 25.2  ;  6 hours: 25.6 ")
print("rband2:  1 hour: 24.3  ;  10 hours: 25.6  ;  3 hours: 24.9  ;  6 hours: 25.3 ")
print("rband3:  1 hour: 24.9  ;  10 hours: 26.1  ;  3 hours: 25.5  ;  6 hours: 23.9 ")
print("rband4:ALL OF THESE NEED TO BE UPDATED")
print("zband:   1 hour: 23.3  ;  10 hours: 24.5  ;  3 hours: 25.9  ;  6 hours: 24.3 ")



filt = input("Plot Schechter Luminosity Function for which filter?  (sample input: zband)")
print("filter = ", filt)

exp_time = input("Calculate for 1 hour or 10 hour exposure time?: (ex: 1hour or 10hour)")
print("exp_time = ",exp_time)

#THIS PART CONTAINS FUNCTIONS THAT WILL THEN BE USED LATER ON IN THE CODE


def lineqLya(z):

	#print("lineqLya has been called")

	#this is used for the Lymanalpha line

	#I AM NOT USING CIARDULLO ANYMORE, IF THIS WORKS
	# #finds linear equation from two points in Ciardullo et al. 2012 paper and interpolates/extraoplates to find Lstar and phistar for different z values

	# z1 = 3.113
	# z2 = 2.063
	# deltaz = z1 - z2

	# log10Lstar1 = 42.76
	# log10Lstar2 = 42.33

	# log10phistar1 = -3.17
	# log10phistar2 = -2.86

	#finds linear equation from two points in UPDATED Sobral et al. 2018 paper and interpolates/extraoplates to find Lstar and phistar for different z values

	#NOTE: NOT SURE IF I ASSUMED THE RIGHT THINGS, TAKEN FROM TABLE C3.

	z1 = 5.4
	z2 = 2.2
	deltaz = z1 - z2

	log10Lstar1 = 43.35
	log10Lstar2 = 42.69

	log10phistar1 = -3.78
	log10phistar2 = -3.33

	mlog10Lstar = (log10Lstar1 - log10Lstar2)/(deltaz)
	mlog10phistar = (log10phistar1 - log10phistar2)/(deltaz)

	blog10Lstar = log10Lstar1 - mlog10Lstar*z1
	blog10phistar = log10phistar1 - mlog10phistar*z1

	#checked to make sure values give same as original
	#Lstartest1 = mLstar*z1 + bLstar
	#Lstartest2 = mLstar*z2 + bLstar
	#phistartest1 = mphistar*z1 + bphistar
	#phistartest2 = mphistar*z2 + bphistar
	#print("Lstar1, Lstar 2 = ",Lstartest1,Lstartest2)
	#print("phistar1, phistar2 = ",phistartest1,phistartest2)

	#print("log10Lstar for Lymanalpha: log10Lstar = z*",mlog10Lstar," + ",blog10Lstar)
	#print("log10phistar for Lymanalpha: log10phistar = z*",mlog10phistar," + ",blog10phistar)

	log10Lstar = mlog10Lstar*z + blog10Lstar
	log10phistar = mlog10phistar*z + blog10phistar

	LstarLya = 10**log10Lstar
	phistarLya = 10**log10phistar

	#print("LstarLya =",LstarLya)
	#print("phistarLya =",phistarLya)

	answers = [LstarLya,phistarLya]

	return answers


def filter_int(filt):

	#print("filter_int has been called")

	if filt=="gband1":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_gband1 #transmission or throughput
		DECamfilter = wav_array_gband1 #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_gband1
		lambdaFWHMright = FWHMhigh_wav_gband1

	if filt=="gband2":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_gband2 #transmission or throughput
		DECamfilter = wav_array_gband2 #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_gband2
		lambdaFWHMright = FWHMhigh_wav_gband2

	if filt=="gband3":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_gband3 #transmission or throughput
		DECamfilter = wav_array_gband3 #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_gband3
		lambdaFWHMright = FWHMhigh_wav_gband3

	if filt=="gband4":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_gband4 #transmission or throughput
		DECamfilter = wav_array_gband4 #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_gband4
		lambdaFWHMright = FWHMhigh_wav_gband4

	if filt=="rband1":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_rband1 #transmission or throughput
		DECamfilter = wav_array_rband1 #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_rband1
		lambdaFWHMright = FWHMhigh_wav_rband1

	if filt=="rband2":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_rband2 #transmission or throughput
		DECamfilter = wav_array_rband2 #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_rband2
		lambdaFWHMright = FWHMhigh_wav_rband2

	if filt=="rband3":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_rband3 #transmission or throughput
		DECamfilter = wav_array_rband3 #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_rband3
		lambdaFWHMright = FWHMhigh_wav_rband3

	if filt=="zband":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_zband #transmission or throughput
		DECamfilter = wav_array_zband #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_zband
		lambdaFWHMright = FWHMhigh_wav_zband

	if filt=="zband_tenth":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_zband_tenth #transmission or throughput
		DECamfilter = wav_array_zband_tenth #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_zband_tenth
		lambdaFWHMright = FWHMhigh_wav_zband_tenth

	if filt=="rband4":
		#use filter transmission functions calculated at beginning of code
		DECam_filter = gaussian_dist_rband4 #transmission or throughput
		DECamfilter = wav_array_rband4 #wavelengths
		DECamlambda = DECamfilter
		#want to find the FWHM of each filter
		#easier to just take directly from gaussian / given in beginning
		lambdaFWHMleft = FWHMlow_wav_rband4
		lambdaFWHMright = FWHMhigh_wav_rband4


	#to find the midpoint, integrate over the entire filter, then find what value would give you half the integrated value
	#see if you need to fix this because should have x=
	stepint_DECamfilter = scipy.integrate.cumtrapz(DECamfilter)
	#print('stepint_DECamfilter')
	#print(stepint_DECamfilter)
	len(stepint_DECamfilter)
	#len gives total length, which is index number + 1
	lastnumber = stepint_DECamfilter[len(stepint_DECamfilter)-1]
	midpoint = lastnumber/2

	#finds closest value by finding minimum difference
	difference = abs(midpoint-stepint_DECamfilter)
	#print("difference = ",difference)
	mindiff = min(difference)
	#print("mindiff = ",mindiff)
	midindex = numpy.where(difference==mindiff)
	#print("midindex = ", midindex)
	#print("type is ",type(midindex))
	#have to change this from tuple to float - first find the actual value in index 0
	midindex = midindex[0]
	midindex = numpy.float(midindex)
	#print("midindex = ", midindex)

	#now I have to figure out if it is really the midpoint of the integrated flux(??)
	#Zfilterleft = zfilter[:658]
	#zfilterright = zfilter[658:]
	#a = scipy.integrate.trapz(zfilterleft)
	#b = scipy.integrate.trapz(zfilterright)
	#print("leftint = ",a)
	#print("rightint = ",b)
	#after testing 656,657,658, I concluded that 658 is the middle - I have to wrap my head around the indices in python and how trapz affects them

	#since 658 = midindex+1
	centerindex = midindex+1
	centerindex = int(centerindex)
	#print("centerindex = ",centerindex)
	#print(type(centerindex))
	centerval = DECam_filter[centerindex]
	#print("centerval = ",centerval)
	#print(type(centerval))
	lambdacenter = centerval#[0]
	#print("lambdacenter = ",lambdacenter)
	lambdacenter = numpy.float(lambdacenter)
	#print("median transmission wavelength of the "+filt+" = ",lambdacenter)

	#the "central" wavelength is called the median transmission wavelength of each band


	#should have output as an array (see LF code)
	#print("returns array: lambdacenter,lambdaFWHMleft,lambdaFWHMright")
	lambdaarray = [lambdacenter,lambdaFWHMleft,lambdaFWHMright]

	return lambdaarray


#the following function calculates the luminosity limit for a certain band with a certain detection limit
def lumlim(z,em,filt,exp_time):

	#print("lumlim has been called")

	#REMINDER TO SELF - z is a function of the emission line wavelength and the filter itself, so it does not have to be adjusted for the FWHM thing

	#print("This function will calculate the luminosity that corresponds to a 5 sigma detection, in erg/s:")
	#print("INFO: 26.2 is AB magnitude for 5 sigma detection limit in the z band") #this was used when I only had the z band
	#first I calculate the flux, then convert to flux density, then find the luminosity limit for the conditions printed above

	#changed if statements (now deleted) to dictionaries below
	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	if exp_time=="1hour":
		ABmag_dict = {"gband1":24.9,"gband2":24.5,"gband3":24.9,"gband4":24.6,"rband1":24.6,"rband2":24.3,"rband3":24.9,"zband":23.3,"zband_tenth":23.3,"rband4":24.6}

	if exp_time=="3hour":
		ABmag_dict = {"gband1":25.5,"gband2":25.1,"gband3":25.5,"gband4":25.2,"rband1":25.2,"rband2":24.9,"rband3":25.5,"zband":23.9,"zband_tenth":0,"rband4":25.2}
	
	if exp_time=="6hour":
		ABmag_dict = {"gband1":25.9,"gband2":25.5,"gband3":25.9,"gband4":25.6,"rband1":25.6,"rband2":25.3,"rband3":25.9,"zband":24.3,"zband_tenth":0,"rband4":25.6}
	
	if exp_time=="7.5hour":
		ABmag_dict = {"rband4":25.7}

	if exp_time=="10hour":
		ABmag_dict = {"gband1":26.2,"gband2":25.7,"gband3":26.2,"gband4":25.9,"rband1":25.8,"rband2":25.6,"rband3":26.1,"zband":24.5,"zband_tenth":24.5,"rband4":25.8}

	if exp_time=="15hour":
		ABmag_dict = {"rband4":26}


	lambdalow_dict = {"gband1":DECam_wav_gband1[0],"gband2":DECam_wav_gband2[0],"gband3":DECam_wav_gband3[0],"gband4":DECam_wav_gband4[0],"rband1":DECam_wav_rband1[0],"rband2":DECam_wav_rband2[0],"rband3":DECam_wav_rband3[0],"zband":DECam_wav_zband[0],"zband_tenth":DECam_wav_zband[0],"rband4":DECam_wav_rband4[0]} #in nm"
	lambdahigh_dict = {"gband1":DECam_wav_gband1[-1],"gband2":DECam_wav_gband2[-1],"gband3":DECam_wav_gband3[-1],"gband4":DECam_wav_gband4[-1],"rband1":DECam_wav_rband1[-1],"rband2":DECam_wav_rband2[-1],"rband3":DECam_wav_rband3[-1],"zband":DECam_wav_zband[-1],"zband_tenth":DECam_wav_zband[-1],"rband4":DECam_wav_rband4[-1]} #in nm
	#these are at the endpoints of the filters, NOT the FWHM

	ABmag = ABmag_dict[filt]
	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]

	#might need the following here:

	if filt=="gband1":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_gband1 #transmission or throughput
		DECamwav = wav_array_gband1 #wavelengths

	if filt=="gband2":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_gband2 #transmission or throughput
		DECamwav = wav_array_gband2 #wavelengths

	if filt=="gband3":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_gband3 #transmission or throughput
		DECamwav = wav_array_gband3 #wavelengths

	if filt=="gband4":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_gband4 #transmission or throughput
		DECamwav = wav_array_gband4 #wavelengths

	if filt=="rband1":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_rband1 #transmission or throughput
		DECamwav = wav_array_rband1 #wavelengths

	if filt=="rband2":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_rband2 #transmission or throughput
		DECamwav = wav_array_rband2 #wavelengths

	if filt=="rband3":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_rband3 #transmission or throughput
		DECamwav = wav_array_rband3 #wavelengths

	if filt=="zband":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_zband #transmission or throughput
		DECamwav = wav_array_zband #wavelengths

	if filt=="zband_tenth":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_zband_tenth #transmission or throughput
		DECamwav = wav_array_zband_tenth

	if filt=="rband4":
		#use filter transmission functions calculated at beginning of code
		DECamfilter = gaussian_dist_rband4 #transmission or throughput
		DECamwav = wav_array_rband4 #wavelengths


	#the following calculates values I use in and plug into the lumlim function

	#need to redo how I find the flux density to take into account the shape of the transmission curve of each filter
	#I will get a (limiting) flux for this (see notes)


	#DECamfilter is the transmission array

	#steps to take: FOR EACH FILTER


	#(1)
	#n_photon(1_microJansky) = int{Transmission(lambda)*[1microJansky/(hc/lambda)]*dlambda}
	#hc/lambda is the energy per photon
	#lambdaem is the emitted wavelength
	#do this for the chosen filter (corresponding to the transmission curve)
	#check units - make sure n ends up as number/((cm^2)*s)

	h_cgs = 6.6261*(10**(-27)) #g*cm^2/s
	c_cgs = 2.9979*(10**10) #in cm/s
	lambdaem_dict_cgs = {"Lymanalpha":121.6/(10**7)} #in cm
	#print("TEST HERE: lambdaem_dict_cgs =",lambdaem_dict_cgs)
	lambdaem_cgs = lambdaem_dict_cgs[em]
	#print("TEST HERE: lambdaem_cgs =",lambdaem_cgs)
	#print("TEST HERE: DECamwav =",DECamwav)
	DECamwav = DECamwav/(10**7)  #nm to cm
	#print("TEST HERE: DECamwav =",DECamwav)

	#print("TEST HERE: DECamfilter: (this may be the issue?",DECamfilter)
	n_photon_1microJansky_array = DECamfilter*(10**(-29))/(h_cgs*DECamwav) #c_cgs/ #fixed the lambda from an incorrect value to the DECamwav from the transmission function array
	#print("TEST HERE: n_photon_1microJansky_array",n_photon_1microJansky_array)
	#dlambda = abs(DECamwav[1]-DECamwav[0])
	#print("TEST HERE: dlambda",dlambda)
	#n_photon_1microJansky = numpy.trapz(n_photon_1microJansky_array,dx=dlambda)
	n_photon_1microJansky = numpy.trapz(n_photon_1microJansky_array,x=DECamwav) #more general
	#print("n_photon_1microJansky",n_photon_1microJansky)
	#print("TEST HERE: n_photon_1microJansky = ",n_photon_1microJansky)


	#(2)
	#lambda_emissionline = lambda_restframe*(1+z_emissionline)

	lambda_emissionline = lambdaem_cgs*(1+z)  #in cm
	#print("TEST HERE: z =",z)
	#print("TEST HERE: lambdaem_cgs =",lambdaem_cgs)
	#print("TEST HERE: lambda_emissionline =",lambda_emissionline)


	#(3)

	#need DECamfilter[lambda_emissionline]
	#will find this using minimum difference between lambda_emissionline and whatever wavelengths are in the array

	#DECamfilter is the transmission
	#DECamwav is the wavelength array
	#ultimately need a value in DECamfilter using an index from DECamwav

	#finds the index that most closely matches the wavelength that is redshifted

	diff_lambda_array = abs(DECamwav-lambda_emissionline) #SHOULD BE IN SAME UNITS - they are both in centimeters in this section of the code
	#do without absolute value
	#print("TEST HERE: DECamwav =",DECamwav)
	#print("TEST HERE: diff_lambda_array =",diff_lambda_array)

	#answer of following should be ~0
	mindiff_lambda = min(diff_lambda_array)
	#print("TEST HERE: mindiff_lambda =", mindiff_lambda)
	mindiff_lambda = numpy.float(mindiff_lambda)
	#print("TEST HERE: mindiff_lambda =",mindiff_lambda)

	#finds index corresponding to value closest to wavelength
	mindiff_lambda_index = numpy.where(diff_lambda_array==mindiff_lambda)

	#find "value" in index 0, then change from tuple to float
	#print("TEST HERE: mindiff_lambda_index =",mindiff_lambda_index)
	mindiff_lambda_index = mindiff_lambda_index[0]
	#print("TEST HERE: mindiff_lambda_index =",mindiff_lambda_index)
	mindiff_lambda_index = numpy.float(mindiff_lambda_index)
	#print("TEST HERE: mindiff_lambda_index =",mindiff_lambda_index)


	#this is the value of the transmission that corresponds to the index found above
	#DECamfilter[mindiff_lambda_index]
	#print("TEST HERE: DECamfilter[mindiff_lambda_index] =",DECamfilter[mindiff_lambda_index])


	#HERE STOP TEST


	#(4)
	#flux_limiting(emissionline) = [F_nu,lim(filt)*n_photon(1_microJansky)/T_filt(lambda_emissionline)]*(hc/lambda_emissionline)
	#F_nu,lim(filt) is the variable I had previously named fluxdens - it is the 5 sigma AB magnitude limit that I found for each DECam filter

	#finds the flux density
	#print("ABmagnitude = -2.5*log10(fluxdensity/(3631 Jansky))")
	#print("consequently:")
	#print("fluxdensity = (10**(ABmagnitude/(-2.5)))*(3631 Janksy)")
	#fluxdens = (10**(ABmag/(-2.5)))*3631 #outputs in Jansky
	#fluxdens = (10**(ABmag/(-2.5)))*3631*(10**(-6)) #outputs in microJansky
	#uses ABmag from earlier in this function
	fluxdens = 10**((ABmag+48.6)/(-2.5)) #outputs in erg/(s*Hz*(cm^2))
	#print("fluxdens",fluxdens)
	#print("flux density =",fluxdens,"erg/(s*Hz*(cm^2))")

	#DECamfilter is the transmission
	#print("f_n_obj / [1_microJansky] = n_photon_object / n_photon_1_microJansky")
	flux_limit = (fluxdens*n_photon_1microJansky/(DECamfilter[int(mindiff_lambda_index)]*(10**(-29))))*(h_cgs*c_cgs/lambda_emissionline)
	#print("flux_limit",flux_limit)

	#print("FLUX LIMIT IS: for z=",z)
	#print(" and em=",em)
	#print(" is:",flux_limit)


	#print("ADAM'S TEST:")
	#print("hc/lambda should be ~1.6e-11 ergs:",h_cgs*c_cgs/lambdaem_cgs)
	#print("T_EL should be of order 0.2:",DECamfilter[mindiff_lambda_index])
	#print("since 1 microJansky is 10^-29, limiting flux density seems good at 1.3*10^-30:",fluxdens)
	#print("n_photon_1microJansky is too small?:",n_photon_1microJansky)

	#THE FOLLOWING IS GOOD:

	#finds the luminosity distance
	lumdist = cosmo.luminosity_distance(z)
	#this outputs a special object that keeps track of units, so first I convert it to cm (cgs units), and then I convert it to a regular number
	lumdist_cgs = lumdist.to('cm')
	lumdist_unitless = lumdist_cgs.value
	#print("the luminosity distance for redshift z =",z,"is lumdist =",lumdist_unitless,"cm")

	#finds the luminosity limit
	#print("Luminosity = 4*pi*(luminositydistance**2)*flux")
	lumlimit = 4*numpy.pi*(lumdist_unitless**2)*flux_limit
	#print("luminosity limit for 5 sigma detection of",em,"in "+filt+" band is",lumlimit,"ergs/s")

	#print("LUMINOSITY LIMIT IS: for z=",z)
	#print(" and em=",em)
	#print(" is:",lumlimit)

	#using return makes the main output of this function the value of lumlimit so that I can use it to calculate other things when I call this function
	return lumlimit


#the following defines and plots the Schechter luminosity function for a chosen emission line along with a separate plot of the number density
#it also calculates the number density above a luminosity limit calculated for a specific band in the previous function
#to use this function, type schechter_LF(z=redshifttoplot) or with additional parameters you want to change inside the ()
def schechter_LF(z,lambdaemitted,alpha,Lstar0,betaL,phistar0,betaphi,zpaper,param,fluxscale,em,filt,style = ""):

	#print("schechter_LF has been called")

	#REMINDER TO SELF - z is a function of the emission line wavelength and the filter itself, so it does not have to be adjusted for the FWHM thing 

	#print("alpha = ", alpha)
	#print("z = ", z)
	#print("emission line = ",em)

	log10Lstararray = numpy.arange(30,55,0.01)
	L = (10**log10Lstararray)

	#I have two different parametrizations for Lstar and phistar; the first one is from Comparat et al 2016, and the second one is from Sobral et al 2015
	#I added a third one when I linearly parametrized Lymanalpha from Ciardullo et al 2012

	# if param == "first":

	# 	Lstar = fluxscale*Lstar0*((1+z)**betaL)
	# 	#print("Lstar = ", Lstar)

	# 	phistar = phistar0*((1+z)**betaphi)
	# 	#print("phistar = ", phistar)


	# if param == "second":

	# 	Lstar = 10**(0.45*z + Lstar0)
	# 	#print("Lstar = ", Lstar)

	# 	phistar = 10**(-0.38*(z**2) + z + phistar0)
	# 	#print("phistar = ", phistar)


	if param == "third":

		answersLya = lineqLya(z=z)
		#print(type(answersLya))

		Lstar = answersLya[0] #this is linearly parametrized from Ciardullo+ 2012
		#print("Lstar = ",Lstar)
		#print(type(Lstar))

		phistar = answersLya[1] #this is linearly parametrized from Ciardullo+ 2012
		#print("phistar = ",phistar)
		#print(type(phistar))

	# if param == "fourth" #this is a cross-check with Gronwall et al. 2007, which is referenced in Gawiser et al. 2007

	if param == "fourth":

		print("Gronwall+2007 cross-check, only valid at z~3.1")

		Lstar = 10**(42.66)

		phistar = 1.28*(10**(-3))
		

	if param == "fifth":

		print("Updated: Sobral+ 2018")

		answersLya = lineqLya(z=z)
		#print(type(answersLya))

		Lstar = answersLya[0]

		phistar = answersLya[1]

	# 	#DO THIS LATER


	#testing - weird results with smaller bands, need to see if follows LF and transmission function
	#decrease Lstar and phistar with redshift - the different wavelength sections of the y-band filter correspond to redshifts for each emission line (THINK!!)

#	if param == "first" or "second" or "third":

#		Lstar = 10**((-0.75)*z + 49.5)
#		print("Lstar = ",Lstar)

#		phistar = 10**((-0.75)*z + 3.5)
#		print("phistar = ",phistar)


	#now need an array of phi for the given z
	phi = phistar*((L/Lstar)**(alpha+1))*(numpy.e**(-L/Lstar))
	#print("HERE IS THE ULTIMATE TEST KINDA:")
	#print("phistar",phistar)
	#print("L",L)
	#print("Lstar",Lstar)
	#print("alpha",alpha)
	#print("phi",phi)

	#this deletes parts of the arrays that are so small python counts them as zero; otherwise, I would not be able to take the logarithm of the array
	L = L[numpy.where(phi!=0)]
	log10Lstararray_nonzero = log10Lstararray[numpy.where(phi!=0)]
	phi = phi[numpy.where(phi!=0)]
	#code suddenly does not know when to stop?  add break:
	#if len(phi)==0:
	#	break
	#this should work?  

	log10L = numpy.log10(L)
	log10phi = numpy.log10(phi)

	#the following calculates values I use in and plug into the lumlim function

	#once again, dictionaries are better-formated than the previous (deleted) if statements
	lambdalow_dict = {"gband1":DECam_wav_gband1[0],"gband2":DECam_wav_gband2[0],"gband3":DECam_wav_gband3[0],"gband4":DECam_wav_gband4[0],"rband1":DECam_wav_rband1[0],"rband2":DECam_wav_rband2[0],"rband3":DECam_wav_rband3[0],"zband":DECam_wav_zband[0],"zband_tenth":DECam_wav_zband[0],"rband4":DECam_wav_rband4[0]} #in nm"
	lambdahigh_dict = {"gband1":DECam_wav_gband1[-1],"gband2":DECam_wav_gband2[-1],"gband3":DECam_wav_gband3[-1],"gband4":DECam_wav_gband4[-1],"rband1":DECam_wav_rband1[-1],"rband2":DECam_wav_rband2[-1],"rband3":DECam_wav_rband3[-1],"zband":DECam_wav_zband[-1],"zband_tenth":DECam_wav_zband[-1],"rband4":DECam_wav_rband4[-1]} #in nm
	#these are at the endpoints of the filters, NOT the FWHM

	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	FWHMlow = lambdaarray[1]
	FWHMhigh = lambdaarray[2]

	#print("wavelengths in nm")
	#print("THIS IS A TEST, lambdacenter",lambdacenter)
	#print("THIS IS A TEST, lambdalow",lambdalow)
	#print("THIS IS A TEST, lambdahigh",lambdahigh)
	#print("THIS IS A TEST, FWHMlow",FWHMlow)
	#print("THIS IS A TEST, FWHMhigh",FWHMhigh)

	#uses previously defined function to get median transmission wavelenth
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	#lambdaemitted = #will be one of the three emission lines #uses nm #these are redshifts
	#filterendlow = (lambdalow/lambdaemitted)-1  #0
	#filterendhigh = (lambdahigh/lambdaemitted)-1  #0
	filtercenter = (lambdacenter/lambdaemitted)-1
	FWHMendlow = (FWHMlow/lambdaemitted)-1
	FWHMendhigh = (FWHMhigh/lambdaemitted)-1


	#beginning of edits

	#print("REDSHIFT",z)

	#luminosity limit OF THE REDSHIFTED EMISSION LINE
	lumlim_emline = lumlim(z = z,em = em,filt = filt,exp_time=exp_time)


	#print("luminosity limit for the redshifted emission line:",lumlim_emline)


	#beginning of edits and temporary commenting out of plot stuff

	#finds luminosity limit in center of each band
	center = lumlim(z = filtercenter,em = em,filt = filt,exp_time=exp_time)

	#FWHMlow_lumlim = lumlim(z = ((908.5/lambdaemitted)-1),em = em,filt = filt)
	#print("THE LUMINOSITY LIMIT FOR "+em+"(at lambda = 908.5 nm) = ",FWHMlow_lumlim)


	#this was the test for the uband

	#305.3, 331.125, 356.95, 382.775, 408.6
	#318.2125
	#344.0375
	#369.8625
	#395.6875


	#quarter1 = lumlim(z = ((318.2125/lambdaemitted)-1),em = em,filt = filt)
	#print("THE LUMINOSITY LIMIT FOR "+em+"(at lambda = 318.2125 nm) = ",quarter1)

	#quarter2 = lumlim(z = ((344.0375/lambdaemitted)-1),em = em,filt = filt)
	#print("THE LUMINOSITY LIMIT FOR "+em+"(at lambda = 344.0375 nm) = ",quarter2)

	#quarter3 = lumlim(z = ((369.8625/lambdaemitted)-1),em = em,filt = filt)
	#print("THE LUMINOSITY LIMIT FOR "+em+"(at lambda = 369.8625 nm) = ",quarter3)

	#quarter4 = lumlim(z = ((395.6875/lambdaemitted)-1),em = em,filt = filt)
	#print("THE LUMINOSITY LIMIT FOR "+em+"(at lambda = 395.6875 nm) = ",quarter4)


	#finds luminosity limits at endpoints of FWHM
	FWHMlowpoint = lumlim(z = FWHMendlow,em = em,filt = filt,exp_time=exp_time)
	FWHMhighpoint = lumlim(z = FWHMendhigh,em = em,filt = filt,exp_time=exp_time)

	#print("luminosity limits")
	#print("OVER HEEEEEEEERE: center",center)
	#print("OVER HEEEEEEEERE: FWHMlowpoint",FWHMlowpoint)
	#print("OVER HEEEEEEEERE: FWHMhighpoint",FWHMhighpoint)

	#sets up an array that I use to plot
	lumarray = numpy.full(15,numpy.log10(center)) #this is an array so that I can plot a line - can use for all lumlim dashed lines
	yarray = numpy.arange(-10,5,1) #used to plot a line - can use for all lumlim dashed lines

	#should replicate this for the FWHM endpoints
	lumarraylow = numpy.full(15,numpy.log10(FWHMlowpoint))
	lumarrayhigh = numpy.full(15,numpy.log10(FWHMhighpoint))

	#lumarray_1 = numpy.full(15,numpy.log10(quarter1))
	#lumarray_2 = numpy.full(15,numpy.log10(quarter2))
	#lumarray_3 = numpy.full(15,numpy.log10(quarter3))
	#lumarray_4 = numpy.full(15,numpy.log10(quarter4))

	#alpha parameter sets the transparency/opacity (from 0 to 1)
	figure(1)
	plot(log10L,log10phi,style,alpha=0.10)#,label = zpaper)
	plot(lumarray,yarray,style+"--")
	plot(lumarraylow,yarray,style+"--")
	plot(lumarrayhigh,yarray,style+"--")
	#plot(lumarray_1,yarray,style+"--")
	#plot(lumarray_2,yarray,style+"--")
	#plot(lumarray_3,yarray,style+"--")
	#plot(lumarray_4,yarray,style+"--")
	xlim(35,46)
	ylim(-7,0)
	legend(loc = "upper right")
	#LaTeX is used with $ signs in the titles below
	xlabel("$\log_{10}(L [ergs/s])$")
	ylabel("$\log_{10}(\phi [\t{Mpc}^{-3}])$")
	title("Schechter Luminosity Function, "+filt)
	#print("z = ", z, " and alpha = ", alpha, " are plotted, ")
	#saves image
	#pyplot.savefig('/home/lanaeid/Desktop/fig1'+filter+'.png',bbox_inches = 'tight')

	show()

	#print("now the number density is calculated by integrating the LF using cumtrapz:")
	#cumptrapz integrates the opposite way than I need to integrate, so I flip it twice in the process
	phiflip = phi[::-1]
	#print(phi)
	#print(len(phiflip),len(log10Lstararray_nonzero))
	#always pay attention to the version of the LF that is being used
	#NOTE TO FIX? : phiflip is backwards, but the L array is not
	phiflipint = scipy.integrate.cumtrapz(phiflip,x=log10Lstararray_nonzero)
	num_dens = phiflipint[::-1]

	Lmin = numpy.delete(L,(len(L)-1))
	#shifted_Lmin = Lmin*fluxscale
	log10Lmin = numpy.log10(Lmin)
	#log10shifted_Lmin = numpy.log10(shifted_Lmin)

	#alpha parameter sets the transparency/opacity (from 0 to 1)
	figure(2)
	plot(log10Lmin,num_dens,style,alpha=0.10)#,label = zpaper)
	plot(lumarray,yarray,style+"--")
	plot(lumarraylow,yarray,style+"--")
	plot(lumarrayhigh,yarray,style+"--")
	#plot(lumarray_1,yarray,style+"--")
	#plot(lumarray_2,yarray,style+"--")
	#plot(lumarray_3,yarray,style+"--")
	#plot(lumarray_4,yarray,style+"--")
	legend(loc = "upper right")
	xlim(35,46)
	ylim(0,0.015)
	#LaTeX is used with $ signs in the titles below
	xlabel("$\log_{10}(L_{min} [ergs/s])$")
	ylabel("$\log_{10}(\phi [\t{Mpc}^{-3}])$")
	title("Number Density, "+filt)
	#saves image
	#pyplot.savefig('/home/lanaeid/Desktop/fig2'+filter+'.png',bbox_inches = 'tight')

	show()

	#end of plot code that I temporarily commented out


	#the following finds the comoving number density above a certain detection limit (shot noise limit),
	#which I calculated in the previous part of the code (also shown below), where I saved the variable named center

	#runs through the LF code for chosen emission line
	#then gets number density above the luminosity limit using the center value


	#THIS IS ALL I NEED: I edited so center is what it should be lumlim_emline

	#first, I have to shorten the phi array to contain only values above the luminosity limit
	#shortens the array to be above the luminosity limit, then integrates to get comoving number density
	philim = phi[numpy.where(L>lumlim_emline)]
	#log10Lstararray = numpy.arange(30,55,0.01)
	#L = (10**log10Lstararray)
	log10Lstararray_lumlim = log10Lstararray_nonzero[numpy.where(L>lumlim_emline)]
	comovingphi = scipy.integrate.trapz(philim,x=log10Lstararray_lumlim)

	#print("comovingphi = ",comovingphi)
	#print("comovingphi =",comovingphi,"Mpc^-3")


#	#don't need for now here, kept just for reference
#	arealphi = totalnumgalaxies/(4*numpy.pi)
#	print("arealphi =",arealphi,"steradian^-1")


	show()

	#print("comovingphi is returned for the input")
	
	return comovingphi






if emline == "allends":



	#changed if statements to dictionaries and deleted them
	lambdalow_dict = {"gband1":DECam_wav_gband1[0],"gband2":DECam_wav_gband2[0],"gband3":DECam_wav_gband3[0],"gband4":DECam_wav_gband4[0],"rband1":DECam_wav_rband1[0],"rband2":DECam_wav_rband2[0],"rband3":DECam_wav_rband3[0],"zband":DECam_wav_zband[0],"zband_tenth":DECam_wav_zband[0],"rband4":DECam_wav_rband4[0]} #in nm"
	lambdahigh_dict = {"gband1":DECam_wav_gband1[-1],"gband2":DECam_wav_gband2[-1],"gband3":DECam_wav_gband3[-1],"gband4":DECam_wav_gband4[-1],"rband1":DECam_wav_rband1[-1],"rband2":DECam_wav_rband2[-1],"rband3":DECam_wav_rband3[-1],"zband":DECam_wav_zband[-1],"zband_tenth":DECam_wav_zband[-1],"rband4":DECam_wav_rband4[-1]} #in nm
	#these are at the endpoints of the filters, NOT the FWHM

	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	FWHMlow = lambdaarray[1]
	FWHMhigh = lambdaarray[2]
	

	lambda_Lymanalpha = 121.6 #in nm

	#finds redshift of FWHM
	zFWHMlow = (FWHMlow/lambda_Lymanalpha)-1
	zFWHMhigh = (FWHMhigh/lambda_Lymanalpha)-1


	#first define array of 100 evenly spaced wavelengths
	if filt=="gband1":
		z_Lymanalpha_array = numpy.linspace(2.59-0.08,2.59+0.08,num=100)

	if filt=="gband2":
		z_Lymanalpha_array = numpy.linspace(2.59-0.036,2.59+0.036,num=100)

	if filt=="gband3":
		z_Lymanalpha_array = numpy.linspace(3.12-0.08,3.12+0.08,num=100)

	if filt=="gband4":
		z_Lymanalpha_array = numpy.linspace(3.12-0.045,3.12+0.045,num=100)

	if filt=="rband1":
		z_Lymanalpha_array = numpy.linspace(4.40-0.08,4.40+0.08,num=100)

	if filt=="rband2":
		z_Lymanalpha_array = numpy.linspace(4.40-0.054,4.40+0.054,num=100)

	if filt=="rband3":
		z_Lymanalpha_array = numpy.linspace(4.45-0.14,4.45+0.14,num=100)

	if filt=="zband":
		z_Lymanalpha_array = numpy.linspace(6.93-0.08,6.93+0.08,num=100)
		# print("z_Lymanalpha_array")

	if filt=="zband_tenth":
		z_Lymanalpha_array = numpy.linspace(6.93-0.08,6.93+0.08,num=100)
		# print("z_Lymanalpha_array")

	if filt=="rband4":
		z_Lymanalpha_array = numpy.linspace(4.54-0.04,4.54+0.04,num=100)

	#don't want negative redshifts
	z_Lymanalpha_array = z_Lymanalpha_array[numpy.where(z_Lymanalpha_array>0)]
	#print("z_Lymanalpha_array",z_Lymanalpha_array)

 	#NEED TO HAVE AN ARRAY OF Z VALUES TO INTEGRATE PHI
 	#set up empty arrays to print out and eventually integrate comovingphi and arealphi
	comovingphiarrayLymanalpha = numpy.zeros(100)


	print("start Lymanalpha")


	if len(z_Lymanalpha_array)>0:

		print("NOTE: FOR THE DECAM VALUES, USE ALPHA = -1.8, AND FOR GRONWALL+ 2007, USE ALPHA = -1.36")

		for l in range(len(z_Lymanalpha_array)):
			zLymanalpha = z_Lymanalpha_array[l]
			#the alpha parameter SHOULD BE updated to the Sobral+ 2018 value, which is valid from z~2.2-5.4, with alpha=-1.93, I just need to figure out how to parametrize the values given in the paper, as I am having trouble finding any set parametrization of Lstar and phistar 
			comovingphiarrayLymanalpha[l] = schechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.8,Lstar0 = 0,betaL = 0,phistar0 = 0,betaphi = 0,param = "fifth",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalpha,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filt = filt,style = "y")
			

		#shortens to use only positive values
		z_Lymanalpha_array = z_Lymanalpha_array[numpy.where(comovingphiarrayLymanalpha>0)]
		comovingphiarrayLymanalpha = comovingphiarrayLymanalpha[numpy.where(comovingphiarrayLymanalpha>0)]
		# print("z_Lymanalpha_array",z_Lymanalpha_array)
		# print("comovingphiarrayLymanalpha",comovingphiarrayLymanalpha)

		zFWHMarray_Lymanalpha_NOT = z_Lymanalpha_array[numpy.where(z_Lymanalpha_array>=zFWHMlow)]
		comovingphiarrayLymanalpha_FWHM_NOT = comovingphiarrayLymanalpha[numpy.where(z_Lymanalpha_array>=zFWHMlow)]
		# print("numpy.where(z_Lymanalpha_array>=zFWHMlow):",numpy.where(z_Lymanalpha_array>=zFWHMlow))

		zFWHMarray_Lymanalpha = zFWHMarray_Lymanalpha_NOT[numpy.where(zFWHMarray_Lymanalpha_NOT<=zFWHMhigh)]
		comovingphiarrayLymanalpha_FWHM = comovingphiarrayLymanalpha_FWHM_NOT[numpy.where(zFWHMarray_Lymanalpha_NOT<=zFWHMhigh)]
		# print("numpy.where(zFWHMarray_Lymanalpha_NOT<=zFWHMhigh):",numpy.where(zFWHMarray_Lymanalpha_NOT<=zFWHMhigh))

		# print("zFWHMarray_Lymanalpha = ",zFWHMarray_Lymanalpha)
		# print("comovingphiarrayLymanalpha_FWHM_NOT = ",comovingphiarrayLymanalpha_FWHM_NOT)
		# print("comovingphiarrayLymanalpha_FWHM = ",comovingphiarrayLymanalpha_FWHM)
		# print("zFWHMlow",zFWHMlow)
		# print("zFWHMhigh",zFWHMhigh)


		#adds legend

		firstz = z_Lymanalpha_array[0]
		lastz = z_Lymanalpha_array[-1]
		zpaper = r"Ly$\alpha$ z = "+str(round(firstz,2))+"-"+str(round(lastz,2))+" Ciardullo+ 2012"

		figure(1)
		plot(0,0,"y",label=zpaper)
		legend(loc="upper right")

		figure(2)
		plot(0,0,"y",label=zpaper)
		legend(loc="upper right")

		show()


		#print("Lymanalpha:")

		#I need the total number of galaxies WITHIN the FWHM
		#1 - integrate comovingphi using trapz and the comoving volume - will get total comoving number
		#2 - get DECam areal/sky density from this - multiply by 3 square degrees to get final answer
		#or skip areal number density and just multiply the result from the comovingphi integration by 3/42000 to get the DECam area

		#first convert z_Lymanalpha_array to a comoving volume array
		comovingvolarrayLymanalpha = numpy.zeros(len(z_Lymanalpha_array))

		#then multiply phi by the comoving volume array to get total number of galaxies expected
		#the following finds the comoving volume array to do so
		for r in range(len(z_Lymanalpha_array)):
			temparray = cosmo.comoving_volume(z_Lymanalpha_array[r]) #units are in Mpc^3
			#change tuple to value
			temparray = temparray.value
			comovingvolarrayLymanalpha[r] = temparray #this is in Mpc^3


		#DECam area of the sky is 3/42000 of the whole sky

		#this is ending up as a weird number


		#print("comovingphiarrayLymanalpha = ",comovingphiarrayLymanalpha)
		#print("comovingvolarrayLymanalpha = ",comovingvolarrayLymanalpha)
		#print("z_Lymanalpha_array = ",z_Lymanalpha_array)


		#integrates comovingphi over whole sky
		total_number_Lymanalpha = numpy.trapz(comovingphiarrayLymanalpha,x=comovingvolarrayLymanalpha)
		#print("total_number_FWHM_Lymanalpha = ",total_number_FWHM_Lymanalpha)
		#DECam area is 3/42000
		total_number_DECamLymanalpha = total_number_Lymanalpha*3/42000


		#want average number density for gaussian-approximated filters, so will use from FWHM of each filter only
		#use FWHMlow_actual and FWHMhigh_actual, which are wavelength values
		average_number_density_FWHM = numpy.sum(comovingphiarrayLymanalpha_FWHM)/len(comovingphiarrayLymanalpha_FWHM)


	else:
		zLymanalpha = 0




	print("end Lymanalpha")



	print("the total expected number of galaxies in the whole area of the sky covered by DECam is:")


	if len(z_Lymanalpha_array)>0:

		print("total_number_DECamLymanalpha = ",total_number_DECamLymanalpha)
		print("average_number_density_FWHM = ",average_number_density_FWHM)
		#arealphiLymanalpha = total_number_DECam_Lymanalpha/(4*numpy.pi)
		#print("arealphiLymanalpha =",arealphiLymanalpha,"steradian^-1")




show()