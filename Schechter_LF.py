print("This code plots Schechter Luminosity Functions for different emission lines in different papers.  ")
print("Comparat et al 2016 plots [OII] 3726/3729, Hbeta 4861, and [OIII] 5007")
print("[OIII] 5007 is always 3 times stronger than [OIII] 4959")
print("Sobral et al 2013 plots Halpha 6563")
print("Ciardullo et al 2012 plots Lymanalpha 1216")

print("(units are in Angstroms):")
print("[OII] 3726/3729 unresolved doublet, Hbeta 4861, [OIII] 4959/5007, Halpha 6563, Lymanalpha 1216")
emline = input("Plot Schechter Luminosity Function for which emission line?  (sample input: Halpha)")
print("emline = ", emline)


import numpy
from numpy import loadtxt
#should reorganize this somehow
from matplotlib.pyplot import *
import scipy
from scipy import integrate
from astropy.cosmology import Planck15 as cosmo


#this is the stuff that I need to do now: (actually updated)
#integrate the comovingphi array wrt the distance using astropy.cosmology (see end)


#the code will be able to use any filter
print("The ugrizy filter options with corresponding final coadded depths (5 sigma) are as follows:")
print("u: 26.1, g: 27.4, r: 27.5, i: 26.8, z: 26.1, y: 24.9")
filt = input("Plot Schechter Luminosity Function for which filter?  (sample input: zband)")
print("filter = ", filt)


#THIS PART CONTAINS FUNCTIONS THAT WILL THEN BE USED LATER ON IN THE CODE


def lineqLya(z):

	#this is used for the Lymanalpha line
	#finds linear equation from two points in Ciardullo et al. 2012 paper and interpolates/extraoplates to find Lstar and phistar for different z values

	z1 = 3.113
	z2 = 2.063
	deltaz = z1 - z2

	log10Lstar1 = 42.76
	log10Lstar2 = 42.33

	log10phistar1 = -3.17
	log10phistar2 = -2.86

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

	print("log10Lstar for Lymanalpha: log10Lstar = z*",mlog10Lstar," + ",blog10Lstar)
	print("log10phistar for Lymanalpha: log10phistar = z*",mlog10phistar," + ",blog10phistar)

	log10Lstar = mlog10Lstar*z + blog10Lstar
	log10phistar = mlog10phistar*z + blog10phistar

	LstarLya = 10**log10Lstar
	phistarLya = 10**log10phistar

	print("LstarLya =",LstarLya)
	print("phistarLya =",phistarLya)

	answers = [LstarLya,phistarLya]

	return answers


def filter_int(filt):

	if filt=="uband":

		#read in each data file
		print("u_filter has wavelengths 305.30 - 408.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		u_filter = loadtxt('ufilteredit.csv')
		print(u_filter)
		#I shorten this to only the second column
		LSST_filter = u_filter
		LSSTfilter = u_filter[:,1]

	if filt=="gband":
		#read in each data file
		print("g_filter has wavelengths 386.30 - 567.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		g_filter = loadtxt('gfilteredit.csv')
		print(g_filter)
		#I shorten this to only the second column
		LSST_filter = g_filter
		LSSTfilter = g_filter[:,1]

	if filt=="rband":

		#read in each data file
		print("r_filter has wavelengths 536.90 - 706.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput 
		r_filter = loadtxt('rfilteredit.csv')
		print(r_filter)
		#I shorten this to only the second column
		LSST_filter = r_filter
		LSSTfilter = r_filter[:,1]

	if filt=="iband":

		#read in each data file
		print("i_filter has wavelengths 675.90 - 833.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		i_filter = loadtxt('ifilteredit.csv')
		print(i_filter)
		#I shorten this to only the second column
		LSST_filter = i_filter
		LSSTfilter = i_filter[:,1]

	if filt=="zband":

		#read in each data file 
		print("z_filter has wavelengths 802.90 - 938.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		z_filter = loadtxt('zfilteredit.csv')
		print(z_filter)
		#I shorten this to only the second column
		LSST_filter = z_filter
		LSSTfilter = z_filter[:,1]

	if filt=="yband":

		#read in each data file
		print("y_filter has wavelengths 908.30 - 1099.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		y_filter = loadtxt('yfilteredit.csv')
		print(y_filter)
		#I shorten this to only the second column
		LSST_filter = y_filter
		LSSTfilter = y_filter[:,1]

	#to find the midpoint, integrate over the entire filter, then find what value would give you half the integrated value
	#see if you need to fix this because should have x=
	stepint_LSSTfilter = scipy.integrate.cumtrapz(LSSTfilter)
	print('stepint_LSSTfilter')
	print(stepint_LSSTfilter)
	len(stepint_LSSTfilter)
	#len gives total length, which is index number + 1
	lastnumber = stepint_LSSTfilter[len(stepint_LSSTfilter)-1]
	midpoint = lastnumber/2

	#finds closest value by finding minimum difference
	difference = abs(midpoint-stepint_LSSTfilter)
	print("difference = ",difference)
	mindiff = min(difference)
	print("mindiff = ",mindiff)
	midindex = numpy.where(difference==mindiff)
	print("midindex = ", midindex)
	print("type is ",type(midindex))
	#have to change this from tuple to float - first find the actual value in index 0
	midindex = midindex[0]
	midindex = numpy.float(midindex)
	print("midindex = ", midindex)

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
	print("centerindex = ",centerindex)
	print(type(centerindex))
	centerval = LSST_filter[centerindex]
	print(type(centerval))
	lambdacenter = centerval[0]
	lambdacenter = numpy.float(lambdacenter)
	print("median transmission wavelength of the "+filt+" = ",lambdacenter)

	#the central wavelength is actually called the median transmission wavelength of each band


	#want to find the FWHM of each filter

	#finds the index of the maximum value and the corresponding wavelength
	maxval = max(LSSTfilter)
	maxindex = numpy.where(LSSTfilter==maxval)
	print(type(maxindex))
	#want to change maxindex from tuple to float - just take value
	maxindex = maxindex[0]
	maxindex = numpy.float(maxindex)
	maxindex = numpy.int(maxindex) #attempt at fixing an error
	LSSTlambda = LSST_filter[:,0]
	maxlambda = LSSTlambda[maxindex]
	print("maxlambda = ",maxlambda)

	#need the two points of the filter whose values give half the maximum value

	#checks first and last fourth of iband filter for half max value because of weird iband shape
	if filt == "iband":
		halfmaxval = maxval/2
		fourth = (len(LSSTfilter))/4
		print(fourth*4,fourth)
		fourthapprox = round(fourth)
		LSSTfilterfourth1 = LSSTfilter[:fourthapprox]
		LSSTfilterfourth4 = LSSTfilter[fourthapprox*2:]
		left = LSSTfilterfourth1
		right = LSSTfilterfourth4
		#finds closest value by finding minimum difference - need left and right parts
		diffleft = abs(halfmaxval-LSSTfilterfourth1)
		diffright = abs(halfmaxval-LSSTfilterfourth4)

	#check each half of the filter for half max value
	if filt != "iband":
		halfmaxval = maxval/2
		half = (len(LSSTfilter))/2
		print(half*2,half)
		halfapprox = round(half)
		LSSTfilterhalf1 = LSSTfilter[:halfapprox]
		LSSTfilterhalf2 = LSSTfilter[halfapprox:]
		left = LSSTfilterhalf1
		right = LSSTfilterhalf2
		#finds closest value by finding minimum difference - need left and right parts
		diffleft = abs(halfmaxval-LSSTfilterhalf1)
		diffright = abs(halfmaxval-LSSTfilterhalf2)

	closestleft = min(diffleft)
	closestright = min(diffright)
	print("closestleft = ",closestleft)
	print("closestright = ",closestright)
	indexL = numpy.where(diffleft==closestleft)
	indexR = numpy.where(diffright==closestright)
	print("indexL = ", indexL)
	print("indexR = ", indexR)
	#have to change this from tuple to int - first find the actual value in index 0
	indexL = indexL[0]
	indexR = indexR[0]
	indexL = numpy.int(indexL)
	indexR = numpy.int(indexR)
	print("indexL = ",indexL)
	print("indexR = ",indexR)

	#find value in left and right sections of the filter
	valueleft = left[indexL]
	valueright = right[indexR]
	print(valueleft,valueright)

	#match value to LSSTfilter
	LSSThalfmaxleft = numpy.where(LSSTfilter==valueleft)
	LSSThalfmaxright = numpy.where(LSSTfilter==valueright)
	print(LSSThalfmaxleft,LSSThalfmaxright)
	print(type(LSSThalfmaxleft))
	print(type(LSSThalfmaxright))
	#want to change from tuple to float - just take value
	indexleftFWHM = LSSThalfmaxleft[0]
	indexrightFWHM = LSSThalfmaxright[0]
	indexleftFWHM = numpy.int(indexleftFWHM)
	indexrightFWHM = numpy.int(indexrightFWHM)

	#match to lambda
	#LSSTlambda = LSST_filter[:,0]
	lambdaFWHMleft = LSSTlambda[indexleftFWHM]
	lambdaFWHMright = LSSTlambda[indexrightFWHM]
	print(lambdaFWHMleft,lambdaFWHMright)

	#should have output as an array (see LF code)
	print("returns array: lambdacenter,lambdaFWHMleft,lambdaFWHMright")
	lambdaarray = [lambdacenter,lambdaFWHMleft,lambdaFWHMright]

	return lambdaarray


#the following function calculates the luminosity limit for a certain band with a certain detection limit
def lumlim(z,em,filt):

	#REMINDER TO SELF - z is a function of the emission line wavelength and the filter itself, so it does not have to be adjusted for the FWHM thing

	print("This function will calculate the luminosity that corresponds to a 5 sigma detection, in erg/s:")
	#print("INFO: 26.2 is AB magnitude for 5 sigma detection limit in the z band") #this was used when I only had the z band
	#first I calculate the flux, then convert to flux density, then find the luminosity limit for the conditions printed above

	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	ABmag_dict = {"uband":26.1,"gband":27.4,"rband":27.5,"iband":26.8,"zband":26.1,"yband":24.9}
	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm
	#these are at the endpoints of the filters, NOT the FWHM

	ABmag = ABmag_dict[filt]
	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]

	#this is how the code previously was before I changed it to using dictionaries

	#if filt=="uband":
	#	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	#	ABmag = 26.1
	#	lambdalow = 305.30 #in nm
	#	lambdahigh = 408.60 #in nm

	#if filt=="gband":
	#	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	#	ABmag = 27.4
	#	lambdalow = 386.30 #in nm
	#	lambdahigh = 567.00 #in nm

	#if filt=="rband":
	#	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	#	ABmag = 27.5
	#	lambdalow = 536.90 #in nm
	#	lambdahigh = 706.00 #in nm

	#if filt=="iband":
	#	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	#	ABmag = 26.8
	#	lambdalow = 675.90 #in nm
	#	lambdahigh = 833.00 #in nm

	#if filt=="zband":
	#	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	#	ABmag = 26.1
	#	lambdalow = 802.90 #in nm
	#	lambdahigh = 938.60 #in nm

	#if filt=="yband":
	#	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	#	ABmag = 24.9
	#	lambdalow = 908.30 #in nm
	#	lambdahigh = 1099.60 #in nm

	#might need the following here:

	if filt=="uband":

		#read in each data file
		print("u_filter has wavelengths 305.30 - 408.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		u_filter = loadtxt('ufilteredit.csv')
		print(u_filter)
		#I shorten this to only the second column
		LSST_filter = u_filter
		LSSTfilter = u_filter[:,1]

	if filt=="gband":
		#read in each data file
		print("g_filter has wavelengths 386.30 - 567.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		g_filter = loadtxt('gfilteredit.csv')
		print(g_filter)
		#I shorten this to only the second column
		LSST_filter = g_filter
		LSSTfilter = g_filter[:,1]

	if filt=="rband":

		#read in each data file
		print("r_filter has wavelengths 536.90 - 706.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput 
		r_filter = loadtxt('rfilteredit.csv')
		print(r_filter)
		#I shorten this to only the second column
		LSST_filter = r_filter
		LSSTfilter = r_filter[:,1]

	if filt=="iband":

		#read in each data file
		print("i_filter has wavelengths 675.90 - 833.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		i_filter = loadtxt('ifilteredit.csv')
		print(i_filter)
		#I shorten this to only the second column
		LSST_filter = i_filter
		LSSTfilter = i_filter[:,1]

	if filt=="zband":

		#read in each data file 
		print("z_filter has wavelengths 802.90 - 938.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		z_filter = loadtxt('zfilteredit.csv')
		print(z_filter)
		#I shorten this to only the second column
		LSST_filter = z_filter
		LSSTfilter = z_filter[:,1]

	if filt=="yband":

		#read in each data file
		print("y_filter has wavelengths 908.30 - 1099.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		y_filter = loadtxt('yfilteredit.csv')
		print(y_filter)
		#I shorten this to only the second column
		LSST_filter = y_filter
		LSSTfilter = y_filter[:,1]



	#THIS PART IS OKAY

	#the following calculates values I use in and plug into the lumlim function

	#finds the flux density
	print("ABmagnitude = -2.5*log10(fluxdensity/(3631 Jansky))")
	print("consequently:")
	#print("fluxdensity = (10**(ABmagnitude/(-2.5)))*(3631 Janksy)")
	#fluxdens = (10**(ABmag/(-2.5)))*3631 #outputs in Jansky
	#uses ABmag from earlier in this function
	fluxdens = 10**((ABmag+48.6)/(-2.5)) #outputs in erg/(s*Hz*(cm^2))
	print("flux density =",fluxdens,"erg/(s*Hz*(cm^2))")



	#THE FOLLOWING NEEDS TO BE MADE MORE ACCURATE:  (everything before/after is okay)

	#finds the flux using the difference between the frequencies at each end of the band
	c = 2.9979*(10**17) #in nm/s
	deltanu = c*((1/lambdalow)-(1/lambdahigh)) #the nm should cancel out
	print("deltanu =",deltanu,"s^-1")
	flux = fluxdens*deltanu#*(10**(-23)) #the extra factor converts from Janskys to ergs/(s*Hz*(cm^2)) #THEN WHAT ARE MY UNITS????
	print("flux =",flux,"erg/(s*(cm^2))")  #IS THIS EVEN RIGHT THEN, IF I LOOK AT THE ABOVE LINE????  REVIEW THE UNITS!!  



	#THE FOLLOWING IS GOOD:

	#finds the luminosity distance
	lumdist = cosmo.luminosity_distance(z)
	#this outputs a special object that keeps track of units, so first I convert it to cm (cgs units), and then I convert it to a regular number
	lumdist_cgs = lumdist.to('cm')
	lumdist_unitless = lumdist_cgs.value
	print("the luminosity distance for redshift z =",z,"is lumdist =",lumdist_unitless,"cm")

	#finds the luminosity limit
	print("Luminosity = 4*pi*(luminositydistance**2)*flux")
	lumlimit = 4*numpy.pi*(lumdist_unitless**2)*flux
	print("luminosity limit for 5 sigma detection of",em,"in "+filt+" band is",lumlimit,"ergs/s")

	#using return makes the main output of this function the value of lumlimit so that I can use it to calculate other things when I call this function
	return lumlimit


#the following defines and plots the Schechter luminosity function for a chosen emission line along with a separate plot of the number density
#it also calculates the number density above a luminosity limit calculated for a specific band in the previous function
#to use this function, type schechter_LF(z=redshifttoplot) or with additional parameters you want to change inside the ()
def schechter_LF(z,lambdaemitted,alpha,Lstar0,betaL,phistar0,betaphi,zpaper,param,fluxscale,em,filt,style = ""):

	#REMINDER TO SELF - z is a function of the emission line wavelength and the filter itself, so it does not have to be adjusted for the FWHM thing 

	print("alpha = ", alpha)
	print("z = ", z)

	log10Lstararray = numpy.arange(30,55,0.01)
	L = (10**log10Lstararray)

	#I have two different parametrizations for Lstar and phistar; the first one is from Comparat et al 2016, and the second one is from Sobral et al 2015
	#I added a third one when I linearly parametrized Lymanalpha from Ciardullo et al 2012, but I am still trying to get it to work

	if param == "first":

		Lstar = fluxscale*Lstar0*((1+z)**betaL)
		print("Lstar = ", Lstar)

		phistar = phistar0*((1+z)**betaphi)
		print("phistar = ", phistar)

	if param == "second":

		Lstar = 10**(0.45*z + Lstar0)
		print("Lstar = ", Lstar)

		phistar = 10**(-0.38*(z**2) + z + phistar0)
		print("phistar = ", phistar)

	if param == "third":

		answersLya = lineqLya(z=z)
		print(type(answersLya))

		Lstar = answersLya[0] #this is linearly parametrized from Ciardullo+ 2012
		print("Lstar = ",Lstar)
		print(type(Lstar))

		phistar = answersLya[1] #this is linearly parametrized from Ciardullo+ 2012
		print("phistar = ",phistar)
		print(type(phistar))


	#now need an array of phi for each z or lambda step in that array within the FWHM
	phi = phistar*((L/Lstar)**(alpha+1))*(numpy.e**(-L/Lstar))

	#this deletes parts of the arrays that are so small python counts them as zero; otherwise, I would not be able to take the logarithm of the array
	L = L[numpy.where(phi!=0)]
	log10Lstararray_nonzero = log10Lstararray[numpy.where(phi!=0)]
	phi = phi[numpy.where(phi!=0)]

	log10L = numpy.log10(L)
	log10phi = numpy.log10(phi) #the error traces back to here

	#the following calculates values I use in and plug into the lumlim function

	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm

	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	FWHMlow = lambdaarray[1]
	FWHMhigh = lambdaarray[2]

	#once again, changed below to dictionaries above

	#if filt=="uband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.1
	#	lambdalow = 305.30 #in nm
	#	lambdahigh = 408.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="gband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 27.4
	#	lambdalow = 386.30 #in nm
	#	lambdahigh = 567.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="rband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 27.5
	#	lambdalow = 536.90 #in nm
	#	lambdahigh = 706.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="iband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.8
	#	lambdalow = 675.90 #in nm
	#	lambdahigh = 833.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="zband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.1
	#	lambdalow = 802.90 #in nm
	#	lambdahigh = 938.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="yband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 24.9
	#	lambdalow = 908.30 #in nm
	#	lambdahigh = 1099.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#uses previously defined function to get median transmission wavelenth
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	#lambdaemitted = #will be one of the three emission lines #uses nm #these are redshifts
	filterendlow = (lambdalow/lambdaemitted)-1
	filterendhigh = (lambdahigh/lambdaemitted)-1
	filtercenter = (lambdacenter/lambdaemitted)-1
	#I don't think I used the following at all
	#deltafilter = filterendhigh-filterendlow
	#print("for",em,"deltafilter =",deltafilter)

	#finds luminosity limit in center of z band
	center = lumlim(z = filtercenter,em = em,filt = filt)

	#sets up an array that I use to plot
	lumarray = numpy.full(15,numpy.log10(center))
	yarray = numpy.arange(-10,5,1)

	figure(1)
	plot(log10L,log10phi,style,label = zpaper)
	plot(lumarray,yarray,style+"--")
	xlim(40,45)
	ylim(-7,0)
	legend(loc = "upper right")
	#LaTeX is used with $ signs in the titles below
	xlabel("$\log_{10}(L [ergs/s])$")
	ylabel("$\log_{10}(\phi [\t{Mpc}^{-3}])$")
	title("Schechter Luminosity Function, "+filt)
	print("z = ", z, " and alpha = ", alpha, " are plotted, ")
	#saves image
	#pyplot.savefig('/home/lanaeid/Desktop/fig1'+filter+'.png',bbox_inches = 'tight')

	print("now the number density is calculated by integrating the LF using cumtrapz:")
	#cumptrapz integrates the opposite way than I need to integrate, so I flip it twice in the process
	phiflip = phi[::-1]
	print(phi)
	print(len(phiflip),len(log10Lstararray_nonzero))
	#always pay attention to the version of the LF that is being used
	phiflipint = scipy.integrate.cumtrapz(phiflip,x=log10Lstararray_nonzero)
	num_dens = phiflipint[::-1]

	Lmin = numpy.delete(L,(len(L)-1))
	#shifted_Lmin = Lmin*fluxscale
	log10Lmin = numpy.log10(Lmin)
	#log10shifted_Lmin = numpy.log10(shifted_Lmin)

	figure(2)
	plot(log10Lmin,num_dens,style,label = zpaper)
	plot(lumarray,yarray,style+"--")
	legend(loc = "upper right")
	xlim(40,44)
	ylim(0,0.03)
	#LaTeX is used with $ signs in the titles below
	xlabel("$\log_{10}(L_{min} [ergs/s])$")
	ylabel("$\log_{10}(\phi [\t{Mpc}^{-3}])$")
	title("Number Density, "+filt)
	#saves image
	#pyplot.savefig('/home/lanaeid/Desktop/fig2'+filter+'.png',bbox_inches = 'tight')

	#the following finds the comoving number density above a certain detection limit (shot noise limit),
	#which I calculated in the previous part of the code (also shown below), where I saved the variable named center

	#runs through the LF code for chosen emission line
	#then gets number density above the luminosity limit using the center value

	#first, I have to shorten the phi array to contain only values above the luminosity limit

	#EDIT HERE
	#to be able to use for not just median transmission wavelength
	#want to find expected number density within FWHM - eventually want to integrate over entire filter, but for now gradually stepping up in improvements

	#OKAY I THINK I DID THIS WRONG
	#what ended up happening is I was essentially integrating over the FWHM twice, resulting in a really weird number
	#I am going to comment out the parts of the code in the Schechter_LF function below that are unnecessary
	#then I am going to see if it fixes the issue - if this was not even used before for the center, I should delete it

	#first define array of 100 evenly spaced wavelengths

#	FWHMlow = lambdaarray[1]
#	FWHMhigh = lambdaarray[2]
#	lambdaFWHMarray = numpy.linspace(FWHMlow,FWHMhigh,num=100)
	#print("lambdaFWHMarray",lambdaFWHMarray)
	#print("LOOK HERE",len(lambdaFWHMarray))

	#use these to get an array of redshifts with which to get the comoving volume, for each of them

	#defines stuff to use in next part
	lambda_OII = 372.7
	lambda_OIII = 500.7
	lambda_Halpha = 656.3
	lambda_Lymanalpha = 121.6

	#finds lambda corresponding to emission line
	#FIX THIS - did not adjust for FWHM calculation, just sets it to the same number - incorrect
	#make sure this does not happen anywhere else
	if z==zOII:
		lambdaem = lambda_OII
	elif z==zOIII:
		lambdaem = lambda_OIII
	elif z==zHalpha:
		lambdaem = lambda_Halpha
	elif z==zLymanalpha:
		lambdaem = lambda_Lymanalpha

#	zFWHMarray = numpy.zeros(100)
	#print("LOOK HERE",len(zFWHMarray))

	#constructs z array from wavelength array within FWHM
#	for l in range(len(lambdaFWHMarray)):
#		zFWHMarray[l] = (lambdaFWHMarray[l]/lambdaem)-1
	#print(zFWHMarray)

	#ASK IF THIS MAKES SENSE TO DO
	#don't want negative redshifts
#	zFWHMarray = zFWHMarray[numpy.where(zFWHMarray>=0)]

	#basically getting the thing I found below but for everything across the FWHM

	#set up empty array for comoving volumes in evenly spaced wavelength intervals
	#comovingvolarray= numpy.zeros(100) unless positive redshifts
#	if len(zFWHMarray)!=0:
#		comovingvolarray = numpy.zeros(len(zFWHMarray))

		#now find comoving volumes
		#need to "integrate" or sum because not linear as function of redshift of lambda
		#can just ignore everything with redshift 0 (bc cosmological volume will be 0)
#		for r in range(len(zFWHMarray)):
#			partr = cosmo.comoving_volume(zFWHMarray[r]) #units are in Mpc^3
			#change tuple to value
#			partr = partr.value
#			comovingvolarray[r] = partr

#	if len(zFWHMarray)==0:
#		print("could not find comoving volume because not in this filter")

	#NOW JUST NEED TO INTEGRATE THE ARRAY
	#comoving vol is a function of z, so i should integrate it wrt dz?  so i need the zFWHMarray
	#ONLY WORKS IF COMOVINGVOLARRAY IS PER DELTA Z
	#JUST SUM COMOVINGVOLARRAY
	#numpy.sum(numpy.diff(array))  ->  replace the next line with this (is commented out bc wrong)
	#comovingvol = numpy.trapz(comovingvolarray,x=zFWHMarray) #integrated
#	comovingvol = numpy.sum(numpy.diff(comovingvolarray))
#	print("comovingvol =",comovingvol,"Mpc^3")

	#COMMENTED THE NEXT PART OUT, WAS INCORRECT ANYWAYS?
	#I initially commented this out, because I did not have the correct values, and I was repurposing the code for the FWHM
	#now I need it again to see if I am calculating the FWHM number correctly
	#uses astropy.cosmology.Planck15 to find comoving volume in shell between redshifts at each end of the z filter
	#note: DO NOT USE FOR FWHM CALCULATION, THIS IS ONLY FOR THE INITIAL ESTIMATE AT THE MEDIAN TRANSMISSION WAVELENGTH
	comovingvolmin = cosmo.comoving_volume(filterendlow) #units are in Mpc^3
	comovingvolmax = cosmo.comoving_volume(filterendhigh) #units are in Mpc^3
	comovingvol = comovingvolmax-comovingvolmin
	comovingvol = comovingvol.value
	print("comovingvol =",comovingvol,"Mpc^3")

	#THIS IS ALL I NEED:

	#shortens the array to be above the luminosity limit, then integrates to get comoving number density
	philim = phi[numpy.where(L>center)]
	#log10Lstararray = numpy.arange(30,55,0.01)
	#L = (10**log10Lstararray)
	log10Lstararray_lumlim = log10Lstararray_nonzero[numpy.where(L>center)]
	comovingphi = scipy.integrate.trapz(philim,x=log10Lstararray_lumlim)
	print("comovingphi =",comovingphi,"Mpc^-3")

	#finds total number of galaxies and areal number density
	#eventually this should be done for the part of the sky that LSST observes, not the entire sky
	#ONCE AGAIN UNCOMMENTING LINES I COMMENTED OUT BEFORE
	totalnumgalaxies = comovingphi*comovingvol

#	totalnumgalaxiesarray = comovingphi*comovingvolarray
	#changed numpy.trapz to numpy.sum, which makes more sense here
#	totalnumgalaxies = numpy.sum(totalnumgalaxiesarray)

	#the total number of galaxies here should be slightly less than the one calculated all the way at the end of the code for the FWHM
	#however..it is much larger..
#	print("OVER HERE LANA totalnumgalaxies = ",totalnumgalaxies)
#	arealphi = totalnumgalaxies/(4*numpy.pi)
#	print("arealphi =",arealphi,"steradian^-1")

	show()

	#expected number of galaxies in LSST area of the sky
	totalnumgalaxiesLSST = totalnumgalaxies*18000./42000.
	print("OVER HERE LANA")
	print(em+"totalnumgalaxiesLSST = ",totalnumgalaxiesLSST)

	print("comovingphi is returned for the input")
	
	return comovingphi



#THIS NEXT PART CONTAINS IF STATEMENTS FOR EACH EMLINE


#functionality of code has changed - commenting these out, bc not using them anymore - also editing part for median transmission wavelength above to use FWHM


#if emline == "[OII]":

#	print("The following plot will contain [OII] LF's from Comparat et al 2016:")

#	print("[OII] values with errors given in the Comparat et al 2016 paper are:")
#	print("log10Lstar =   41.1  (+0.02 -0.03)")
#	print("betaL =         2.33 (+0.14 -0.16)")
#	print("log10phistar = -2.4  (+0.03 -0.03)")
#	print("betaphi =      -0.73 (+0.25 -0.29)")
#	print("alpha =        -1.46 (+0.06 -0.05)")

#	schechter_LF(z=0.75,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "z = 0.75 Comparat+ 2016",fluxscale = 1)
#	schechter_LF(z=1.5,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "z = 1.5 Comparat+ 2016",fluxscale = 1)

	#The Khostovan paper is not useful, as it does not actually parametrize Lstar and phistar

	#print("[OII] values with errors given in the Khostovan et al 2015 paper are:")
	#print("log10Lstar =   41.86 (+0.03 -0.03)")
	#print("log10phistar = -2.25 (+0.04 -0.04)")
	#print("alpha =        -1.3")

	#schechter_LF(z=1.47,alpha = -1.3,Lstar0 = 10**41.86,betaL = 0,phistar0 = 10**(-2.25),betaphi = 0,param = "first",zpaper = "z = 1.47 Khostovan+ 2015",fluxscale = 1)


#if emline == "Hbeta":

#	print("The following plot will contain Hbeta LF's from Comparat et al 2016:")

#	print("Hbeta values with errors given in the Comparat et al 2016 paper are:")
#	print("log10Lstar =   40.88 (+0.05 -0.07)")
#	print("betaL =         2.19 (+0.25 -0.32)")
#	print("log10phistar = -3.34 (+0.09 -0.12)")
#	print("betaphi =       2.7  (+0.44 -0.57)")
#	print("alpha =        -1.51 (+0.27 -0.2)")

#	schechter_LF(z=0.45,alpha = -1.51,Lstar0 = 10**40.88,betaL = 2.19,phistar0 = 10**(-3.34),betaphi = 2.7,param = "first",zpaper = "z = 0.45 Comparat+ 2016",fluxscale = 1)
#	schechter_LF(z=0.9,alpha = -1.51,Lstar0 = 10**40.88,betaL = 2.19,phistar0 = 10**(-3.34),betaphi = 2.7,param = "first",zpaper = "z = 0.9 Comparat+ 2016",fluxscale = 1)


#if emline == "[OIII]":

#	print("The following plot will contain [OIII] LF's from Comparat et al 2016:")

#	print("[OIII] values with errors given in the Comparat et al 2016 paper are:")
#	print("log10Lstar =   41.42 (+0.07 -0.09)")
#	print("betaL =         3.91 (+0.32 -0.4)")
#	print("log10phistar = -3.41 (+0.08 -0.01)")
#	print("betaphi =      -0.76 (+0.39 -0.49)")
#	print("alpha =        -1.83 (+0.1 -0.08)")

#	schechter_LF(z=0.45,alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "z = 0.45 Comparat+ 2016",fluxscale = 1)
#	schechter_LF(z=0.9,alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "z = 0.9 Comparat+ 2016",fluxscale = 1)


#if emline == "Halpha":

#	print("This following plot will contain an Halpha LF from Sobral et al 2013.  ")

#	print("z = 0.4 and alpha = -1.6 are plotted, ")
#	print("but the values with errors given in the Sobral et al 2013 paper are:")
#	print("z = 0.40 +_ 0.01")
#	print("alpha =       -1.75 (+0.12 -0.08)")
#	print("Lstar_exp =   41.95 (+0.47 -0.12)")
#	print("phistar_exp = -3.12 (+0.10 -0.34)")

#	schechter_LF(z = 0.4,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = "z = 0.4 Sobral+ 2013",fluxscale = 1)
#	schechter_LF(z = 0.45,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = "z = 0.4 Sobral+ 2013",fluxscale = 1)	

#	print("This will now plot an Halpha LF scaled from the Hbeta LF in this code.  ")
#	print("Halpha/Hbeta = 2.86")

#	schechter_LF(z=0.45,alpha = -1.51,Lstar0 = 10**40.88,betaL = 2.19,phistar0 = 10**(-3.34),betaphi = 2.7,param = "first",zpaper = "z = 0.45 Comparat+ 2016",fluxscale = 2.86,style = "--")
#	schechter_LF(z=0.9,alpha = -1.51,Lstar0 = 10**40.88,betaL = 2.19,phistar0 = 10**(-3.34),betaphi = 2.7,param = "first",zpaper = "z = 0.9 Comparat+ 2016",fluxscale = 2.86,style = "--")


#if emline == "Lymanalpha":

#	print("This following plot will contain a Lymanalpha LF from Ciardullo et al 2012.  ")
 	#either this needs to be combined with the "test" option, or I need to edit the rest of these emline options to be usable with all the different parameters

#	lambda_OII = 372.7
#	lambda_OIII = 500.7
#	lambda_Halpha = 656.3
#	lambda_Lymanalpha = 121.6

	#the following calculates values I use in and plug into the lumlim function
#	lambdalow = 818.95 #in nm
#	lambdahigh = 921.15 #in nm
#	lambdacenter = (lambdalow+lambdahigh)/2 #in nm #NO, THIS IS NOR CORRECT, FIX IT
	#lambdaemitted = #will be one of the three emission lines #in nm
	#zendlow = (lambdalow/lambdaemitted)-1
	#zendhigh = (lambdahigh/lambdaemitted)-1
	#zcenter = (lambdacenter/lambdaemitted)-1

	#finds redshift of emission line at center of z band -> use to find LF
#	zOII = (lambdacenter/lambda_OII)-1
#	zOIII = (lambdacenter/lambda_OIII)-1
#	zHalpha = (lambdacenter/lambda_Halpha)-1
#	zLymanalpha = (lambdacenter/lambda_Lymanalpha)-1

	#for the first redshift below
#	print("z = 3.113 and alpha = -1.65 are plotted, ")
#	print("the values with errors given in the Ciardullo et al 2012 paper are:")
#	print("Lstar_exp =   42.76 (+0.10 -0.10)")
#	print("phistar_exp = -3.17 (+0.05 -0.05)")

	#for the second redshift below
#	print("z = 2.063 and alpha = -1.65 are plotted, ")
#	print("the values with errors given in the Ciardullo et al 2012 paper are:")
#	print("Lstar_exp =   42.33 (+0.12 -0.12)")
#	print("phistar_exp = -2.86 (+0.05 -0.05)")

	#schechter_LF(z=3.113,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**42.76,betaL = 0,phistar0 = 10**(-3.17),betaphi = 0,param = "first",zpaper = r"Ly$\alpha$ z = 3.113 Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",style = "c")
	#schechter_LF(z=2.063,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**42.33,betaL = 0,phistar0 = 10**(-2.86),betaphi = 0,param = "first",zpaper = r"Ly$\alpha$ z = 2.063 Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",style = "m")
	#the following is also from the separate linearequation code that I tried to put in this one, but it throws back errors every time I use the "third" option, so I have to fix that later
	#THE OPTION "first" ACTUALLY WORKS, BUT ONLY IF YOU ALSO INPUT THE VALUES FOR LSTAR0 AND PHISTAR0 TO BYPASS THE ACTUAL THING I WAS TRYING TO DO
#	schechter_LF(z=6.155016447368421,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**44.0057781641604,betaL = 0,phistar0 = 10**(-4.068119141604011),betaphi = 0,param = "first",zpaper = r"Ly$\alpha$ z = 6.155 Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filter = "zband",style = "y")


#if emline == "testHalpha":

#	print("This following plot will contain an Halpha LF from Sobral et al 2013.  ")

#	print("z = 0.4 and alpha = -1.6 are plotted, ")
#	print("but the values with errors given in the Sobral et al 2013 paper are:")
#	print("z = 0.40 +_ 0.01")
#	print("alpha =       -1.75 (+0.12 -0.08)")
#	print("Lstar_exp =   41.95 (+0.47 -0.12)")
#	print("phistar_exp = -3.12 (+0.10 -0.34)")

#	schechter_LF(z = 0.4,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = 10**(-3.18),betaphi = 0,param = "second",zpaper = "z = 0.4 Sobral+ 2013",fluxscale = 1)
	
#I have no idea why I wrote the following three lines to print out, so I'm commenting them out for now unless I have a reason to do otherwise
#print("Comparat et al 2016 plots [OII] 3726/3729, Hbeta 4861, and [OIII] 5007")
#print("[OIII] 5007 is always 3 times stronger than [OIII] 4959")
#print("Sobral et al 2013 plots Halpha 6563")


#might not be able to use this while editing the Schechter_LF routine for the allFWHM option, bc integrating in different order
#don't try it
if emline == "all":

	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm

	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	FWHMlow = lambdaarray[1]
	FWHMhigh = lambdaarray[2]

	#once again, changed below to dictionaries above

	#if filt=="uband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.1
	#	lambdalow = 305.30 #in nm
	#	lambdahigh = 408.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="gband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 27.4
	#	lambdalow = 386.30 #in nm
	#	lambdahigh = 567.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="rband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 27.5
	#	lambdalow = 536.90 #in nm
	#	lambdahigh = 706.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="iband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.8
	#	lambdalow = 675.90 #in nm
	#	lambdahigh = 833.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="zband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.1
	#	lambdalow = 802.90 #in nm
	#	lambdahigh = 938.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="yband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 24.9
	#	lambdalow = 908.30 #in nm
	#	lambdahigh = 1099.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#this was used to print out data and plots for the LSST-DESC conference poster

	print("This option finds comoving and areal number densities for:")
	print("[OII] 3726/3729 and [OIII] 5007 from Comparat et al 2016")
	print("Halpha 6563 from Sobral et al 2013")

	lambda_OII = 372.7
	lambda_OIII = 500.7
	lambda_Halpha = 656.3
	lambda_Lymanalpha = 121.6

	#the following calculates values I use in and plug into the lumlim function
	#THESE WERE ONLY FOR THE Z BAND - I AM TRYING TO MAKE IT MORE GENERAL - DELETE THESE LATER
	#lambdalow = 818.95 #in nm
	#lambdahigh = 921.15 #in nm
	#lambdacenter = (lambdalow+lambdahigh)/2 #in nm #NO, THIS IS NOT CORRECT, FIX THIS
	#lambdaemitted = #will be one of the three emission lines #in nm
	#zendlow = (lambdalow/lambdaemitted)-1
	#zendhigh = (lambdahigh/lambdaemitted)-1
	#zcenter = (lambdacenter/lambdaemitted)-1

	#finds redshift of emission line at center of z band -> use to find LF
	zOII = (lambdacenter/lambda_OII)-1
	zOIII = (lambdacenter/lambda_OIII)-1
	zHalpha = (lambdacenter/lambda_Halpha)-1
	zLymanalpha = (lambdacenter/lambda_Lymanalpha)-1

	#left side of FWHM
	#zOIIlow = (FWHMlow/lambda_OII)-1
	#zOIIIlow = (FWHMlow/lambda_OIII)-1
	#zHalphalow = (FWHMlow/lambda_Halpha)-1
	#zLymanalphalow = (FWHMlow/lambda_Lymanalpha)-1

	#right side of FWHM
	#zOIIhigh = (FWHMhigh/lambda_OII)-1
	#zOIIIhigh = (FWHMhigh/lambda_OIII)-1
	#zHalphahigh = (FWHMhigh/lambda_Halpha)-1
	#zLymanalphahigh = (FWHMhigh/lambda_Lymanalpha)-1

	#schechter_LF(z=zOII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "z = "+str(zOII)+" Comparat+ 2016",fluxscale = 1,style = "r")
	#schechter_LF(z=zOIII,alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "z = "+str(zOIII)+" Comparat+ 2016",fluxscale = 1,style = "g")
	#schechter_LF(z = zHalpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = "z = "+str(zHalpha)+" Sobral+ 2013",fluxscale = 1,style = "b")

	#THESE WERE THE ORIGINAL ONES: 
	#schechter_LF(z=zOII,lambdaemitted = lambda_OII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "[OII] z = 1.33 Comparat+ 2016",fluxscale = 1,em = "[OII]",filter = "zband",style = "r")
	#schechter_LF(z=zOIII,lambdaemitted = lambda_OIII, alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "[OIII] z = 0.74 Comparat+ 2016",fluxscale = 1,em = "[OIII]",filter = "zband",style = "g")
	#schechter_LF(z=zHalpha,lambdaemitted = lambda_Halpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = r"H$\alpha$ z = 0.33 Sobral+ 2013",fluxscale = 1, em = "Halpha",filter = "zband",style = "b")
	#the following is also from the separate linearequation code that I tried to put in this one, but it throws back errors every time I use the "third" option, so I have to fix that later
	#schechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**44.0057781641604,betaL = 0,phistar0 = 10**(-4.068119141604011),betaphi = 0,param = "first",zpaper = r"Ly$\alpha$ z = 6.155 Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filter = "zband",style = "y")
 	#used LaTex above for legend
 	#python recognizes LaTeX instead of thinking they're escape characters if I write it as a "raw" string, denoting it with an r in the beginning of the line

 	#CANNOT HAVE NEGATIVE REDSHIFTS - UNPHYSICAL FOR WHAT WE'RE DOING
 	
 	#NEED TO HAVE AN ARRAY OF Z VALUES TO INTEGRATE PHI????  

	if zOII>0:
		comovingphiOII = schechter_LF(z=zOII,lambdaemitted = lambda_OII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "[OII] z = "+str(round(zOII,2))+" Comparat+ 2016",fluxscale = 1,em = "[OII]",filt = filt,style = "r")
	
	if zOIII>0:
		comovingphiOIII = schechter_LF(z=zOIII,lambdaemitted = lambda_OIII, alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "[OIII] z = "+str(round(zOIII,2))+" Comparat+ 2016",fluxscale = 1,em = "[OIII]",filt = filt,style = "g")
	
	if zHalpha>0:
		comovingphiHalpha = schechter_LF(z=zHalpha,lambdaemitted = lambda_Halpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = r"H$\alpha$ z = "+str(round(zHalpha,2))+" Sobral+ 2013",fluxscale = 1, em = "Halpha",filt = filt,style = "b")
	
	if zLymanalpha>0:
		#schechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**44.0057781641604,betaL = 0,phistar0 = 10**(-4.068119141604011),betaphi = 0,param = "first",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalpha,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filter = filter,style = "y")
		#the following is also from the separate linearequation code that I tried to put in this one, but it throws back errors every time I use the "third" option, so I have to fix that later
 		#NOW THIS WORKS YESSSSS
 		comovingphiLymanalpha = schechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 0,betaL = 0,phistar0 = 0,betaphi = 0,param = "third",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalpha,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filt = filt,style = "y")
 	
 	#used LaTex above for legends
 	#python recognizes LaTeX instead of thinking they're escape characters if I write it as a "raw" string, denoting it with an r in the beginning of the line

 	#fix with if statements
 	#commented out bc not the way I set it up earlier
 	#print("comovingphiOII = ",comovingphiOII)
 	#print("comovingphiOIII = ",comovingphiOIII)
 	#print("comovingphiHalpha = ",comovingphiHalpha)
 	#print("comovingphiLymanalpha = ",comovingphiLymanalpha)

 	#now get exp num using comovingvol for ??????????????

 	#low end of FWHM

	#if zOIIlow>0:
	#	schechter_LF(z=zOIIlow,lambdaemitted = lambda_OII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "[OII] z = "+str(round(zOIIlow,2))+" Comparat+ 2016",fluxscale = 1,em = "[OII]",filter = filter,style = "r")
	
	#if zOIIIlow>0:
	#	schechter_LF(z=zOIIIlow,lambdaemitted = lambda_OIII, alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "[OIII] z = "+str(round(zOIIIlow,2))+" Comparat+ 2016",fluxscale = 1,em = "[OIII]",filter = filter,style = "g")
	
	#if zHalphalow>0:
	#	schechter_LF(z=zHalphalow,lambdaemitted = lambda_Halpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = r"H$\alpha$ z = "+str(round(zHalphalow,2))+" Sobral+ 2013",fluxscale = 1, em = "Halpha",filter = filter,style = "b")
	
	#if zLymanalphalow>0:
		#the following is also from the separate linearequation code that I tried to put in this one, but it throws back errors every time I use the "third" option, so I have to fix that later
	#	schechter_LF(z=zLymanalphalow,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**44.0057781641604,betaL = 0,phistar0 = 10**(-4.068119141604011),betaphi = 0,param = "first",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalphalow,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filter = filter,style = "y")
 		#schechter_LF(z=zLymanalphalow,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 0,betaL = 0,phistar0 = 0,betaphi = 0,param = "third",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalphalow,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filter = filter,style = "y")
 	

	#high end of FWHM

	#if zOIIhigh>0:
	#	schechter_LF(z=zOIIhigh,lambdaemitted = lambda_OII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "[OII] z = "+str(round(zOIIhigh,2))+" Comparat+ 2016",fluxscale = 1,em = "[OII]",filter = filter,style = "r")
	
	#if zOIIIhigh>0:
	#	schechter_LF(z=zOIIIhigh,lambdaemitted = lambda_OIII, alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "[OIII] z = "+str(round(zOIIIhigh,2))+" Comparat+ 2016",fluxscale = 1,em = "[OIII]",filter = filter,style = "g")
	
	#if zHalphahigh>0:
	#	schechter_LF(z=zHalphahigh,lambdaemitted = lambda_Halpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = r"H$\alpha$ z = "+str(round(zHalphahigh,2))+" Sobral+ 2013",fluxscale = 1, em = "Halpha",filter = filter,style = "b")
	
	#if zLymanalphahigh>0:
		#the following is also from the separate linearequation code that I tried to put in this one, but it throws back errors every time I use the "third" option, so I have to fix that later
	#	schechter_LF(z=zLymanalphahigh,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**44.0057781641604,betaL = 0,phistar0 = 10**(-4.068119141604011),betaphi = 0,param = "first",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalphahigh,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filter = filter,style = "y")
 		#schechter_LF(z=zLymanalphahigh,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 0,betaL = 0,phistar0 = 0,betaphi = 0,param = "third",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalphahigh,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filter = filter,style = "y")
 	

if emline == "allFWHM":

	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm

	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	FWHMlow = lambdaarray[1]
	FWHMhigh = lambdaarray[2]

	#once again, changed below to dictionaries above

	#if filt=="uband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.1
	#	lambdalow = 305.30 #in nm
	#	lambdahigh = 408.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="gband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 27.4
	#	lambdalow = 386.30 #in nm
	#	lambdahigh = 567.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="rband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 27.5
	#	lambdalow = 536.90 #in nm
	#	lambdahigh = 706.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="iband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.8
	#	lambdalow = 675.90 #in nm
	#	lambdahigh = 833.00 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="zband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 26.1
	#	lambdalow = 802.90 #in nm
	#	lambdahigh = 938.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]

	#if filt=="yband":
		#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
		#ABmag = 24.9
	#	lambdalow = 908.30 #in nm
	#	lambdahigh = 1099.60 #in nm
	#	lambdaarray = filter_int(filt = filt)
	#	lambdacenter = lambdaarray[0]
	#	FWHMlow = lambdaarray[1]
	#	FWHMhigh = lambdaarray[2]
	
	lambda_OII = 372.7 #in nm
	lambda_OIII = 500.7 #in nm
	lambda_Halpha = 656.3 #in nm
	lambda_Lymanalpha = 121.6 #in nm

	#finds redshift of emission lines

	#first define array of 100 evenly spaced wavelengths
	lambdaFWHMarray = numpy.linspace(FWHMlow,FWHMhigh,num=100)

	#use these to get an array of redshifts
	zFWHMarrayOII = numpy.zeros(100)
	zFWHMarrayOIII = numpy.zeros(100)
	zFWHMarrayHalpha = numpy.zeros(100)
	zFWHMarrayLymanalpha = numpy.zeros(100)

	#constructs z array from wavelength array within FWHM
	#do this for each emline
	for l in range(len(lambdaFWHMarray)):
		zFWHMarrayOII[l] = (lambdaFWHMarray[l]/lambda_OII)-1
		zFWHMarrayOIII[l] = (lambdaFWHMarray[l]/lambda_OIII)-1
		zFWHMarrayHalpha[l] = (lambdaFWHMarray[l]/lambda_Halpha)-1
		zFWHMarrayLymanalpha[l] = (lambdaFWHMarray[l]/lambda_Lymanalpha)-1

	#don't want negative redshifts
	zFWHMarrayOII = zFWHMarrayOII[numpy.where(zFWHMarrayOII>0)]
	zFWHMarrayOIII = zFWHMarrayOIII[numpy.where(zFWHMarrayOIII>0)]
	zFWHMarrayHalpha = zFWHMarrayHalpha[numpy.where(zFWHMarrayHalpha>0)]
	zFWHMarrayLymanalpha = zFWHMarrayLymanalpha[numpy.where(zFWHMarrayLymanalpha>0)]
	#print("zFWHMarray",zFWHMarray)

 	#NEED TO HAVE AN ARRAY OF Z VALUES TO INTEGRATE PHI
 	#set up empty arrays to print out and eventually integrate comovingphi and arealphi
	comovingphiarrayOII = numpy.zeros(100)
	comovingphiarrayOIII = numpy.zeros(100)
	comovingphiarrayHalpha = numpy.zeros(100)
	comovingphiarrayLymanalpha = numpy.zeros(100)


	#NEED TO CHANGE ALL THE FOR LOOPS TO JUST ARRAY CALCULTIONS - FASTER

	if len(zFWHMarrayOII)>0:

		for l in range(len(zFWHMarrayOII)):
			zOII = zFWHMarrayOII[l]
			comovingphiarrayOII[l] = schechter_LF(z=zOII,lambdaemitted = lambda_OII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "[OII] z = "+str(round(zOII,2))+" Comparat+ 2016",fluxscale = 1,em = "[OII]",filt = filt,style = "r")
		
		#shortens to use only positive values
		zFWHMarrayOII = zFWHMarrayOII[numpy.where(comovingphiarrayOII>0)]
		comovingphiarrayOII = comovingphiarrayOII[numpy.where(comovingphiarrayOII>0)]

		print("OII:")

		
		#first convert zFWHMarrayHalpha to a comoving volume array
		comovingvolarrayOII = numpy.zeros(len(zFWHMarrayOII))

		#then multiply phi by the comoving volume array to get total number of galaxies expected
		#the following finds the comoving volume array to do so
		for r in range(len(zFWHMarrayOII)):
			temparray = cosmo.comoving_volume(zFWHMarrayOII[r]) #units are in Mpc^3
			#change tuple to value
			temparray = temparray.value
			comovingvolarrayOII[r] = temparray #this is in Mpc^3

		#THIS NEEDS TO BE FIXED
		#1 steradian = (180/pi)^2 degrees


		#this is ending up as a weird number
		#integrates comovingphi over whole sky
		total_number_FWHM_OII = numpy.trapz(comovingphiarrayOII,x=comovingvolarrayOII)
		print("total_number_FWHM_OII = ",total_number_FWHM_OII)
		#LSST area
		total_number_FWHM_LSSTOII = total_number_FWHM_OII*18000./42000.		

		print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		print("total_number_FWHM_LSSTOII = ",total_number_FWHM_LSSTOII)

	else:
		zOII = 0


	if len(zFWHMarrayOIII)>0:

		for l in range(len(zFWHMarrayOIII)):
			zOIII = zFWHMarrayOIII[l]
			comovingphiarrayOIII[l] = schechter_LF(z=zOIII,lambdaemitted = lambda_OIII, alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "[OIII] z = "+str(round(zOIII,2))+" Comparat+ 2016",fluxscale = 1,em = "[OIII]",filt = filt,style = "g")
		
		#shortens to use only positive values
		zFWHMarrayOIII = zFWHMarrayOIII[numpy.where(comovingphiarrayOIII>0)]
		comovingphiarrayOIII = comovingphiarrayOIII[numpy.where(comovingphiarrayOIII>0)]

		print("OIII:")

		
		#first convert zFWHMarrayOIII to a comoving volume array
		comovingvolarrayOIII = numpy.zeros(len(zFWHMarrayOIII))

		#then multiply phi by the comoving volume array to get total number of galaxies expected
		#the following finds the comoving volume array to do so
		for r in range(len(zFWHMarrayOIII)):
			temparray = cosmo.comoving_volume(zFWHMarrayOIII[r]) #units are in Mpc^3
			#change tuple to value
			temparray = temparray.value
			comovingvolarrayOIII[r] = temparray #this is in Mpc^3

		#THIS NEEDS TO BE FIXED
		#1 steradian = (180/pi)^2 degrees


		#this is ending up as a weird number
		#integrates comovingphi over whole sky
		total_number_FWHM_OIII = numpy.trapz(comovingphiarrayOIII,x=comovingvolarrayOIII)
		print("total_number_FWHM_OIII = ",total_number_FWHM_OIII)
		#LSST area
		total_number_FWHM_LSSTOIII = total_number_FWHM_OIII*18000./42000.		

		print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		print("total_number_FWHM_LSSTOIII = ",total_number_FWHM_LSSTOIII)

	else:
		zOIII = 0


	if len(zFWHMarrayHalpha)>0:

		for l in range(len(zFWHMarrayHalpha)):
			zHalpha = zFWHMarrayHalpha[l]
			comovingphiarrayHalpha[l] = schechter_LF(z=zHalpha,lambdaemitted = lambda_Halpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = r"H$\alpha$ z = "+str(round(zHalpha,2))+" Sobral+ 2013",fluxscale = 1, em = "Halpha",filt = filt,style = "b")
		
		#shortens to use only positive values
		zFWHMarrayHalpha = zFWHMarrayHalpha[numpy.where(comovingphiarrayHalpha>0)]
		comovingphiarrayHalpha = comovingphiarrayHalpha[numpy.where(comovingphiarrayHalpha>0)]

		print("Halpha:")

		
		#first convert zFWHMarrayHalpha to a comoving volume array
		comovingvolarrayHalpha = numpy.zeros(len(zFWHMarrayHalpha))

		#then multiply phi by the comoving volume array to get total number of galaxies expected
		#the following finds the comoving volume array to do so
		for r in range(len(zFWHMarrayHalpha)):
			temparray = cosmo.comoving_volume(zFWHMarrayHalpha[r]) #units are in Mpc^3
			#change tuple to value
			temparray = temparray.value
			comovingvolarrayHalpha[r] = temparray #this is in Mpc^3

		#THIS NEEDS TO BE FIXED
		#1 steradian = (180/pi)^2 degrees


		#this is ending up as a weird number
		#integrates comovingphi over whole sky
		total_number_FWHM_Halpha = numpy.trapz(comovingphiarrayHalpha,x=comovingvolarrayHalpha)
		print("total_number_FWHM_Halpha = ",total_number_FWHM_Halpha)
		#LSST area
		total_number_FWHM_LSSTHalpha = total_number_FWHM_Halpha*18000./42000.

		print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		print("total_number_FWHM_LSSTHalpha = ",total_number_FWHM_LSSTHalpha)

	else:
		zHalpha = 0


	if len(zFWHMarrayLymanalpha)>0:

		for l in range(len(zFWHMarrayLymanalpha)):
			zLymanalpha = zFWHMarrayLymanalpha[l]
			#comovingphiarrayLymanalpha[l] = schechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**44.0057781641604,betaL = 0,phistar0 = 10**(-4.068119141604011),betaphi = 0,param = "first",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalpha,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filt = filt,style = "y")
			#fixed - the following is also from the separate linearequation code that I tried to put in this one, but it throws back errors every time I use the "third" option, so I have to fix that later
			#NOW THIS WORKS YESSSSS
			comovingphiarrayLymanalpha[l] = schechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 0,betaL = 0,phistar0 = 0,betaphi = 0,param = "third",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalpha,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filt = filt,style = "y")

		#shortens to use only positive values
		zFWHMarrayLymanalpha = zFWHMarrayLymanalpha[numpy.where(comovingphiarrayLymanalpha>0)]
		comovingphiarrayLymanalpha = comovingphiarrayLymanalpha[numpy.where(comovingphiarrayLymanalpha>0)]

		print("Lymanalpha:")

		#ok figure out what I am doing
		#I need the total number of galaxies WITHIN the FWHM
		#1 - integrate comovingphi using trapz and the comoving volume - will get total comoving number
		#2 - get LSST areal/sky density from this - multiply by 18000 square degrees to get final answer
		#or skip areal number density and just multiply the result from the comovingphi integration by 18000/42000 to get the LSST area

		#first convert zFWHMarrayLymanalpha to a comoving volume array
		comovingvolarrayLymanalpha = numpy.zeros(len(zFWHMarrayLymanalpha))

		#then multiply phi by the comoving volume array to get total number of galaxies expected
		#the following finds the comoving volume array to do so
		for r in range(len(zFWHMarrayLymanalpha)):
			temparray = cosmo.comoving_volume(zFWHMarrayLymanalpha[r]) #units are in Mpc^3
			#change tuple to value
			temparray = temparray.value
			comovingvolarrayLymanalpha[r] = temparray #this is in Mpc^3




		#THIS NEEDS TO BE FIXED (updated)
		#1 steradian = (180/pi)^2 degrees - ok check what on earth I am doing here and why
		#the answer is off by a factor of ~10^3 => str~10^3degrees - could it just be that the units are (very) off?  
		#what about the difference between the two methods though?  

		#LSST area of the sky is 18000/42000 of the whole sky

		#this is ending up as a weird number
		#I tried this two ways, and they are both coming out to odd numbers that are nowhere near the value obtained for the median transmission wavelength
		#now they are both the same weird number


		print("comovingphiarrayLymanalpha = ",comovingphiarrayLymanalpha)
		print("comovingvolarrayLymanalpha = ",comovingvolarrayLymanalpha)
		print("zFWHMarrayLymanalpha = ",zFWHMarrayLymanalpha)

		#FIRST TRY

		#integrates comovingphi over whole sky
		total_number_FWHM_Lymanalpha = numpy.trapz(comovingphiarrayLymanalpha,x=comovingvolarrayLymanalpha)
		print("total_number_FWHM_Lymanalpha = ",total_number_FWHM_Lymanalpha)
		#LSST area
		total_number_FWHM_LSSTLymanalpha = total_number_FWHM_Lymanalpha*18000./42000. #rename as sky fraction bc *not* total anymore

		print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		print("total_number_FWHM_LSSTLymanalpha = ",total_number_FWHM_LSSTLymanalpha)


		#SECOND TRY

		#numarray = numpy.zeros(len(zFWHMarrayLymanalpha))
		#comovingvoldiffLymanalpha = numpy.diff(comovingvolarrayLymanalpha)
		#print(len(comovingphiarrayLymanalpha)) #is 100
		#print(len(comovingvolarrayLymanalpha)) #is 100
		#print(len(comovingvoldiffLymanalpha)) #is 99
		#for p in range(len(comovingphiarrayLymanalpha)-1):
		#	numarray[p] = comovingphiarrayLymanalpha[p]*comovingvoldiffLymanalpha[p]
		#delete first number corresponding to index 0 to match voldiff array
		#total_number_FWHM_Lymanalpha = numpy.sum(numarray)

		#print("total_number_FWHM_Lymanalpha = ",total_number_FWHM_Lymanalpha)

		#LSST area
		#total_number_FWHM_LSSTLymanalpha = total_number_FWHM_Lymanalpha*18000./42000.


		#print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		#print("total_number_FWHM_LSSTLymanalpha = ",total_number_FWHM_LSSTLymanalpha)


	else:
		zLymanalpha = 0


	#print("comovingphiarrayOII:",comovingphiarrayOII)
	#print("zFWHMarrayOII:",zFWHMarrayOII)
	#print("comovingphiarrayOIII:",comovingphiarrayOIII)
	#print("zFWHMarrayOIII:",zFWHMarrayOIII)
	#print("comovingphiarrayHalpha:",comovingphiarrayHalpha)
	#print("zFWHMarrayHalpha:",zFWHMarrayHalpha)
	#print("comovingphiarrayLymanalpha:",comovingphiarrayLymanalpha)
	#print("zFWHMarrayLymanalpha:",zFWHMarrayLymanalpha)


	print("the total expected number of galaxies in the LSST area (18000/42000) is:")
	print("total_number_FWHM_LSSTOII = ",total_number_FWHM_LSSTOII)
	print("total_number_FWHM_LSSTOIII = ",total_number_FWHM_LSSTOIII)
	print("total_number_FWHM_LSSTHalpha = ",total_number_FWHM_LSSTHalpha)
	print("total_number_FWHM_LSSTLymanalpha = ",total_number_FWHM_LSSTLymanalpha)