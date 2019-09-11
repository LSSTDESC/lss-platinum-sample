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
print("u: 26.3, g: 27.5, r: 27.7, i: 27.0, z: 26.2, y: 24.9")
filt = input("Plot Schechter Luminosity Function for which filter?  (sample input: zband)")
print("filter = ", filt)


#THIS PART CONTAINS FUNCTIONS THAT WILL THEN BE USED LATER ON IN THE CODE


def lineqLya(z):

	#print("lineqLya has been called")

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

	if filt=="uband":

		#read in each data file
		#print("u_filter has wavelengths 305.30 - 408.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		u_filter = loadtxt('ufilteredit.csv')
		#print(u_filter)
		#I shorten this to only the second column
		LSST_filter = u_filter
		LSSTfilter = u_filter[:,1]

	if filt=="gband":
		#read in each data file
		#print("g_filter has wavelengths 386.30 - 567.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		g_filter = loadtxt('gfilteredit.csv')
		#print(g_filter)
		#I shorten this to only the second column
		LSST_filter = g_filter
		LSSTfilter = g_filter[:,1]

	if filt=="rband":

		#read in each data file
		#print("r_filter has wavelengths 536.90 - 706.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput 
		r_filter = loadtxt('rfilteredit.csv')
		#print(r_filter)
		#I shorten this to only the second column
		LSST_filter = r_filter
		LSSTfilter = r_filter[:,1]

	if filt=="iband":

		#read in each data file
		#print("i_filter has wavelengths 675.90 - 833.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		i_filter = loadtxt('ifilteredit.csv')
		#print(i_filter)
		#I shorten this to only the second column
		LSST_filter = i_filter
		LSSTfilter = i_filter[:,1]

	if filt=="zband":

		#read in each data file 
		#print("z_filter has wavelengths 802.90 - 938.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		z_filter = loadtxt('zfilteredit.csv')
		#print(z_filter)
		#I shorten this to only the second column
		LSST_filter = z_filter
		LSSTfilter = z_filter[:,1]

	if filt=="yband":

		#read in each data file
		#print("y_filter has wavelengths 908.30 - 1099.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		y_filter = loadtxt('yfilteredit.csv')
		#print(y_filter)
		#I shorten this to only the second column
		LSST_filter = y_filter
		LSSTfilter = y_filter[:,1]

	#to find the midpoint, integrate over the entire filter, then find what value would give you half the integrated value
	#see if you need to fix this because should have x=
	stepint_LSSTfilter = scipy.integrate.cumtrapz(LSSTfilter)
	#print('stepint_LSSTfilter')
	#print(stepint_LSSTfilter)
	len(stepint_LSSTfilter)
	#len gives total length, which is index number + 1
	lastnumber = stepint_LSSTfilter[len(stepint_LSSTfilter)-1]
	midpoint = lastnumber/2

	#finds closest value by finding minimum difference
	difference = abs(midpoint-stepint_LSSTfilter)
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
	centerval = LSST_filter[centerindex]
	#print(type(centerval))
	lambdacenter = centerval[0]
	lambdacenter = numpy.float(lambdacenter)
	#print("median transmission wavelength of the "+filt+" = ",lambdacenter)

	#the "central" wavelength is called the median transmission wavelength of each band

	#want to find the FWHM of each filter

	#finds the index of the maximum value and the corresponding wavelength
	maxval = max(LSSTfilter)
	maxindex = numpy.where(LSSTfilter==maxval)
	#print(type(maxindex))
	#want to change maxindex from tuple to float - just take value
	maxindex = maxindex[0]
	maxindex = numpy.float(maxindex)
	maxindex = numpy.int(maxindex) #attempt at fixing an error
	LSSTlambda = LSST_filter[:,0]
	maxlambda = LSSTlambda[maxindex]
	#print("maxlambda = ",maxlambda)

	#need the two points of the filter whose values give half the maximum value

	#checks first and last fourth of iband filter for half max value because of weird iband shape
	if filt == "iband":
		halfmaxval = maxval/2
		fourth = (len(LSSTfilter))/4
		#print(fourth*4,fourth)
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
		#print(half*2,half)
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
	#print("closestleft = ",closestleft)
	#print("closestright = ",closestright)
	indexL = numpy.where(diffleft==closestleft)
	indexR = numpy.where(diffright==closestright)
	#print("indexL = ", indexL)
	#print("indexR = ", indexR)
	#have to change this from tuple to int - first find the actual value in index 0
	indexL = indexL[0]
	indexR = indexR[0]
	indexL = numpy.int(indexL)
	indexR = numpy.int(indexR)
	#print("indexL = ",indexL)
	#print("indexR = ",indexR)

	#find value in left and right sections of the filter
	valueleft = left[indexL]
	valueright = right[indexR]
	#print(valueleft,valueright)

	#match value to LSSTfilter
	LSSThalfmaxleft = numpy.where(LSSTfilter==valueleft)
	LSSThalfmaxright = numpy.where(LSSTfilter==valueright)
	#print(LSSThalfmaxleft,LSSThalfmaxright)
	#print(type(LSSThalfmaxleft))
	#print(type(LSSThalfmaxright))
	#want to change from tuple to float - just take value
	indexleftFWHM = LSSThalfmaxleft[0]
	indexrightFWHM = LSSThalfmaxright[0]
	indexleftFWHM = numpy.int(indexleftFWHM)
	indexrightFWHM = numpy.int(indexrightFWHM)

	#match to lambda
	#LSSTlambda = LSST_filter[:,0]
	lambdaFWHMleft = LSSTlambda[indexleftFWHM]
	lambdaFWHMright = LSSTlambda[indexrightFWHM]
	#print(lambdaFWHMleft,lambdaFWHMright)

	#should have output as an array (see LF code)
	#print("returns array: lambdacenter,lambdaFWHMleft,lambdaFWHMright")
	lambdaarray = [lambdacenter,lambdaFWHMleft,lambdaFWHMright]

	return lambdaarray


#the following function calculates the luminosity limit for a certain band with a certain detection limit
def lumlim(z,em,filt):

	#print("lumlim has been called")

	#REMINDER TO SELF - z is a function of the emission line wavelength and the filter itself, so it does not have to be adjusted for the FWHM thing

	#print("This function will calculate the luminosity that corresponds to a 5 sigma detection, in erg/s:")
	#print("INFO: 26.2 is AB magnitude for 5 sigma detection limit in the z band") #this was used when I only had the z band
	#first I calculate the flux, then convert to flux density, then find the luminosity limit for the conditions printed above

	#changed if statements (now deleted) to dictionaries below
	#ABmag is the coadded depth for a 5 sigma magnitude limit in this filter
	ABmag_dict = {"uband":26.3,"gband":27.5,"rband":27.7,"iband":27.0,"zband":26.2,"yband":24.9}
	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm
	#these are at the endpoints of the filters, NOT the FWHM

	ABmag = ABmag_dict[filt]
	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]

	#might need the following here:

	if filt=="uband":

		#read in each data file
		#print("u_filter has wavelengths 305.30 - 408.60 nm")
		#this is in two columns; the left is wavelength in nm, the right is throughput
		u_filter = loadtxt('ufilteredit.csv')		#print(u_filter)
		#I shorten this to only the second column
		LSST_filter = u_filter
		LSSTfilter = u_filter[:,1] #transmission
		LSSTwav = u_filter[:,0] #the wavelengths are the first column

	if filt=="gband":
		#read in each data file
		#print("g_filter has wavelengths 386.30 - 567.00 nm")
		#this is in two columns; the left is wavelength in nm, the right is throughput
		g_filter = loadtxt('gfilteredit.csv')
		#print(g_filter)
		#I shorten this to only the second column
		LSST_filter = g_filter
		LSSTfilter = g_filter[:,1] #transmission
		LSSTwav = g_filter[:,0] #the wavelengths are the first column

	if filt=="rband":

		#read in each data file
		#print("r_filter has wavelengths 536.90 - 706.00 nm")
		#this is in two columns; the left is wavelength in nm, the right is throughput 
		r_filter = loadtxt('rfilteredit.csv')
		#print(r_filter)
		#I shorten this to only the second column
		LSST_filter = r_filter
		LSSTfilter = r_filter[:,1] #transmission
		LSSTwav = r_filter[:,0] #the wavelengths are the first column

	if filt=="iband":

		#read in each data file
		#print("i_filter has wavelengths 675.90 - 833.00 nm")
		#this is in two columns; the left is wavelength in nm, the right is throughput
		i_filter = loadtxt('ifilteredit.csv')
		#print(i_filter)
		#I shorten this to only the second column
		LSST_filter = i_filter
		LSSTfilter = i_filter[:,1] #transmission
		LSSTwav = i_filter[:,0] #the wavelengths are the first column

	if filt=="zband":

		#read in each data file 
		#print("z_filter has wavelengths 802.90 - 938.60 nm")
		#this is in two columns; the left is wavelength in nm, the right is throughput
		z_filter = loadtxt('zfilteredit.csv')
		#print(z_filter)
		#I shorten this to only the second column
		LSST_filter = z_filter
		LSSTfilter = z_filter[:,1] #transmission
		LSSTwav = z_filter[:,0] #the wavelengths are the first column

	if filt=="yband":

		#read in each data file
		#print("y_filter has wavelengths 908.30 - 1099.60 nm")
		#this is in two columns; the left is wavelength in nm, the right is throughput
		y_filter = loadtxt('yfilteredit.csv')
		#print(y_filter)
		#I shorten this to only the second column
		LSST_filter = y_filter
		LSSTfilter = y_filter[:,1] #transmission
		LSSTwav = y_filter[:,0] #the wavelengths are the first column


	#the following calculates values I use in and plug into the lumlim function

	#need to redo how I find the flux density to take into account the shape of the transmission curve of each filter
	#I will get a (limiting) flux for this (see notes)


	#LSSTfilter is the transmission array

	#steps to take: FOR EACH FILTER


	#(1)
	#n_photon(1_microJansky) = int{Transmission(lambda)*[1microJansky/(hc/lambda)]*dlambda}
	#hc/lambda is the energy per photon
	#lambdaem is the emitted wavelength
	#do this for the chosen filter (corresponding to the transmission curve)
	#check units - make sure n ends up as number/((cm^2)*s)

	h_cgs = 6.6261*(10**(-27)) #g*cm^2/s
	c_cgs = 2.9979*(10**10) #in cm/s
	lambdaem_dict_cgs = {"[OII]":372.7/(10**7),"[OIII]":500.7/(10**7),"Halpha":656.3/(10**7),"Lymanalpha":121.6/(10**7)} #in cm
	#print("TEST HERE: lambdaem_dict_cgs =",lambdaem_dict_cgs)
	lambdaem_cgs = lambdaem_dict_cgs[em]
	#print("TEST HERE: lambdaem_cgs =",lambdaem_cgs)
	#print("TEST HERE: LSSTwav =",LSSTwav)
	LSSTwav = LSSTwav/(10**7)  #nm to cm
	#print("TEST HERE: LSSTwav =",LSSTwav)

	#print("TEST HERE: LSSTfilter: (this may be the issue?",LSSTfilter)
	n_photon_1microJansky_array = LSSTfilter*(10**(-29))/(h_cgs*LSSTwav) #c_cgs/ #fixed the lambda from an incorrect value to the LSSTwav from the transmission function array
	#print("TEST HERE: n_photon_1microJansky_array",n_photon_1microJansky_array)
	#dlambda = abs(LSSTwav[1]-LSSTwav[0])
	#print("TEST HERE: dlambda",dlambda)
	#n_photon_1microJansky = numpy.trapz(n_photon_1microJansky_array,dx=dlambda)
	n_photon_1microJansky = numpy.trapz(n_photon_1microJansky_array,x=LSSTwav) #more general
	#print("n_photon_1microJansky",n_photon_1microJansky)
	#print("TEST HERE: n_photon_1microJansky = ",n_photon_1microJansky)


	#(2)
	#lambda_emissionline = lambda_restframe*(1+z_emissionline)

	lambda_emissionline = lambdaem_cgs*(1+z)  #in cm
	#print("TEST HERE: z =",z)
	#print("TEST HERE: lambdaem_cgs =",lambdaem_cgs)
	#print("TEST HERE: lambda_emissionline =",lambda_emissionline)


	#(3)

	#need LSSTfilter[lambda_emissionline]
	#will find this using minimum difference between lambda_emissionline and whatever wavelengths are in the array

	#LSSTfilter is the transmission
	#LSSTwav is the wavelength array
	#ultimately need a value in LSSTfilter using an index from LSSTwav

	#finds the index that most closely matches the wavelength that is redshifted

	diff_lambda_array = abs(LSSTwav-lambda_emissionline) #SHOULD BE IN SAME UNITS - they are both in centimeters in this section of the code
	#do without absolute value
	#print("TEST HERE: LSSTwav =",LSSTwav)
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
	#LSSTfilter[mindiff_lambda_index]
	#print("TEST HERE: LSSTfilter[mindiff_lambda_index] =",LSSTfilter[mindiff_lambda_index])


	#HERE STOP TEST


	#(4)
	#flux_limiting(emissionline) = [F_nu,lim(filt)*n_photon(1_microJansky)/T_filt(lambda_emissionline)]*(hc/lambda_emissionline)
	#F_nu,lim(filt) is the variable I had previously named fluxdens - it is the 5 sigma AB magnitude limit that I found for each LSST filter

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

	#LSSTfilter is the transmission
	#print("f_n_obj / [1_microJansky] = n_photon_object / n_photon_1_microJansky")
	flux_limit = (fluxdens*n_photon_1microJansky/(LSSTfilter[int(mindiff_lambda_index)]*(10**(-29))))*(h_cgs*c_cgs/lambda_emissionline)
	#print("flux_limit",flux_limit)

	#print("FLUX LIMIT IS: for z=",z)
	#print(" and em=",em)
	#print(" is:",flux_limit)


	#print("ADAM'S TEST:")
	#print("hc/lambda should be ~1.6e-11 ergs:",h_cgs*c_cgs/lambdaem_cgs)
	#print("T_EL should be of order 0.2:",LSSTfilter[mindiff_lambda_index])
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

	if param == "first":

		Lstar = fluxscale*Lstar0*((1+z)**betaL)
		#print("Lstar = ", Lstar)

		phistar = phistar0*((1+z)**betaphi)
		#print("phistar = ", phistar)


	if param == "second":

		Lstar = 10**(0.45*z + Lstar0)
		#print("Lstar = ", Lstar)

		phistar = 10**(-0.38*(z**2) + z + phistar0)
		#print("phistar = ", phistar)


	if param == "third":

		answersLya = lineqLya(z=z)
		#print(type(answersLya))

		Lstar = answersLya[0] #this is linearly parametrized from Ciardullo+ 2012
		#print("Lstar = ",Lstar)
		#print(type(Lstar))

		phistar = answersLya[1] #this is linearly parametrized from Ciardullo+ 2012
		#print("phistar = ",phistar)
		#print(type(phistar))



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
	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm

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
	lumlim_emline = lumlim(z = z,em = em,filt = filt)


	#print("luminosity limit for the redshifted emission line:",lumlim_emline)


	#beginning of edits and temporary commenting out of plot stuff

	#finds luminosity limit in center of each band
	center = lumlim(z = filtercenter,em = em,filt = filt)

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
	FWHMlowpoint = lumlim(z = FWHMendlow,em = em,filt = filt)
	FWHMhighpoint = lumlim(z = FWHMendhigh,em = em,filt = filt)

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



#THIS NEXT PART CONTAINS IF STATEMENTS FOR EACH EMLINE


#functionality of code has changed - previously had separate statements (now..hesitantly deleted) for each emission line as I was beginning, but now I combined them all into allFWHM

#the following seemed important, so I kept them for reference

#	print("[OII] values with errors given in the Comparat et al 2016 paper are:")
#	print("log10Lstar =   41.1  (+0.02 -0.03)")
#	print("betaL =         2.33 (+0.14 -0.16)")
#	print("log10phistar = -2.4  (+0.03 -0.03)")
#	print("betaphi =      -0.73 (+0.25 -0.29)")
#	print("alpha =        -1.46 (+0.06 -0.05)")

#	print("[OIII] values with errors given in the Comparat et al 2016 paper are:")
#	print("log10Lstar =   41.42 (+0.07 -0.09)")
#	print("betaL =         3.91 (+0.32 -0.4)")
#	print("log10phistar = -3.41 (+0.08 -0.01)")
#	print("betaphi =      -0.76 (+0.39 -0.49)")
#	print("alpha =        -1.83 (+0.1 -0.08)")

#	print("z = 0.4 and alpha = -1.6 are plotted, ")
#	print("Halpha values with errors given in the Sobral et al 2013 paper are:")
#	print("z = 0.40 +_ 0.01")
#	print("alpha =       -1.75 (+0.12 -0.08)")
#	print("Lstar_exp =   41.95 (+0.47 -0.12)")
#	print("phistar_exp = -3.12 (+0.10 -0.34)")

	#for the first redshift below
#	print("z = 3.113 and alpha = -1.65 are plotted, ")
#	print("Lymanalpha values with errors given in the Ciardullo et al 2012 paper are:")
#	print("Lstar_exp =   42.76 (+0.10 -0.10)")
#	print("phistar_exp = -3.17 (+0.05 -0.05)")

	#for the second redshift below
#	print("z = 2.063 and alpha = -1.65 are plotted, ")
#	print("Lymanalpha values with errors given in the Ciardullo et al 2012 paper are:")
#	print("Lstar_exp =   42.33 (+0.12 -0.12)")
#	print("phistar_exp = -2.86 (+0.05 -0.05)")


#might not be able to use this while editing the Schechter_LF routine for the allFWHM option, bc integrating in different order
#I am keeing it for reference
#Anakin don't try it
if emline == "all":

	#if statements deleted again changed to dictionaries
	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm

	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	FWHMlow = lambdaarray[1]
	FWHMhigh = lambdaarray[2]

	#this was used to print out data and plots for the LSST-DESC conference poster

	print("This option finds comoving and areal number densities for:")
	print("[OII] 3726/3729 and [OIII] 5007 from Comparat et al 2016")
	print("Halpha 6563 from Sobral et al 2013")

	lambda_OII = 372.7
	lambda_OIII = 500.7
	lambda_Halpha = 656.3
	lambda_Lymanalpha = 121.6

	#finds redshift of emission line at center of z band -> use to find LF
	zOII = (lambdacenter/lambda_OII)-1
	zOIII = (lambdacenter/lambda_OIII)-1
	zHalpha = (lambdacenter/lambda_Halpha)-1
	zLymanalpha = (lambdacenter/lambda_Lymanalpha)-1

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


if emline == "allFWHM":

	#changed if statements to dictionaries and deleted them
	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm

	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	FWHMlow = lambdaarray[1]
	FWHMhigh = lambdaarray[2]
	
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


		#adds legend

		firstz = zFWHMarrayOII[0]
		lastz = zFWHMarrayOII[-1]
		zpaper = "[OII] z = "+str(round(firstz,2))+"-"+str(round(lastz,2))+" Comparat+ 2016"

		figure(1)
		plot(0,0,"r",label=zpaper)
		legend(loc="upper right")

		figure(2)
		plot(0,0,"r",label=zpaper)
		legend(loc="upper right")

		show()


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


		#adds legend

		firstz = zFWHMarrayOIII[0]
		lastz = zFWHMarrayOIII[-1]
		zpaper = "[OIII] z = "+str(round(firstz,2))+"-"+str(round(lastz,2))+" Comparat+ 2016"

		figure(1)
		plot(0,0,"g",label=zpaper)
		legend(loc="upper right")

		figure(2)
		plot(0,0,"g",label=zpaper)
		legend(loc="upper right")

		show()


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


		#adds legend

		firstz = zFWHMarrayHalpha[0]
		lastz = zFWHMarrayHalpha[-1]
		zpaper = r"H$\alpha$ z = "+str(round(firstz,2))+"-"+str(round(lastz,2))+" Sobral+ 2013"

		figure(1)
		plot(0,0,"b",label=zpaper)
		legend(loc="upper right")

		figure(2)
		plot(0,0,"b",label=zpaper)
		legend(loc="upper right")

		show()


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
			#OKAY THIS NUMBER IS small?
			comovingphiarrayLymanalpha[l] = schechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 0,betaL = 0,phistar0 = 0,betaphi = 0,param = "third",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalpha,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filt = filt,style = "y")
			

		#shortens to use only positive values
		zFWHMarrayLymanalpha = zFWHMarrayLymanalpha[numpy.where(comovingphiarrayLymanalpha>0)]
		comovingphiarrayLymanalpha = comovingphiarrayLymanalpha[numpy.where(comovingphiarrayLymanalpha>0)]


		#adds legend

		firstz = zFWHMarrayLymanalpha[0]
		lastz = zFWHMarrayLymanalpha[-1]
		zpaper = r"Ly$\alpha$ z = "+str(round(firstz,2))+"-"+str(round(lastz,2))+" Ciardullo+ 2012"

		figure(1)
		plot(0,0,"y",label=zpaper)
		legend(loc="upper right")

		figure(2)
		plot(0,0,"y",label=zpaper)
		legend(loc="upper right")

		show()


		print("Lymanalpha:")

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


		#LSST area of the sky is 18000/42000 of the whole sky

		#this is ending up as a weird number


		print("comovingphiarrayLymanalpha = ",comovingphiarrayLymanalpha)
		print("comovingvolarrayLymanalpha = ",comovingvolarrayLymanalpha)
		print("zFWHMarrayLymanalpha = ",zFWHMarrayLymanalpha)


		#integrates comovingphi over whole sky
		total_number_FWHM_Lymanalpha = numpy.trapz(comovingphiarrayLymanalpha,x=comovingvolarrayLymanalpha)
		print("total_number_FWHM_Lymanalpha = ",total_number_FWHM_Lymanalpha)
		#LSST area
		total_number_FWHM_LSSTLymanalpha = total_number_FWHM_Lymanalpha*18000./42000. #rename as sky fraction bc *not* total anymore


		print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		print("total_number_FWHM_LSSTLymanalpha = ",total_number_FWHM_LSSTLymanalpha)


	else:
		zLymanalpha = 0


	print("the total expected number of galaxies in the LSST area (18000/42000) is:")
	if len(zFWHMarrayOII)>0:
		print("total_number_FWHM_LSSTOII = ",total_number_FWHM_LSSTOII)
	if len(zFWHMarrayOIII)>0:
		print("total_number_FWHM_LSSTOIII = ",total_number_FWHM_LSSTOIII)
	if len(zFWHMarrayHalpha)>0:
		print("total_number_FWHM_LSSTHalpha = ",total_number_FWHM_LSSTHalpha)
	if len(zFWHMarrayLymanalpha)>0:
		print("total_number_FWHM_LSSTLymanalpha = ",total_number_FWHM_LSSTLymanalpha)


if emline == "allends":  #WHAT IS THE POINT OF THIS?  is it the same as allFWHM?????????????????????

	#PAY ATTENTION TO THE FOLLOWING COMMENT / PRINT STATEMENT
	print("NOTE: PERHAPS TERRIBLY, I HAVE REDEFINED THE FIRST VALUES LABELED FWHM TO ACTUALLY THE FILTER ENDPOINTS..WHICH I NEED TO FIX LATER.  I JUST WANT TO CHECK STUFF")

	#now check the fake lf and see if good and yayyyyyyyyyy

	#changed if statements to dictionaries and deleted them
	lambdalow_dict = {"uband":305.30,"gband":386.30,"rband":536.90,"iband":675.90,"zband":802.90,"yband":908.30} #in nm
	lambdahigh_dict = {"uband":408.60,"gband":567.00,"rband":706.00,"iband":833.00,"zband":938.60,"yband":1099.60} #in nm

	lambdalow = lambdalow_dict[filt]
	lambdahigh = lambdahigh_dict[filt]
	lambdaarray = filter_int(filt = filt)
	lambdacenter = lambdaarray[0]
	FWHMlow = lambdaarray[1]  #not here
	FWHMhigh = lambdaarray[2]  #not here
	FWHMlow = lambdalow
	FWHMhigh = lambdahigh
	


	#plots the filter transmission curve

	figure(3)

	#read in each data file
	print("u_filter has wavelengths 305.30 - 408.60 nm")
	print("g_filter has wavelengths 386.30 - 567.00 nm")
	print("r_filter has wavelengths 536.90 - 706.00 nm")
	print("i_filter has wavelengths 675.90 - 833.00 nm")
	print("z_filter has wavelengths 802.90 - 938.60 nm")
	print("y_filter has wavelengths 908.30 - 1099.60 nm")

	#need two things - wavelength array and throughput array
	#this is in two columns; the left is wavelength in nm, the right is throughput

	u_filter = loadtxt('ufilteredit.csv')
	u_throughput = u_filter[:,1] #transmission
	u_wav = u_filter[:,0] #the wavelengths are the first column
	
	g_filter = loadtxt('gfilteredit.csv')
	g_throughput = g_filter[:,1] #transmission
	g_wav = g_filter[:,0] #the wavelengths are the first column
	
	r_filter = loadtxt('rfilteredit.csv')
	r_throughput = r_filter[:,1] #transmission
	r_wav = r_filter[:,0] #the wavelengths are the first column
	
	i_filter = loadtxt('ifilteredit.csv')
	i_throughput = i_filter[:,1] #transmission
	i_wav = i_filter[:,0] #the wavelengths are the first column
	
	z_filter = loadtxt('zfilteredit.csv')
	z_throughput = z_filter[:,1] #transmission
	z_wav = z_filter[:,0] #the wavelengths are the first column
	
	y_filter = loadtxt('yfilteredit.csv')
	y_throughput = y_filter[:,1] #transmission
	y_wav = y_filter[:,0] #the wavelengths are the first column

	plot(u_wav,u_throughput,"m")
	plot(g_wav,g_throughput,"b")
	plot(r_wav,r_throughput,"c")
	plot(i_wav,i_throughput,"g")
	plot(z_wav,z_throughput,"y")
	plot(y_wav,y_throughput,"r")

	xlim(300,1100)
	ylim(0,1.0)
	ylabel("Throughput")
	xlabel("Wavelengths (nm)")
	title("LSST Transmission Functions, ugrizy bands")

	show()


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


	#testing individual parts of the LFs

	# zFWHMarrayOII = numpy.zeros(2) 
	# zFWHMarrayOIII = numpy.zeros(2)
	# zFWHMarrayHalpha = numpy.zeros(2)
	# zFWHMarrayLymanalpha = numpy.zeros(2)


	# zFWHMarrayOII[0] = 0.5#(318.2125/lambda_OII)-1
	# zFWHMarrayOII[1] = 1#(344.0375/lambda_OII)-1
	# # #zFWHMarrayOII[2] = (369.8625/lambda_OII)-1
	# # #zFWHMarrayOII[3] = (395.6875/lambda_OII)-1

	# zFWHMarrayOIII[0] = 0.1#(318.2125/lambda_OIII)-1
	# zFWHMarrayOIII[1] = 1#(344.0375/lambda_OIII)-1
	# # #zFWHMarrayOIII[2] = (369.8625/lambda_OIII)-1
	# # #zFWHMarrayOIII[3] = (395.6875/lambda_OIII)-1

	# zFWHMarrayHalpha[0] = 0.1#(318.2125/lambda_Halpha)-1
	# zFWHMarrayHalpha[1] = 1#(344.0375/lambda_Halpha)-1
	# # #zFWHMarrayHalpha[2] = (369.8625/lambda_Halpha)-1
	# # #zFWHMarrayHalpha[3] = (395.6875/lambda_Halpha)-1

	# zFWHMarrayLymanalpha[0] = 0.1#(318.2125/lambda_Lymanalpha)-1
	# zFWHMarrayLymanalpha[1] = 1#(344.0375/lambda_Lymanalpha)-1
	# # #zFWHMarrayLymanalpha[2] = (369.8625/lambda_Lymanalpha)-1
	# # #zFWHMarrayLymanalpha[3] = (395.6875/lambda_Lymanalpha)-1

	# # #print("zFWHMarrayOII:",zFWHMarrayOII)
	# # #print("zFWHMarrayOIII:",zFWHMarrayOIII)
	# # #print("zFWHMarrayHalpha:",zFWHMarrayHalpha)
	# # #print("zFWHMarrayLymanalpha:",zFWHMarrayLymanalpha)

	# #don't want negative redshifts
	# zFWHMarrayOII = zFWHMarrayOII[numpy.where(zFWHMarrayOII>0)]
	# zFWHMarrayOIII = zFWHMarrayOIII[numpy.where(zFWHMarrayOIII>0)]
	# zFWHMarrayHalpha = zFWHMarrayHalpha[numpy.where(zFWHMarrayHalpha>0)]
	# zFWHMarrayLymanalpha = zFWHMarrayLymanalpha[numpy.where(zFWHMarrayLymanalpha>0)]

	#print("zFWHMarrayOII:",zFWHMarrayOII)
	#print("zFWHMarrayOIII:",zFWHMarrayOIII)
	#print("zFWHMarrayHalpha:",zFWHMarrayHalpha)
	#print("zFWHMarrayLymanalpha:",zFWHMarrayLymanalpha)


	#zFWHMarrayOII[0] = 1.47
	#zFWHMarrayOII[1] = 1.525
	#zFWHMarrayOII[2] = 1.585
	#zFWHMarrayOII[3] = 1.6425
	#zFWHMarrayOII[4] = 1.7

	#zFWHMarrayOIII[0] = 0.84
	#zFWHMarrayOIII[1] = 0.8825
	#zFWHMarrayOIII[2] = 0.925
	#zFWHMarrayOIII[3] = 0.9675
	#zFWHMarrayOIII[4] = 1.01

	#zFWHMarrayHalpha[0] = 0.4
	#zFWHMarrayHalpha[1] = 0.4325
	#zFWHMarrayHalpha[2] = 0.465
	#zFWHMarrayHalpha[3] = 0.4975
	#zFWHMarrayHalpha[4] = 0.53

	#zFWHMarrayLymanalpha[0] = 6.57
	#zFWHMarrayLymanalpha[1] = 6.745
	#zFWHMarrayLymanalpha[2] = 6.92
	#zFWHMarrayLymanalpha[3] = 7.095
	#zFWHMarrayLymanalpha[4] = 7.27


	# comovingphiarrayOII = numpy.zeros(2)
	# comovingphiarrayOIII = numpy.zeros(2)
	# comovingphiarrayHalpha = numpy.zeros(2)
	# comovingphiarrayLymanalpha = numpy.zeros(2)


	#NEED TO CHANGE ALL THE FOR LOOPS TO JUST ARRAY CALCULTIONS - FASTER

	print("start [OII]")

	if len(zFWHMarrayOII)>0:

		for l in range(len(zFWHMarrayOII)):
			zOII = zFWHMarrayOII[l]
			comovingphiarrayOII[l] = schechter_LF(z=zOII,lambdaemitted = lambda_OII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = "first",zpaper = "[OII] z = "+str(round(zOII,2))+" Comparat+ 2016",fluxscale = 1,em = "[OII]",filt = filt,style = "r")
		
		#shortens to use only positive values
		zFWHMarrayOII = zFWHMarrayOII[numpy.where(comovingphiarrayOII>0)]
		comovingphiarrayOII = comovingphiarrayOII[numpy.where(comovingphiarrayOII>0)]


		#adds legend

		firstz = zFWHMarrayOII[0]
		lastz = zFWHMarrayOII[-1]
		zpaper = "[OII] z = "+str(round(firstz,2))+"-"+str(round(lastz,2))+" Comparat+ 2016"

		figure(1)
		plot(0,0,"r",label=zpaper)
		legend(loc="upper right")

		figure(2)
		plot(0,0,"r",label=zpaper)
		legend(loc="upper right")

		show()


		#print("OII:")

		
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
		#print("total_number_FWHM_OII = ",total_number_FWHM_OII)
		#LSST area
		total_number_FWHM_LSSTOII = total_number_FWHM_OII*18000./42000.		

		#print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		#print("total_number_FWHM_LSSTOII = ",total_number_FWHM_LSSTOII)


	else:
		zOII = 0


	print("end [OII]")

	print("start [OIII]")

	if len(zFWHMarrayOIII)>0:

		for l in range(len(zFWHMarrayOIII)):
			zOIII = zFWHMarrayOIII[l]
			comovingphiarrayOIII[l] = schechter_LF(z=zOIII,lambdaemitted = lambda_OIII, alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = "first",zpaper = "[OIII] z = "+str(round(zOIII,2))+" Comparat+ 2016",fluxscale = 1,em = "[OIII]",filt = filt,style = "g")
		
		#shortens to use only positive values
		zFWHMarrayOIII = zFWHMarrayOIII[numpy.where(comovingphiarrayOIII>0)]
		comovingphiarrayOIII = comovingphiarrayOIII[numpy.where(comovingphiarrayOIII>0)]


		#adds legend

		firstz = zFWHMarrayOIII[0]
		lastz = zFWHMarrayOIII[-1]
		zpaper = "[OIII] z = "+str(round(firstz,2))+"-"+str(round(lastz,2))+" Comparat+ 2016"

		figure(1)
		plot(0,0,"g",label=zpaper)
		legend(loc="upper right")

		figure(2)
		plot(0,0,"g",label=zpaper)
		legend(loc="upper right")

		show()


		#print("OIII:")

		
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
		#print("total_number_FWHM_OIII = ",total_number_FWHM_OIII)
		#LSST area
		total_number_FWHM_LSSTOIII = total_number_FWHM_OIII*18000./42000.		

		#print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		#print("total_number_FWHM_LSSTOIII = ",total_number_FWHM_LSSTOIII)

	else:
		zOIII = 0


	print("end [OIII]")

	print("start Halpha")


	if len(zFWHMarrayHalpha)>0:

		for l in range(len(zFWHMarrayHalpha)):
			zHalpha = zFWHMarrayHalpha[l]
			comovingphiarrayHalpha[l] = schechter_LF(z=zHalpha,lambdaemitted = lambda_Halpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = "second",zpaper = r"H$\alpha$ z = "+str(round(zHalpha,2))+" Sobral+ 2013",fluxscale = 1, em = "Halpha",filt = filt,style = "b")
		
		#shortens to use only positive values
		zFWHMarrayHalpha = zFWHMarrayHalpha[numpy.where(comovingphiarrayHalpha>0)]
		comovingphiarrayHalpha = comovingphiarrayHalpha[numpy.where(comovingphiarrayHalpha>0)]


		#adds legend

		firstz = zFWHMarrayHalpha[0]
		lastz = zFWHMarrayHalpha[-1]
		zpaper = r"H$\alpha$ z = "+str(round(firstz,2))+"-"+str(round(lastz,2))+" Sobral+ 2013"

		figure(1)
		plot(0,0,"b",label=zpaper)
		legend(loc="upper right")

		figure(2)
		plot(0,0,"b",label=zpaper)
		legend(loc="upper right")

		show()


		#print("Halpha:")

		
		#first convert zFWHMarrayHalpha to a comoving volume array
		comovingvolarrayHalpha = numpy.zeros(len(zFWHMarrayHalpha))

		#can input array into comoving vol calculation - note for future

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
		#print("total_number_FWHM_Halpha = ",total_number_FWHM_Halpha)
		#LSST area
		total_number_FWHM_LSSTHalpha = total_number_FWHM_Halpha*18000./42000.

		#print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		#print("total_number_FWHM_LSSTHalpha = ",total_number_FWHM_LSSTHalpha)

	else:
		zHalpha = 0


	print("end Halpha")

	print("start Lymanalpha")


	if len(zFWHMarrayLymanalpha)>0:

		for l in range(len(zFWHMarrayLymanalpha)):
			zLymanalpha = zFWHMarrayLymanalpha[l]
			#OKAY THIS NUMBER IS small?
			comovingphiarrayLymanalpha[l] = schechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 0,betaL = 0,phistar0 = 0,betaphi = 0,param = "third",zpaper = r"Ly$\alpha$ z = "+str(round(zLymanalpha,2))+" Ciardullo+ 2012",fluxscale = 1,em = "Lymanalpha",filt = filt,style = "y")
			

		#shortens to use only positive values
		zFWHMarrayLymanalpha = zFWHMarrayLymanalpha[numpy.where(comovingphiarrayLymanalpha>0)]
		comovingphiarrayLymanalpha = comovingphiarrayLymanalpha[numpy.where(comovingphiarrayLymanalpha>0)]


		#adds legend

		firstz = zFWHMarrayLymanalpha[0]
		lastz = zFWHMarrayLymanalpha[-1]
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


		#LSST area of the sky is 18000/42000 of the whole sky

		#this is ending up as a weird number


		#print("comovingphiarrayLymanalpha = ",comovingphiarrayLymanalpha)
		#print("comovingvolarrayLymanalpha = ",comovingvolarrayLymanalpha)
		#print("zFWHMarrayLymanalpha = ",zFWHMarrayLymanalpha)


		#integrates comovingphi over whole sky
		total_number_FWHM_Lymanalpha = numpy.trapz(comovingphiarrayLymanalpha,x=comovingvolarrayLymanalpha)
		#print("total_number_FWHM_Lymanalpha = ",total_number_FWHM_Lymanalpha)
		#LSST area
		total_number_FWHM_LSSTLymanalpha = total_number_FWHM_Lymanalpha*18000./42000. #rename as sky fraction bc *not* total anymore


		#print("the total expected number of galaxies in the LSST area (18000/42000) is:")
		#print("total_number_FWHM_LSSTLymanalpha = ",total_number_FWHM_LSSTLymanalpha)


	else:
		zLymanalpha = 0


	print("end Lymanalpha")


	print("the total expected number of galaxies in the LSST area (18000/42000) is:")

	if len(zFWHMarrayOII)>0:
		print("total_number_FWHM_LSSTOII = ",total_number_FWHM_LSSTOII)
		arealphiOII = total_number_FWHM_OII/(4*numpy.pi)
		print("arealphiOII =",arealphiOII,"steradian^-1")

	if len(zFWHMarrayOIII)>0:
		print("total_number_FWHM_LSSTOIII = ",total_number_FWHM_LSSTOIII)
		arealphiOIII = total_number_FWHM_OIII/(4*numpy.pi)
		print("arealphiOIII =",arealphiOIII,"steradian^-1")

	if len(zFWHMarrayHalpha)>0:
		print("total_number_FWHM_LSSTHalpha = ",total_number_FWHM_LSSTHalpha)
		arealphiHalpha = total_number_FWHM_Halpha/(4*numpy.pi)
		print("arealphiHalpha =",arealphiHalpha,"steradian^-1")

	if len(zFWHMarrayLymanalpha)>0:
		print("total_number_FWHM_LSSTLymanalpha = ",total_number_FWHM_LSSTLymanalpha)
		arealphiLymanalpha = total_number_FWHM_Lymanalpha/(4*numpy.pi)
		print("arealphiLymanalpha =",arealphiLymanalpha,"steradian^-1")


#print("These are the endpoints, labeled incorrectly:")
#print("FWHMlow = ",FWHMlow)
#print("FWHMhigh = ",FWHMhigh)


show()