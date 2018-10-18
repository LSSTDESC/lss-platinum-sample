import numpy
from matplotlib.pyplot import *
import scipy
from scipy import integrate
from astropy.cosmology import Planck15 as cosmo


def filter_int(filter):

	if filter=="uband":

		#read in each data file
		print("u_filter has wavelengths 305.30 - 408.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		u_filter = loadtxt('ufilteredit.csv')
		print(u_filter)
		#I shorten this to only the second column
		LSST_filter = u_filter
		LSSTfilter = u_filter[:,1]

	if filter=="gband":
		#read in each data file
		print("g_filter has wavelengths 386.30 - 567.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		g_filter = loadtxt('gfilteredit.csv')
		print(g_filter)
		#I shorten this to only the second column
		LSST_filter = g_filter
		LSSTfilter = g_filter[:,1]

	if filter=="rband":

		#read in each data file
		print("r_filter has wavelengths 536.90 - 706.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput 
		r_filter = loadtxt('rfilteredit.csv')
		print(r_filter)
		#I shorten this to only the second column
		LSST_filter = r_filter
		LSSTfilter = r_filter[:,1]

	if filter=="iband":

		#read in each data file
		print("i_filter has wavelengths 675.90 - 833.00 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		i_filter = loadtxt('ifilteredit.csv')
		print(i_filter)
		#I shorten this to only the second column
		LSST_filter = i_filter
		LSSTfilter = i_filter[:,1]

	if filter=="zband":

		#read in each data file 
		print("z_filter has wavelengths 802.90 - 938.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		z_filter = loadtxt('zfilteredit.csv')
		print(z_filter)
		#I shorten this to only the second column
		LSST_filter = z_filter
		LSSTfilter = z_filter[:,1]

	if filter=="yband":

		#read in each data file
		print("y_filter has wavelengths 908.30 - 1099.60 nm")
		#this is in two columns; the left is wavelength, the right is throughput
		y_filter = loadtxt('yfilteredit.csv')
		print(y_filter)
		#I shorten this to only the second column
		LSST_filter = y_filter
		LSSTfilter = y_filter[:,1]

	#to find the midpoint, integrate over the entire filter, then find what value would give you half the integrated value
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
	test = LSST_filter[centerindex]
	print(type(test))
	lambdacenter = test[0]
	lambdacenter = numpy.float(lambdacenter)
	print("median transmission wavelength of the "+filter+" = ",lambdacenter)

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
	if filter == "iband":
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
	if filter != "iband":
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
	#have to change this from tuple to float - first find the actual value in index 0
	indexL = indexL[0]
	indexR = indexR[0]
	indexL = numpy.float(indexL)
	indexR = numpy.float(indexR)
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
	indexleftFWHM = numpy.float(indexleftFWHM)
	indexrightFWHM = numpy.float(indexrightFWHM)

	#match to lambda
	#LSSTlambda = LSST_filter[:,0]
	lambdaFWHMleft = LSSTlambda[indexleftFWHM]
	lambdaFWHMright = LSSTlambda[indexrightFWHM]
	print(lambdaFWHMleft,lambdaFWHMright)

	#should have output as an array (see LF code)
	print("returns array: lambdacenter,lambdaFWHMleft,lambdaFWHMright")
	lambdaarray = [lambdacenter,lambdaFWHMleft,lambdaFWHMright]

	return lambdaarray


#print(filter_int(filter = "uband"))
#print(filter_int(filter = "gband"))
#print(filter_int(filter = "rband"))
#print(filter_int(filter = "iband"))
print(filter_int(filter = "zband"))
#print(filter_int(filter = "yband"))


#first define array of 100 evenly spaced wavelengths

lambdaFWHMleft = lambdaarray[1]
lambdaFWHMright = lambdaarray[2]
lambdaFWHMarray = numpy.linspace(lambdaFWHMleft,lambdaFWHMright,num=100)
print("lambdaFWHMarray",lambdaFWHMarray)

#use these to get an array of redshifts with which to get the comoving volume, for each of them
#basically getting the thing I found below but for everything across the FWHM





