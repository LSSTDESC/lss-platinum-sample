{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy\n",
    "from matplotlib.pyplot import *\n",
    "%matplotlib inline\n",
    "import scipy\n",
    "from scipy import integrate\n",
    "from astropy.cosmology import Planck15 as cosmo\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [],
   "source": [
    "#This code plots Schechter Luminosity Functions for different emission lines in different papers.\n",
    "#Comparat et al 2016 plots [OII] 3726/3729, Hbeta 4861, and [OIII] 5007\n",
    "#[OIII] 5007 is always 3 times stronger than [OIII] 4959\n",
    "#Sobral et al 2013 plots Halpha 6563\n",
    "#Ciardullo et al 2012 plots Lymanalpha 1216\n",
    "#units are in Angstroms\n",
    "#[OII] 3726/3729 unresolved doublet, Hbeta 4861, [OIII] 4959/5007, Halpha 6563, Lymanalpha 1216\n",
    "#the code will be able to use any filter\n",
    "#The ugrizy filter options with corresponding final coadded depths (5 sigma) are as follows:\n",
    "#u: 26.1, g: 27.4, r: 27.5, i: 26.8, z: 26.1, y: 24.9"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Plot Schechter Luminosity Function for which emission line?  (Halpha,Lymanalpha)Lymanalpha\n",
      "emline =  Lymanalpha\n",
      "Plot Schechter Luminosity Function for which filt?  (uband,gband,rband,iband,zband,yband)zband\n",
      "filt =  zband\n",
      "This following plot will contain a Lymanalpha LF from Ciardullo et al 2012.  \n",
      "z = 3.113 and alpha = -1.65 are plotted, \n",
      "the values with errors given in the Ciardullo et al 2012 paper are:\n",
      "Lstar_exp =   42.76 (+0.10 -0.10)\n",
      "phistar_exp = -3.17 (+0.05 -0.05)\n",
      "z = 2.063 and alpha = -1.65 are plotted, \n",
      "the values with errors given in the Ciardullo et al 2012 paper are:\n",
      "Lstar_exp =   42.33 (+0.12 -0.12)\n",
      "phistar_exp = -2.86 (+0.05 -0.05)\n",
      "alpha =  -1.65\n",
      "z =  3.113\n",
      "Lstar =  5.754399373371543e+42\n",
      "phistar =  0.0006760829753919819\n",
      "Note: Each paper uses a slightly different convention; I decided to consolidate them with the following:\n",
      "using the LF equation with the alpha+1 exponent as follows: \n",
      "phi = phistar*((L/Lstar)**(alpha+1))*(numpy.e**(-L/Lstar))\n"
     ]
    },
    {
     "ename": "IndexError",
     "evalue": "boolean index did not match indexed array along dimension 0; dimension is 2500 but corresponding boolean dimension is 1563",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mIndexError\u001b[0m                                Traceback (most recent call last)",
      "\u001b[0;32m<ipython-input-14-b3ad9fa413f3>\u001b[0m in \u001b[0;36m<module>\u001b[0;34m()\u001b[0m\n\u001b[1;32m    500\u001b[0m         \u001b[0mprint\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m\"phistar_exp = -2.86 (+0.05 -0.05)\"\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    501\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 502\u001b[0;31m         \u001b[0mschechter_LF\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mz\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0;36m3.113\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mlambdaemitted\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mlambda_Lymanalpha\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0malpha\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m-\u001b[0m\u001b[0;36m1.65\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mLstar0\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m10\u001b[0m\u001b[0;34m**\u001b[0m\u001b[0;36m42.76\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mbetaL\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m0\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mphistar0\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m10\u001b[0m\u001b[0;34m**\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m-\u001b[0m\u001b[0;36m3.17\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mbetaphi\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m0\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mparam\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m\"first\"\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mzpaper\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34mr\"Ly$\\alpha$ z = 3.113 Ciardullo+ 2012\"\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mfluxscale\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m1\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mem\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m\"Lymanalpha\"\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mfilt\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mfilt\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mstyle\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m\"c\"\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    503\u001b[0m         \u001b[0mschechter_LF\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mz\u001b[0m\u001b[0;34m=\u001b[0m\u001b[0;36m2.063\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mlambdaemitted\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mlambda_Lymanalpha\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0malpha\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m-\u001b[0m\u001b[0;36m1.65\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mLstar0\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m10\u001b[0m\u001b[0;34m**\u001b[0m\u001b[0;36m42.33\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mbetaL\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m0\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mphistar0\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m10\u001b[0m\u001b[0;34m**\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m-\u001b[0m\u001b[0;36m2.86\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mbetaphi\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m0\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mparam\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m\"first\"\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mzpaper\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34mr\"Ly$\\alpha$ z = 2.063 Ciardullo+ 2012\"\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mfluxscale\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;36m1\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mem\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m\"Lymanalpha\"\u001b[0m\u001b[0;34m,\u001b[0m\u001b[0mfilt\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mfilt\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0mstyle\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0;34m\"m\"\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    504\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;32m<ipython-input-14-b3ad9fa413f3>\u001b[0m in \u001b[0;36mschechter_LF\u001b[0;34m(z, lambdaemitted, alpha, Lstar0, betaL, phistar0, betaphi, zpaper, param, fluxscale, em, filt, style)\u001b[0m\n\u001b[1;32m    296\u001b[0m         \u001b[0mL\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mL\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mphi\u001b[0m\u001b[0;34m!=\u001b[0m\u001b[0;36m0\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    297\u001b[0m         \u001b[0mphi\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mphi\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mphi\u001b[0m\u001b[0;34m!=\u001b[0m\u001b[0;36m0\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m--> 298\u001b[0;31m         \u001b[0mtestarray\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mtest\u001b[0m\u001b[0;34m[\u001b[0m\u001b[0mphi\u001b[0m\u001b[0;34m!=\u001b[0m\u001b[0;36m0\u001b[0m\u001b[0;34m]\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m    299\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m    300\u001b[0m         \u001b[0mlog10L\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mnumpy\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mlog10\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mL\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;31mIndexError\u001b[0m: boolean index did not match indexed array along dimension 0; dimension is 2500 but corresponding boolean dimension is 1563"
     ]
    }
   ],
   "source": [
    "emline = input(\"Plot Schechter Luminosity Function for which emission line?  (Halpha,Lymanalpha)\")\n",
    "print(\"emline = \", emline)\n",
    "filt = input(\"Plot Schechter Luminosity Function for which filt?  (uband,gband,rband,iband,zband,yband)\")\n",
    "print(\"filt = \", filt)\n",
    "\n",
    "#lineqLya - interpolates/extrapolates values for the Lyman alpha line\n",
    "\n",
    "#filt_int - finds the central wavelength, actually called the median transmission wavelength, of each filt\n",
    "\n",
    "#lumlim -  calculates the luminosity limit for a each band with a 5 sigma detection limit\n",
    "\n",
    "#schechter_LF - which basically (among calling the other functions to get useful values) plots the Schechter luminosity function for\n",
    "\n",
    "#an Extreme Emission Line galaxy sample and then integrates it to get a number density for specific redshifts (then plots both)\n",
    "\n",
    "\n",
    "def lineqLya(z):\n",
    "\n",
    "\t#this is used for the Lymanalpha line\n",
    "\t#finds linear equation from two points in Ciardullo et al. 2012 paper and interpolates/extraoplates to find Lstar and phistar for different z values\n",
    "\n",
    "\tz1 = 3.113\n",
    "\tz2 = 2.063\n",
    "\tdeltaz = z1 - z2\n",
    "\n",
    "\tLstar1 = 42.76\n",
    "\tLstar2 = 42.33\n",
    "\n",
    "\tphistar1 = -3.17\n",
    "\tphistar2 = -2.86\n",
    "\n",
    "\tmLstar = (Lstar1 - Lstar2)/(deltaz)\n",
    "\tmphistar = (phistar1 - phistar2)/(deltaz)\n",
    "\n",
    "\tbLstar = Lstar1 - mLstar*z1\n",
    "\tbphistar = phistar1 - mphistar*z1\n",
    "\n",
    "\t#checked to make sure values give same as original\n",
    "\t#Lstartest1 = mLstar*z1 + bLstar\n",
    "\t#Lstartest2 = mLstar*z2 + bLstar\n",
    "\t#phistartest1 = mphistar*z1 + bphistar\n",
    "\t#phistartest2 = mphistar*z2 + bphistar\n",
    "\t#print(\"Lstar1, Lstar 2 = \",Lstartest1,Lstartest2)\n",
    "\t#print(\"phistar1, phistar2 = \",phistartest1,phistartest2)\n",
    "\n",
    "\tprint(\"Lstar for Lymanalpha: Lstar = z*\",mLstar,\" + \",bLstar)\n",
    "\tprint(\"phistar for Lymanalpha: phistar = z*\",mphistar,\" + \",bphistar)\n",
    "\n",
    "\tLstarLya = mLstar*z + bLstar\n",
    "\tphistarLya = mphistar*z + bphistar\n",
    "\n",
    "\tprint(\"LstarLya =\",LstarLya)\n",
    "\tprint(\"phistarLya =\",phistarLya)\n",
    "\n",
    "\tanswers = [LstarLya,phistarLya]\n",
    "\n",
    "\treturn answers\n",
    "\n",
    "\n",
    "def filt_int(filt):\n",
    "\n",
    "\tif filt==\"uband\":\n",
    "\n",
    "\t\t#read in each data file\n",
    "\t\tprint(\"u_filt has wavelengths 305.30 - 408.60 nm\")\n",
    "\t\t#this is in two columns; the left is wavelength, the right is throughput\n",
    "\t\tu_filt = numpy.loadtxt('ufilteredit.csv')\n",
    "\t\tprint(u_filt)\n",
    "\t\t#I shorten this to only the second column\n",
    "\t\tLSST_filt = u_filt\n",
    "\t\tLSSTfilt = u_filt[:,1]\n",
    "\n",
    "\tif filt==\"gband\":\n",
    "\t\t#read in each data file\n",
    "\t\tprint(\"g_filt has wavelengths 386.30 - 567.00 nm\")\n",
    "\t\t#this is in two columns; the left is wavelength, the right is throughput\n",
    "\t\tg_filt = numpy.loadtxt('gfilteredit.csv')\n",
    "\t\tprint(g_filt)\n",
    "\t\t#I shorten this to only the second column\n",
    "\t\tLSST_filt = g_filt\n",
    "\t\tLSSTfilt = g_filt[:,1]\n",
    "\n",
    "\tif filt==\"rband\":\n",
    "\n",
    "\t\t#read in each data file\n",
    "\t\tprint(\"r_filt has wavelengths 536.90 - 706.00 nm\")\n",
    "\t\t#this is in two columns; the left is wavelength, the right is throughput \n",
    "\t\tr_filt = numpy.loadtxt('rfilteredit.csv')\n",
    "\t\tprint(r_filt)\n",
    "\t\t#I shorten this to only the second column\n",
    "\t\tLSST_filt = r_filt\n",
    "\t\tLSSTfilt = r_filt[:,1]\n",
    "\n",
    "\tif filt==\"iband\":\n",
    "\n",
    "\t\t#read in each data file\n",
    "\t\tprint(\"i_filt has wavelengths 675.90 - 833.00 nm\")\n",
    "\t\t#this is in two columns; the left is wavelength, the right is throughput\n",
    "\t\ti_filt = numpy.loadtxt('ifilteredit.csv')\n",
    "\t\tprint(i_filt)\n",
    "\t\t#I shorten this to only the second column\n",
    "\t\tLSST_filt = i_filt\n",
    "\t\tLSSTfilt = i_filt[:,1]\n",
    "\n",
    "\tif filt==\"zband\":\n",
    "\n",
    "\t\t#read in each data file \n",
    "\t\tprint(\"z_filt has wavelengths 802.90 - 938.60 nm\")\n",
    "\t\t#this is in two columns; the left is wavelength, the right is throughput\n",
    "\t\tz_filt = numpy.loadtxt('zfilteredit.csv')\n",
    "\t\tprint(z_filt)\n",
    "\t\t#I shorten this to only the second column\n",
    "\t\tLSST_filt = z_filt\n",
    "\t\tLSSTfilt = z_filt[:,1]\n",
    "\n",
    "\tif filt==\"yband\":\n",
    "\n",
    "\t\t#read in each data file\n",
    "\t\tprint(\"y_filt has wavelengths 908.30 - 1099.60 nm\")\n",
    "\t\t#this is in two columns; the left is wavelength, the right is throughput\n",
    "\t\ty_filt = numpy.loadtxt('yfilteredit.csv')\n",
    "\t\tprint(y_filt)\n",
    "\t\t#I shorten this to only the second column\n",
    "\t\tLSST_filt = y_filt\n",
    "\t\tLSSTfilt = y_filt[:,1]\n",
    "\n",
    "\t#to find the midpoint, integrate over the entire filt, then find what value would give you half the integrated value\n",
    "\tstepint_LSSTfilt = scipy.integrate.cumtrapz(LSSTfilt)\n",
    "\tprint('stepint_LSSTfilt')\n",
    "\tprint(stepint_LSSTfilt)\n",
    "\tlen(stepint_LSSTfilt)\n",
    "\t#len gives total length, which is index number + 1\n",
    "\tlastnumber = stepint_LSSTfilt[len(stepint_LSSTfilt)-1]\n",
    "\tmidpoint = lastnumber/2\n",
    "\n",
    "\tdifference = abs(midpoint-stepint_LSSTfilt)\n",
    "\tprint(\"difference = \",difference)\n",
    "\tmindiff = min(difference)\n",
    "\tprint(\"mindiff = \",mindiff)\n",
    "\tmidindex = numpy.where(difference==mindiff)\n",
    "\tprint(\"midindex = \", midindex)\n",
    "\tprint(\"type is \",type(midindex))\n",
    "\t#have to change this from tuple to float - first find the actual value in index 0\n",
    "\tmidindex = midindex[0]\n",
    "\tmidindex = numpy.float(midindex)\n",
    "\tprint(\"midindex = \", midindex)\n",
    "\n",
    "\t#now I have to figure out if it is really the midpoint of the integrated throughput (tested for the z band)\n",
    "\t#Zfiltleft = zfilt[:658]\n",
    "\t#zfiltright = zfilt[658:]\n",
    "\t#a = scipy.integrate.trapz(zfiltleft)\n",
    "\t#b = scipy.integrate.trapz(zfiltright)\n",
    "\t#print(\"leftint = \",a)\n",
    "\t#print(\"rightint = \",b)\n",
    "\t#after testing 656,657,658, I concluded that 658 is the middle - I have to wrap my head around the indices in python and how trapz affects them\n",
    "\t#cumtrapz integrates from a point all the way to the right of the array\n",
    "\n",
    "\t#since 658 = midindex+1\n",
    "\tcenterindex = midindex+1\n",
    "\tcenterindex = int(centerindex)\n",
    "\tprint(\"centerindex = \",centerindex)\n",
    "\tprint(type(centerindex))\n",
    "\ttest = LSST_filt[centerindex]\n",
    "\tprint(type(test))\n",
    "\tlambdacenter = test[0]\n",
    "\tlambdacenter = numpy.float(lambdacenter)\n",
    "\tprint(\"median transmission wavelength of the \"+filt+\" = \",lambdacenter)\n",
    "\n",
    "\t#the central wavelength is actually called the median transmission wavelength of each band\n",
    "\n",
    "\treturn lambdacenter\n",
    "\n",
    "\n",
    "#the following function calculates the luminosity limit for a certain band with a certain detection limit\n",
    "def lumlim(z,em,filt):\n",
    "\n",
    "\tprint(\"This function will calculate the luminosity that corresponds to a 5 sigma detection, in erg/s:\")\n",
    "\t#print(\"INFO: 26.2 is AB magnitude for 5 sigma detection limit in the z band\") #this was used when I only had the z band\n",
    "\t#first I calculate the flux, then convert to flux density, then find the luminosity limit for the conditions printed above\n",
    "\n",
    "\tif filt==\"uband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 26.1\n",
    "\t\tlambdalow = 305.30 #in nm\n",
    "\t\tlambdahigh = 408.60 #in nm\n",
    "\n",
    "\tif filt==\"gband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 27.4\n",
    "\t\tlambdalow = 386.30 #in nm\n",
    "\t\tlambdahigh = 567.00 #in nm\n",
    "\n",
    "\tif filt==\"rband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 27.5\n",
    "\t\tlambdalow = 536.90 #in nm\n",
    "\t\tlambdahigh = 706.00 #in nm\n",
    "\n",
    "\tif filt==\"iband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 26.8\n",
    "\t\tlambdalow = 675.90 #in nm\n",
    "\t\tlambdahigh = 833.00 #in nm\n",
    "\n",
    "\tif filt==\"zband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 26.1\n",
    "\t\tlambdalow = 802.90 #in nm\n",
    "\t\tlambdahigh = 938.60 #in nm\n",
    "\n",
    "\tif filt==\"yband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 24.9\n",
    "\t\tlambdalow = 908.30 #in nm\n",
    "\t\tlambdahigh = 1099.60 #in nm\n",
    "\n",
    "\t#the following calculates values I use in and plug into the lumlim function\n",
    "\n",
    "\t#finds the flux density\n",
    "\tprint(\"ABmagnitude = -2.5*log10(fluxdensity/(3631 Jansky))\")\n",
    "\tprint(\"consequently:\")\n",
    "\t#print(\"fluxdensity = (10**(ABmagnitude/(-2.5)))*(3631 Janksy)\")\n",
    "\t#fluxdens = (10**(ABmag/(-2.5)))*3631 #outputs in Jansky\n",
    "\t#uses ABmag from earlier in this function\n",
    "\tfluxdens = 10**((ABmag+48.6)/(-2.5)) #outputs in erg/(s*Hz*(cm^2))\n",
    "\tprint(\"flux density =\",fluxdens,\"erg/(s*Hz*(cm^2))\")\n",
    "\n",
    "\t#finds the flux using the difference between the frequencies at each end of the band\n",
    "\tc = 2.9979*(10**17) #in nm/s\n",
    "\tdeltanu = c*((1/lambdalow)-(1/lambdahigh)) #the nm should cancel out\n",
    "\tprint(\"deltanu =\",deltanu,\"s^-1\")\n",
    "\tflux = fluxdens*deltanu#*(10**(-23)) #the extra factor converts from Janskys to ergs/(s*Hz*(cm^2))\n",
    "\tprint(\"flux =\",flux,\"erg/(s*(cm^2))\")\n",
    "\n",
    "\t#finds the luminosity distance\n",
    "\tlumdist = cosmo.luminosity_distance(z)\n",
    "\t#this outputs a special object that keeps track of units, so first I convert it to cm (cgs units), and then I convert it to a regular number\n",
    "\tlumdist_cgs = lumdist.to('cm')\n",
    "\tlumdist_unitless = lumdist_cgs.value\n",
    "\tprint(\"the luminosity distance for redshift z =\",z,\"is lumdist =\",lumdist_unitless,\"cm\")\n",
    "\n",
    "\t#finds the luminosity limit\n",
    "\tprint(\"Luminosity = 4*pi*(luminositydistance**2)*flux\")\n",
    "\tlumlimit = 4*numpy.pi*(lumdist_unitless**2)*flux\n",
    "\tprint(\"luminosity limit for 5 sigma detection of\",em,\"in \"+filt+\" band is\",lumlimit,\"ergs/s\")\n",
    "\n",
    "\t#using return makes the main output of this function the value of lumlimit so that I can use it to calculate other things when I call this function\n",
    "\treturn lumlimit\n",
    "\n",
    "\n",
    "#the following defines and plots the Schechter luminosity function for a chosen emission line along with a separate plot of the number density\n",
    "#it also calculates the number density above a luminosity limit calculated for a specific band in the previous function\n",
    "#to use this function, type schechter_LF(z=redshifttoplot) or with additional parameters you want to change inside the ()\n",
    "def schechter_LF(z,lambdaemitted,alpha,Lstar0,betaL,phistar0,betaphi,zpaper,param,fluxscale,em,filt,style = \"\"):\n",
    "\n",
    "\tprint(\"alpha = \", alpha)\n",
    "\tprint(\"z = \", z)\n",
    "\n",
    "\ttest = numpy.arange(30,55,0.01)\n",
    "\tL = (10**test)\n",
    "\n",
    "\t#I have two different parametrizations for Lstar and phistar; the first one is from Comparat et al 2016, and the second one is from Sobral et al 2015\n",
    "\n",
    "\tif param == \"first\":\n",
    "\n",
    "\t\tLstar = fluxscale*Lstar0*((1+z)**betaL)\n",
    "\t\tprint(\"Lstar = \", Lstar)\n",
    "\n",
    "\t\tphistar = phistar0*((1+z)**betaphi)\n",
    "\t\tprint(\"phistar = \", phistar)\n",
    "\n",
    "\tif param == \"second\":\n",
    "\n",
    "\t\tLstar = 10**(0.45*z + Lstar0)\n",
    "\t\tprint(\"Lstar = \", Lstar)\n",
    "\n",
    "\t\tphistar = 10**(-0.38*(z**2) + z + phistar0)\n",
    "\t\tprint(\"phistar = \", phistar)\n",
    "\n",
    "\tif param == \"third\":  #THIS STILL NEEDS TO BE FIXED\n",
    "\n",
    "\t\tanswersLya = lineqLya(z)\n",
    "\n",
    "\t\tLstar = answersLya[0] #this is linearly parametrized from Ciardullo+ 2012\n",
    "\t\tprint(\"Lstar = \",Lstar)\n",
    "\n",
    "\t\tphistar = answersLya[1] #this is linearly parametrized from Ciardullo+ 2012\n",
    "\t\tprint(\"phistar = \",phistar)\n",
    "\n",
    "\tphi = phistar*((L/Lstar)**(alpha+1))*(numpy.e**(-L/Lstar))\n",
    "\tprint(\"Note: Each paper uses a slightly different convention; I decided to consolidate them with the following:\")\n",
    "\tprint(\"using the LF equation with the alpha+1 exponent as follows: \")\n",
    "\tprint(\"phi = phistar*((L/Lstar)**(alpha+1))*(numpy.e**(-L/Lstar))\")\n",
    "\n",
    "\t#this deletes parts of the arrays that are so small python counts them as zero; otherwise, I would not be able to take the logarithm of the array\n",
    "\tL = L[phi!=0]\n",
    "\tphi = phi[phi!=0]\n",
    "\ttestarray = test[phi!=0]\n",
    "\n",
    "\tlog10L = numpy.log10(L)\n",
    "\tlog10phi = numpy.log10(phi)\n",
    "\n",
    "\t#the following calculates values I use in and plug into the lumlim function\n",
    "\n",
    "\tif filt==\"uband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 26.1\n",
    "\t\tlambdalow = 305.30 #in nm\n",
    "\t\tlambdahigh = 408.60 #in nm\n",
    "\t\tlambdacenter = filt_int(filt = filt)\n",
    "\n",
    "\tif filt==\"gband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 27.4\n",
    "\t\tlambdalow = 386.30 #in nm\n",
    "\t\tlambdahigh = 567.00 #in nm\n",
    "\t\tlambdacenter = filt_int(filt = filt)\n",
    "\n",
    "\tif filt==\"rband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 27.5\n",
    "\t\tlambdalow = 536.90 #in nm\n",
    "\t\tlambdahigh = 706.00 #in nm\n",
    "\t\tlambdacenter = filt_int(filt = filt)\n",
    "\n",
    "\tif filt==\"iband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 26.8\n",
    "\t\tlambdalow = 675.90 #in nm\n",
    "\t\tlambdahigh = 833.00 #in nm\n",
    "\t\tlambdacenter = filt_int(filt = filt)\n",
    "\n",
    "\tif filt==\"zband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 26.1\n",
    "\t\tlambdalow = 802.90 #in nm\n",
    "\t\tlambdahigh = 938.60 #in nm\n",
    "\t\tlambdacenter = filt_int(filt = filt)\n",
    "\n",
    "\tif filt==\"yband\":\n",
    "\t\t#ABmag is the coadded depth for a 5 sigma magnitude limit in this filt\n",
    "\t\tABmag = 24.9\n",
    "\t\tlambdalow = 908.30 #in nm\n",
    "\t\tlambdahigh = 1099.60 #in nm\n",
    "\t\tlambdacenter = filt_int(filt = filt)\n",
    "\n",
    "\t#lambdaemitted = #will be one of the three emission lines #in nm\n",
    "\tfiltendlow = (lambdalow/lambdaemitted)-1\n",
    "\tfiltendhigh = (lambdahigh/lambdaemitted)-1\n",
    "\tfiltcenter = (lambdacenter/lambdaemitted)-1\n",
    "\tdeltafilt = filtendhigh-filtendlow\n",
    "\tprint(\"for\",em,\"deltafilt =\",deltafilt)\n",
    "\n",
    "\t#finds luminosity limit in center of z band\n",
    "\tcenter = lumlim(z = filtcenter,em = em,filt = filt)\n",
    "\n",
    "\tlumarray = numpy.full(15,numpy.log10(center))\n",
    "\tyarray = numpy.arange(-10,5,1)\n",
    "\n",
    "\tfigure(1)\n",
    "\tplot(log10L,log10phi,style,label = zpaper)\n",
    "\tplot(lumarray,yarray,style+\"--\")\n",
    "\txlim(40,45)\n",
    "\tylim(-7,0)\n",
    "\tlegend(loc = \"upper right\")\n",
    "\t#LaTeX is used with $ signs in the titles below\n",
    "\txlabel(\"$\\log_{10}(L [ergs/s])$\")\n",
    "\tylabel(\"$\\log_{10}(\\phi [\\t{Mpc}^{-3}])$\")\n",
    "\ttitle(\"Schechter Luminosity Function\") #emline+\n",
    "\tprint(\"z = \", z, \" and alpha = \", alpha, \" are plotted, \")\n",
    "\t#saves image\n",
    "\n",
    "\tprint(\"now the number density is calculated by integrating the LF using cumtrapz:\")\n",
    "\t#cumptrapz integrates the opposite way than I need to integrate, so I flip it twice in the process\n",
    "\tphiflip = phi[::-1]\n",
    "\tphiflipint = scipy.integrate.cumtrapz(phiflip,x=testarray)\n",
    "\tnum_dens = phiflipint[::-1]\n",
    "\n",
    "\tLmin = numpy.delete(L,(len(L)-1))\n",
    "\t#shifted_Lmin = Lmin*fluxscale\n",
    "\tlog10Lmin = numpy.log10(Lmin)\n",
    "\t#log10shifted_Lmin = numpy.log10(shifted_Lmin)\n",
    "\n",
    "\tfigure(2)\n",
    "\tplot(log10Lmin,num_dens,style,label = zpaper)\n",
    "\tplot(lumarray,yarray,style+\"--\")\n",
    "\tlegend(loc = \"upper right\")\n",
    "\txlim(40,44)\n",
    "\tylim(0,0.03)\n",
    "\t#LaTeX is used with $ signs in the titles below\n",
    "\txlabel(\"$\\log_{10}(L_{min} [ergs/s])$\")\n",
    "\tylabel(\"$\\log_{10}(\\phi [\\t{Mpc}^{-filt3}])$\")\n",
    "\ttitle(\"Number Density\") #emline+\n",
    "\t#saves image\n",
    "\tpyplot.savefig('/home/lanaeid/Desktop/fig2.png',bbox_inches = 'tight')\n",
    "\n",
    "\t#the following finds the comoving number density above a certain detection limit (shot noise limit),\n",
    "\t#which I calculated in the previous part of the code (also shown below), where I saved the variable named center\n",
    "\n",
    "\t#runs through the LF code for chosen emission line\n",
    "\t#then gets number density above the luminosity limit using the center value\n",
    "\n",
    "\t#first, I have to shorten the phi array to contain only values above the luminosity limit\n",
    "\n",
    "\t#uses astropy.cosmology.Planck15 to find comoving volume in shell between redshifts at each end of the z filt\n",
    "\tcomovingvolmin = cosmo.comoving_volume(filtendlow) #units are in Mpc^3\n",
    "\tcomovingvolmax = cosmo.comoving_volume(filtendhigh) #units are in Mpc^3\n",
    "\tcomovingvol = comovingvolmax-comovingvolmin\n",
    "\tcomovingvol = comovingvol.value\n",
    "\tprint(\"comovingvol =\",comovingvol,\"Mpc^3\")\n",
    "\n",
    "\t#shortens the array to be above the luminosity limit, then integrates to get comoving number density\n",
    "\tphilim = phi[L>center]\n",
    "\ttestarraydos = test[L>center]\n",
    "\tcomovingphi = scipy.integrate.trapz(philim,x=testarraydos)\n",
    "\tprint(\"comovingphi =\",comovingphi,\"Mpc^-3\") #units?  \n",
    "\n",
    "\t#finds total number of galaxies and areal number density\n",
    "\ttotalnumgalaxies = comovingphi*comovingvol\n",
    "\tarealphi = totalnumgalaxies/(4*numpy.pi)\n",
    "\tprint(\"arealphi =\",arealphi,\"steradian^-1\")\n",
    "\n",
    "\tshow()\n",
    "if emline == \"test\":\n",
    "\n",
    "\t#this was used to print out data and plots for the LSST-DESC conference poster\n",
    "\n",
    "\tprint(\"This option finds comoving and areal number densities for:\")\n",
    "\tprint(\"[OII] 3726/3729 and [OIII] 5007 from Comparat et al 2016\")\n",
    "\tprint(\"Halpha 6563 from Sobral et al 2013\")\n",
    "\n",
    "\tlambda_OII = 372.7\n",
    "\tlambda_OIII = 500.7\n",
    "\tlambda_Halpha = 656.3\n",
    "\tlambda_Lymanalpha = 121.6\n",
    "\n",
    "\t#the following calculates values I use in and plug into the lumlim function\n",
    "\tlambdalow = 818.95 #in nm\n",
    "\tlambdahigh = 921.15 #in nm\n",
    "\tlambdacenter = (lambdalow+lambdahigh)/2 #in nm #NO, THIS IS NOT CORRECT, FIX THIS\n",
    "\t#lambdaemitted = #will be one of the three emission lines #in nm\n",
    "\t#zendlow = (lambdalow/lambdaemitted)-1\n",
    "\t#zendhigh = (lambdahigh/lambdaemitted)-1\n",
    "\t#zcenter = (lambdacenter/lambdaemitted)-1\n",
    "\n",
    "\t#finds redshift of emission line at center of z band -> use to find LF\n",
    "\tzOII = (lambdacenter/lambda_OII)-1\n",
    "\tzOIII = (lambdacenter/lambda_OIII)-1\n",
    "\tzHalpha = (lambdacenter/lambda_Halpha)-1\n",
    "\tzLymanalpha = (lambdacenter/lambda_Lymanalpha)-1\n",
    "\n",
    "\t#schechter_LF(z=zOII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = \"first\",zpaper = \"z = \"+str(zOII)+\" Comparat+ 2016\",fluxscale = 1,style = \"r\")\n",
    "\t#schechter_LF(z=zOIII,alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = \"first\",zpaper = \"z = \"+str(zOIII)+\" Comparat+ 2016\",fluxscale = 1,style = \"g\")\n",
    "\t#schechter_LF(z = zHalpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = \"second\",zpaper = \"z = \"+str(zHalpha)+\" Sobral+ 2013\",fluxscale = 1,style = \"b\")\n",
    "\n",
    "\tschechter_LF(z=zOII,lambdaemitted = lambda_OII,alpha = -1.46,Lstar0 = 10**41.1,betaL = 2.33,phistar0 = 10**(-2.4),betaphi = -0.73,param = \"first\",zpaper = \"[OII] z = 1.33 Comparat+ 2016\",fluxscale = 1,em = \"[OII]\",filt = filt,style = \"r\")\n",
    "\tschechter_LF(z=zOIII,lambdaemitted = lambda_OIII, alpha = -1.83,Lstar0 = 10**41.42,betaL = 3.91,phistar0 = 10**(-3.41),betaphi = -0.76,param = \"first\",zpaper = \"[OIII] z = 0.74 Comparat+ 2016\",fluxscale = 1,em = \"[OIII]\",filt = filt,style = \"g\")\n",
    "\tschechter_LF(z=zHalpha,lambdaemitted = lambda_Halpha,alpha = -1.6,Lstar0 = 41.87,betaL = 0,phistar0 = -3.18,betaphi = 0,param = \"second\",zpaper = r\"H$\\alpha$ z = 0.33 Sobral+ 2013\",fluxscale = 1, em = \"Halpha\",filt = filt,style = \"b\")\n",
    "\t#the following is also from the separate linearequation code that I tried to put in this one, but it throws back errors every time I use the \"third\" option, so I have to fix that later\n",
    "\tschechter_LF(z=zLymanalpha,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**44.0057781641604,betaL = 0,phistar0 = 10**(-4.068119141604011),betaphi = 0,param = \"first\",zpaper = r\"Ly$\\alpha$ z = 6.155 Ciardullo+ 2012\",fluxscale = 1,em = \"Lymanalpha\",filt = filt,style = \"y\")\n",
    " \t#used LaTex above for legend\n",
    " \t#python recognizes LaTeX instead of thinking they're escape characters if I write it as a \"raw\" string, denoting it with an r in the beginning of the line\n",
    "\n",
    "\n",
    "if emline == \"Lymanalpha\":\n",
    "\n",
    "\tprint(\"This following plot will contain a Lymanalpha LF from Ciardullo et al 2012.  \")\n",
    " \t#either this needs to be combined with the \"test\" option, or I need to edit the rest of these emline options to be usable with all the different parameters\n",
    "\n",
    "\tlambda_OII = 372.7\n",
    "\tlambda_OIII = 500.7\n",
    "\tlambda_Halpha = 656.3\n",
    "\tlambda_Lymanalpha = 121.6\n",
    "\n",
    "\t#the following calculates values I use in and plug into the lumlim function\n",
    "\tlambdalow = 818.95 #in nm\n",
    "\tlambdahigh = 921.15 #in nm\n",
    "\tlambdacenter = (lambdalow+lambdahigh)/2 #in nm #NO, THIS IS NOR CORRECT, FIX IT\n",
    "\t#lambdaemitted = #will be one of the three emission lines #in nm\n",
    "\t#zendlow = (lambdalow/lambdaemitted)-1\n",
    "\t#zendhigh = (lambdahigh/lambdaemitted)-1\n",
    "\t#zcenter = (lambdacenter/lambdaemitted)-1\n",
    "\n",
    "\t#finds redshift of emission line at center of z band -> use to find LF\n",
    "\tzOII = (lambdacenter/lambda_OII)-1\n",
    "\tzOIII = (lambdacenter/lambda_OIII)-1\n",
    "\tzHalpha = (lambdacenter/lambda_Halpha)-1\n",
    "\tzLymanalpha = (lambdacenter/lambda_Lymanalpha)-1\n",
    "\n",
    "\t#for the first redshift below\n",
    "\tprint(\"z = 3.113 and alpha = -1.65 are plotted, \")\n",
    "\tprint(\"the values with errors given in the Ciardullo et al 2012 paper are:\")\n",
    "\tprint(\"Lstar_exp =   42.76 (+0.10 -0.10)\")\n",
    "\tprint(\"phistar_exp = -3.17 (+0.05 -0.05)\")\n",
    "\n",
    "\t#for the second redshift below\n",
    "\tprint(\"z = 2.063 and alpha = -1.65 are plotted, \")\n",
    "\tprint(\"the values with errors given in the Ciardullo et al 2012 paper are:\")\n",
    "\tprint(\"Lstar_exp =   42.33 (+0.12 -0.12)\")\n",
    "\tprint(\"phistar_exp = -2.86 (+0.05 -0.05)\")\n",
    "    \n",
    "\tschechter_LF(z=3.113,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**42.76,betaL = 0,phistar0 = 10**(-3.17),betaphi = 0,param = \"first\",zpaper = r\"Ly$\\alpha$ z = 3.113 Ciardullo+ 2012\",fluxscale = 1,em = \"Lymanalpha\",filt = filt, style = \"c\")\n",
    "\tschechter_LF(z=2.063,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**42.33,betaL = 0,phistar0 = 10**(-2.86),betaphi = 0,param = \"first\",zpaper = r\"Ly$\\alpha$ z = 2.063 Ciardullo+ 2012\",fluxscale = 1,em = \"Lymanalpha\",filt = filt, style = \"m\")\n",
    "\n",
    "\t#the following is also from the separate linearequation code that I tried to put in this one, but it throws back errors every time I use the \"third\" option, so I have to fix that later\n",
    "\t#THE OPTION \"first\" ACTUALLY WORKS, BUT ONLY IF YOU ALSO INPUT THE VALUES FOR LSTAR0 AND PHISTAR0 TO BYPASS THE ACTUAL THING I WAS TRYING TO DO\n",
    "\tschechter_LF(z=6.155016447368421,lambdaemitted = lambda_Lymanalpha,alpha = -1.65,Lstar0 = 10**44.0057781641604,betaL = 0,phistar0 = 10**(-4.068119141604011),betaphi = 0,param = \"first\",zpaper = r\"Ly$\\alpha$ z = 6.155 Ciardullo+ 2012\",fluxscale = 1,em = \"Lymanalpha\",filt = \"zband\",style = \"y\")\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.6.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
