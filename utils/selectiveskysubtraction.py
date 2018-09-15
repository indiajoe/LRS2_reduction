#!/usr/bin/env python
""" This script is to selectively subtract sky from certain lenslets
The sky is then alligned and subtracted from other lenslets """

from astropy.io import fits
import numpy as np
from scipy import optimize
from scipy import interpolate
from scipy import signal
import sys

FeFitsFile = sys.argv[1] #'FeRpses20170624T064517.2_066_sci_L.fits'#'FeRpses20170624T063753.9_066_sci_L.fits'
hdulist = fits.open(FeFitsFile)

OriginalData = hdulist[0].data.copy()
OriginalHeader = hdulist[0].header.copy()
NewData = hdulist[0].data.copy()
# Top Half
IntensitySortedIndex1 = np.argsort(np.percentile(NewData[:140,:],10,axis=1))
IntensitySortedIndex2 = np.argsort(np.percentile(NewData[140:,:],10,axis=1))
NoOfFaintSky = 10

SkyTemplate1 = np.mean(NewData[:140,:][IntensitySortedIndex1[:NoOfFaintSky],:],axis=0)
SkyTemplate2 = np.mean(NewData[140:,:][IntensitySortedIndex2[:NoOfFaintSky],:],axis=0)
# Now fit and subtract the sky template form every lenslet

# First do a simple subtract
hdulist[0].data[:140,:] = NewData[:140,:] - SkyTemplate1
hdulist[0].data[140:,:] = NewData[140:,:] - SkyTemplate2
hdulist[0].header.add_history('Direct Sky Subtraction')
hdulist.writeto('SkySubDirect'+FeFitsFile)

# Scale the sky and subtract


def subtractcontinuum(spec):
    """ medianfilter out the continuum """
    continuum = signal.medfilt(spec,33)
    diff = spec - continuum
    # set netgative values to zero so that fit is insensitive to absroption lines
    diff[diff<0] = 0
    return diff

def errorfuncwithscale(p,sky,spec, removecontinuum = True):
    """ Returns residue of spec - sky """
    if removecontinuum:
        sky = subtractcontinuum(sky)
        spec = subtractcontinuum(spec)
    return spec - sky*p[0]

def errorfuncwithshift(p,sky,spec, removecontinuum = True):
    """ Returns residue of spec - sky """
    if removecontinuum:
        sky = subtractcontinuum(sky)
        spec = subtractcontinuum(spec)
    tcksky = interpolate.splrep(np.arange(len(sky)),sky)
    shiftedsky = interpolate.splev(np.arange(len(spec))+p[0], tcksky)
    return spec - shiftedsky*p[1]


# Sky removal with just scaling
for i in range(OriginalData.shape[0]):
    SkyTemplate = SkyTemplate1 if (i < 140) else SkyTemplate2
    p,ier = optimize.leastsq(errorfuncwithscale,np.array([1.]),args=(SkyTemplate,OriginalData[i,:]))
    NewData[i,:] = errorfuncwithscale(p,SkyTemplate,OriginalData[i,:], removecontinuum=False)

hdulist[0].data = NewData
hdulist[0].header = OriginalHeader
hdulist[0].header.add_history('Scaled Sky Subtraction')
hdulist.writeto('SkySubScaled'+FeFitsFile)


# Sky removal with just scaling and shifting
for i in range(OriginalData.shape[0]):
    SkyTemplate = SkyTemplate1 if (i < 140) else SkyTemplate2
    p,ier = optimize.leastsq(errorfuncwithshift,np.array([0.,1.]),args=(SkyTemplate,OriginalData[i,:]))
    NewData[i,:] = errorfuncwithshift(p,SkyTemplate,OriginalData[i,:], removecontinuum=False)

hdulist[0].data = NewData
hdulist[0].header = OriginalHeader
hdulist[0].header.add_history('Shifted Scaled Sky Subtraction')
hdulist.writeto('SkySubShifted'+FeFitsFile)



    
