#!/usr/bin/env python
""" This script is to plot a spacially averaged plot of spectrum """
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from plot_IFUmap import LoadIFUskymappings
import numpy as np
import sys

NoNN = 7 # Number of Nearest neighbours to average

Fitsfilename = sys.argv[1]

# Coordinates to average 0,0 1,2 etc..
XYlist = [(float(coo.split(',')[0]), float(coo.split(',')[1])) for coo in sys.argv[2:]]


MappingFile = '/home/joe/OldHome/GitHub/LRS2_reduction/lrs2_config/mapping_files/LRS2_R_NR_mapping.txt'

PosXY = LoadIFUskymappings(MappingFile)

XY = np.array([PosXY[i] for i in sorted(PosXY)])

IFUPosTree = cKDTree(XY)


hdulist = fits.open(Fitsfilename)
WFirst = hdulist[0].header['CRVAL1']
DeltaW = hdulist[0].header['CDELT1']
Unit = hdulist[0].header['CUNIT1']
Warray = WFirst + DeltaW*np.arange(hdulist[0].data.shape[1])

for X,Y in XYlist:
    dd, ii = IFUPosTree.query(np.array([X,Y]),NoNN)


    Data2DArray = hdulist[0].data[ii,:]

    AverageSpectrum = np.average(Data2DArray,axis=0,weights=1/dd)

    plt.plot(Warray,AverageSpectrum,label='{0},{1}'.format(X,Y))

plt.legend()
plt.ylabel('Counts')
plt.xlabel(Unit)
plt.show()
