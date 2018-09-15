#!/usr/bin/env python
""" This script is to plot the hexagonal IFU grid map of the reduced data
It collapses data in wavelegnth region to create 2D honeycomb grid
                                              -indiajoe   
"""

import argparse
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import numpy as np
from astropy.io import fits

MappingFile = '/home/joe/OldHome/GitHub/LRS2_reduction/lrs2_config/mapping_files/LRS2_R_NR_mapping.txt'

#Hex Polygon centered on (0,0)
apothem = 0.59/2.
HexVertices = np.array([[apothem,apothem/np.sqrt(3)],
                        [0,2*apothem/np.sqrt(3)],
                        [-apothem,apothem/np.sqrt(3)],
                        [-apothem,-apothem/np.sqrt(3)],
                        [0,-2*apothem/np.sqrt(3)],
                        [apothem,-apothem/np.sqrt(3)]])



def generate_IFUplot(axis,IFUPatchcollection,datatocolor):
    """ returns the IFU image with colors as per datatocolor vector values """
    colormap = plt.get_cmap('viridis')

    scaleddata = datatocolor-np.min(datatocolor)
    scaleddata = np.round(255*(scaleddata/np.max(scaleddata))).astype(np.int)

    colors = colormap(scaleddata) 
    IFUPatchcollection.set_color(colors)

    axis.add_collection(IFUPatchcollection)
    return axis

def CreateIFUPolygonPatchCollection(PosXY,pAngleDeg=0):
    """ Returns the matplotlib polygon path collection """
    pAngle = np.deg2rad(pAngleDeg)
    RotationM = np.array([[np.cos(pAngle),-np.sin(pAngle)],
                          [np.sin(pAngle),np.cos(pAngle)]])
    IFUs = []
    for i in sorted(PosXY):
        vertices = HexVertices + PosXY[i]
        rotatedvert = np.dot(vertices,RotationM)
        hexIFU = Polygon(rotatedvert,closed=True)
        IFUs.append(hexIFU)

    IFUPatchcollection = PatchCollection(IFUs)
    return IFUPatchcollection

def LoadIFUskymappings(MappingFile):
    """ Returns the IFU positions on sky by loading the MappingFile """
    PosXY = {}
    with open(MappingFile) as mfile:
        # Skip first 10 lines
        for i in range(10): _ = mfile.next()
        for line in mfile:
            if line[0] == '#' : continue
            line = line.rstrip().split()
            PosXY[int(line[4])] = np.array([float(line[1]),float(line[2])])
    return PosXY

def get_dataslice(Fitsfilename,Wstart,Wend):
    """ Returns the data portion sliced in the range Wstart to Wend """
    hdulist = fits.open(Fitsfilename)
    WFirst = hdulist[0].header['CRVAL1']
    DeltaW = hdulist[0].header['CDELT1']
    StartIndx = int((Wstart - WFirst)/DeltaW)
    EndIndx = int((Wend - WFirst)/DeltaW)
    return hdulist[0].data[:,StartIndx:EndIndx]

def parse_args():
    """ Parses the command line input arguments """
    parser = argparse.ArgumentParser(description="Plot LRS2 IFU map")
    parser.add_argument('FibreExtractedFits', type=str,
                        help="The Fibre extracted Wavelenght Calibrated Fits File (Fe*)")
    parser.add_argument('WavelengthStart', type=float,
                        help="Starting Wavelength to collapse")
    parser.add_argument('WavelengthEnd', type=float,
                        help="Ending Wavelength to collapse")
    args = parser.parse_args()
    return args
    
def main():
    """ Standalone Interactive Line Identify Tool """
    args = parse_args()    
    FibreExtractedFits_fname = args.FibreExtractedFits
    WStart = args.WavelengthStart
    WEnd = args.WavelengthEnd
    print('Generating map of {0}..'.format(FibreExtractedFits_fname))
    print('Collapsing data in wavelgnth range {0} to {1}'.format(WStart,WEnd))
    PAngleDeg = fits.getval(FibreExtractedFits_fname,'PARANGLE')
    print('Rotating by Paraletic angle {0}'.format(PAngleDeg))

    dataslice = get_dataslice(FibreExtractedFits_fname,WStart,WEnd)
    collapsedata = np.sum(dataslice,axis=1)
    PosXY = LoadIFUskymappings(MappingFile)
    IFUPatch = CreateIFUPolygonPatchCollection(PosXY,pAngleDeg=-PAngleDeg)

    fig,axis = plt.subplots(1)
    axis = generate_IFUplot(axis,IFUPatch,collapsedata)
    axis.autoscale_view()
    plt.xlim(plt.xlim()[::-1])
    plt.xlabel('RA (")')
    plt.ylabel('Dec (")')
    plt.show()

if __name__ == "__main__":
    main()
