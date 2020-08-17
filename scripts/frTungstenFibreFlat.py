#!/usr/bin/python3

import os, argparse, sys

import numpy as np
from astropy.io import fits
from scipy.signal import savgol_filter


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Apply fibre flat in Frodo L2 DpOL.')
    parser.add_argument('targetfile', action='store', help='Extracted science target FITS from frextract, mandatory')
    parser.add_argument('tungstenfile', action='store', help='Extracted tungsten continuum FITS from frextract, mandatory')
    parser.add_argument('outputfilename', action='store', help='Output file name, mandatory')

    parser.add_argument('-l', dest='low', action='store_true', help='Apply extra smoothing that applies only to the low res grating. (default: Off)')
    parser.add_argument('-e', dest='equalize', action='store_true', help='Do fibre throughput equalisation. Only needed if you have not run frcorrect on the wlamp input. (default: Off)')
    parser.add_argument('-d', dest='debug', action='store_true', help='Turn on debug mode (default: Off)')
    parser.add_argument('-v', dest='verbose', action='store_true', help='Turn on verbose mode (default: Off)')
    # -h is always added automatically by argparse1

    args = parser.parse_args()

    if args.verbose:
      print (args)


# Read the 144 extracted fibres from rhe Tungsten file
hdul = fits.open(args.tungstenfile)

if args.equalize: 
  # First the entire array in collapsed in the spectral direction to get a total counts in each
  # fibre and the image scaled by the inverse of that in order to match the total throughput
  # of each fibre
  spectralCollapse = np.zeros((144,1))
  spectralCollapse[:,0] = np.mean(hdul[0].data,axis=1)
  hdul[0].data = hdul[0].data / spectralCollapse
  if args.debug :
    outname = 'equalized_%s' % (args.tungstenfile)
    fits.writeto(outname, hdul[0].data, header=hdul[0].header, overwrite=True)

# Collapse it in the spatial direction to get the average spectral shape of the tungsten.
# Normalise to a smoothed version of itself in order to get the fibre flat that shows only
# small scale deviations from the mean tungsten shape.
spatialCollapsed = np.mean(hdul[0].data,axis=0)
spatialCollapsedSmooth = savgol_filter(spatialCollapsed, 101, 5)

if (args.low) : 
  # This patch is for Low only. It cleans the LOW config a bit better but not a huge
  # difference and not mandatory
  spatialCollapsedSmooth[0:1250] = savgol_filter(spatialCollapsed[0:1250], 301, 3)

tungstenNormalized = hdul[0].data / spatialCollapsedSmooth

# For debug, save off a copy of the collapsed and collapsedSmoothed versions.
if args.debug :
  outname = 'collased_%s' % (args.tungstenfile)
  fits.writeto(outname, spatialCollapsed, header=hdul[0].header, overwrite=True)
  outname = 'smooth_collased_%s' % (args.tungstenfile)
  fits.writeto(outname, spatialCollapsedSmooth, header=hdul[0].header, overwrite=True)
  outname = 'fibreflat_%s' % (args.tungstenfile)
  fits.writeto(outname, tungstenNormalized, header=hdul[0].header, overwrite=True)

# Close the tungsten
hdul.close()

# Flat field the science data with the tungsten fibre flat just created
hdul = fits.open(args.targetfile)
hdul[0].data = hdul[0].data / tungstenNormalized

#outname = 'tungstened_%s' % (args.targetfile)
fits.writeto(args.outputfilename, hdul[0].data, header=hdul[0].header, overwrite=True)

hdul.close()

