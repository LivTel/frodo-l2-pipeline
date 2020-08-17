#!/usr/bin/python3

# Error exit codes
# 0	success
# 1	tungsten file not found
# 2	target file not found
# 3	tungsten file failed to be read
# 4	target file failed to be read

import os, argparse, sys

import numpy as np
from astropy.io import fits
from scipy.signal import savgol_filter


# Command line parser
parser = argparse.ArgumentParser(description='Apply fibre flat in Frodo L2 DpOL.')
parser.add_argument('targetfile', action='store', help='Extracted science target FITS from frextract, mandatory')
parser.add_argument('tungstenfile', action='store', help='Extracted tungsten continuum FITS from frextract, mandatory')
parser.add_argument('outputfilename', action='store', help='Output file name, mandatory')

parser.add_argument('-l', dest='low', action='store_true', help='Apply extra smoothing that applies only to the low res grating. (default: Off)')
parser.add_argument('-e', dest='equalize', action='store_true', help='Do fibre throughput equalisation. Only needed if you have not run frcorrect on the wlamp input. (default: Off)')
parser.add_argument('-d', dest='debug', action='store_true', help='Turn on debug mode. Saves interim products to disk. (default: Off)')
parser.add_argument('-v', dest='verbose', action='store_true', help='Turn on verbose mode (default: Off)')
# -h is always added automatically by argparse1

args = parser.parse_args()

if args.verbose:
  print ("Command line options")
  print (args)

if not os.path.isfile(args.tungstenfile) :
  print ("Requested tungsten spectrum not found")
  sys.exit(1)

if not os.path.isfile(args.targetfile) :
  print ("Requested target spectrum not found")
  sys.exit(2)


# Read the 144 extracted fibres from rhe Tungsten file
try: 
  hdul = fits.open(args.tungstenfile)
except:
  print ("Failed to read tungsten spectrum file")
  sys.exit(3)

if args.equalize: 
  print ("Applying fibre throughput equalisation")
  # First the entire array in collapsed in the spectral direction to get a total counts in each
  # fibre and the image scaled by the inverse of that in order to match the total throughput
  # of each fibre
  #
  # The frodo L2 pipeline already has a function to do this (correct_throughput) and 
  # generally I would expect you to have used that rather than this trivial implementation.
  spectralCollapse = np.zeros((144,1))
  spectralCollapse[:,0] = np.mean(hdul[0].data,axis=1)
  try:
    hdul[0].data = hdul[0].data / spectralCollapse
  except:
    print("Error equalizing fibre throughoputs. Maybe divide-by-zero, divide-by-nan?")
  if args.debug :
    outname = 'equalized_%s' % (args.tungstenfile)
    fits.writeto(outname, hdul[0].data, header=hdul[0].header, overwrite=True)

# Collapse it in the spatial direction to get the average spectral shape of the tungsten.
# Normalise to a smoothed version of itself in order to get the fibre flat that shows only
# small scale deviations from the mean tungsten shape.
print ("Filtering of tungsten spectral shape hard coded to Savitzky-Golay(101 pix, order 5)")
spatialCollapsed = np.mean(hdul[0].data,axis=0)
spatialCollapsedSmooth = savgol_filter(spatialCollapsed, 101, 5)

if (args.low) : 
  # This patch is for Low only. It cleans the LOW config a bit better but not a huge
  # difference and not mandatory
  print ("Applying second low pass filter to spectrum. Recommended for Red-Low only")
  spatialCollapsedSmooth[0:1250] = savgol_filter(spatialCollapsed[0:1250], 301, 3)

try:
  tungstenNormalized = hdul[0].data / spatialCollapsedSmooth
except:
  print("Error normalizing out the spectral shape. Maybe divide-by-zero, divide-by-nan?")

# For debug, save off a copy of the collapsed and collapsedSmoothed versions.
if args.debug :
  print ("Write debug spectra into working directory");
  outname = 'collased_%s' % (args.tungstenfile)
  fits.writeto(outname, spatialCollapsed, header=hdul[0].header, overwrite=True)
  outname = 'smooth_collased_%s' % (args.tungstenfile)
  fits.writeto(outname, spatialCollapsedSmooth, header=hdul[0].header, overwrite=True)
  outname = 'fibreflat_%s' % (args.tungstenfile)
  fits.writeto(outname, tungstenNormalized, header=hdul[0].header, overwrite=True)

# Close the tungsten
hdul.close()

#
# All the above created the flat. Finally, apply it to the target file.
# Flat field the science data with the tungsten fibre flat just created
#
try:
  hdul = fits.open(args.targetfile)
except:
  print ("Failed to read target spectrum file")
  sys.exit(4)

try:
  hdul[0].data = hdul[0].data / tungstenNormalized
except:
  print("Error applying tungsten correction to data file")

fits.writeto(args.outputfilename, hdul[0].data, header=hdul[0].header, overwrite=True)

hdul.close()

sys.exit(0)
