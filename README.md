frodo-l2-pipeline
=================

# Overview

This package provides a reduction pipeline to reduce data taken with the FRODOSpec spectrograph. 

The core procedures (`src/`) have been written in C, and are mostly ported from the FRODOSpec codebase. The wrapper for the 
pipeline, `scripts/L2_exec`, is written in cshell. A more in depth description of the methodology is given in Barnsley, 2012.

# Installation

1. Clone the repository

2. Edit `scripts/L2_setup` and set the `L2_BASE_DIR` variable to the root directory of the repository.

3. Set up the environment by sourcing the `scripts/L2_setup` file

4. Edit `src/Makefile` and set the `LIBS` and `INCLUDES` parameters accordingly. These paths need to include locations of the GSL and CFITSIO headers/libraries. 

5. Make the binaries and library: `src/make all`

6. Red arm requires Python 3, `astropy.io` (for FITS), `numpy` and `scipy.signal` (for Savizky-Golay filter). Python is only used on the red arm and only used once in order to call `scripts/frTungstenFibreFlat.py`. This is used to generate a per fibre flat field and is only requires because the red arm CCD has a highly complicated and detailed flat field. For a typical CCD (including the blue arm on LT),the call to this function can be completely eliminated from `scripts/L2_exec` and then Python is no longer required.

n.b. for LT operations, `L2_BASE_DIR`, `LIBS` and `INCLUDES` should already have a commented option for `lt-qc`.

# Arc Calibration Solutions

**Unless the instrumental setup has significantly changed and the positions of the arc lines have shifted on the detector dramatically, this procedure should not be required. The pipeline has a degree of built-in robustness to temperature dependent changes in arc line position.**

Arc calibration file lookup tables (`arc.tab`) are kept in `config/lookup_tables` with the subdirectories `blue_g`, `blue_vph`, `red_g` and `red_vph` for each optimised configuration. This file provides a record of: arc calibration file location (relative to `reference_arcs/[red_g||red_vph||blue_g||blue_vph]`), date active from, time active from, date active to and time active to. The last two fields can be replaced by a single value, "now", signifying that the file can be used up to the current date/time, e.g.

`240912/arc.lis	24/09/2012	12:00:00	now`

This table should be appended to for any change in configuration.

The actual configuration files themselves are normally kept in `config/reference_arcs/[red_g||red_vph||blue_g||blue_vph]/[**DATE**]/` but obviously could be kept anywhere, as long as the location is correctly specified in the corresponding lookup table. 

The file itself contains a tab-separated list of identified arc lines as `x	lambda` where [**x**] is in pixels, and [**lambda**] is in Angstrom. The x positions should be as they would be post-trimming of the spectrum. As such, you may need to run the pipeline with the arc as the target and use the intermediate file created (`_target_tr_cor.fits`) to generate a suitable frame to work from. The process of identifiying arc lines is out of the scope of this README, but in short, you will need to collapse the spectrum along the spatial axis, plot it, and cross-match lines visually with that of a known arc lamp spectrum.

# Invoking the Pipeline

The pipeline is invoked specifying the target, arc and flat frame on the command line, in that order, e.g.

`[rmb@rmb-tower scripts]$ ./L2_exec [TARGET_FRAME] [ARC_FRAME] [FLAT_FRAME]`


