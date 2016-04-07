#!/usr/bin/python

from numpy import * 

from pylab import *
import pylab as plt
from scipy import optimize, stats, interpolate
import subprocess
from FITSFile import *
import time
from errors import *
import optparse
import copy
from lmfit import Model

def gaussian(x, height, center_x, width_x):
    """Returns the value of a gaussian at x with the given parameters"""
    width_x = float(width_x)
    return height*np.exp(-(center_x-x)**2/(2*width_x**2))

def two_gaussians(x, h1, c1, w1, h2, c2, w2):
    return gaussian(x, h1, c1, w1) + gaussian(x, h2, c2, w2)

def moments(x, y) :
    """Returns (height, x, width_x)
    the gaussian parameters of a 2D distribution by calculating its
    moments"""
    mu 	= sum(x*y)/sum(y)                         # first moment
    var   = abs(sum((x-mu)**2*y)/sum(y))          # second moment
    sigma = sqrt(var)
    max   = y.max()
    return max, mu, sigma

def fitgaussian(x, y, p, fwhm, useBounds=True):
    """Returns (height, x, width_x)
    the gaussian parameters of a 2D distribution found by leastsq"""
    params = moments(x, y)

    gmod = Model(gaussian)
    if useBounds: 
        # use last FWHM
        if len(np.where(fwhm != None)) > 0:
            fwhm_predict = fwhm[-1]
            fwhm_predict_lower = fwhm_predict - (3*np.std(fwhm))
            fwhm_predict_upper = fwhm_predict + (3*np.std(fwhm))
        else:
            fwhm_predict = p['minFWHM'] + ((p['maxFWHM'] - p['minFWHM'])/2)
            fwhm_predict_lower = p['minFWHM']
            fwhm_predict_upper = p['maxFWHM']

        gmod.set_param_hint('h1', min=p['minPeakVal'],  max=p['maxPeakVal'])
        gmod.set_param_hint('c1', min=params[1]-0.1, max=params[1]+0.1)				# make sure two gaussians are coincident within 0.1px of the centroid moment
        gmod.set_param_hint('w1', min=fwhm_predict_lower/2.35, max=fwhm_predict_upper/2.35)

    p = gmod.make_params(height=params[0], center_x=params[1], width_x=params[2]) 
    res = gmod.fit(y, x=x, params=p)

    logging.debug("(__fitgaussian__) I predicted " + str(round(fwhm_predict, 2)) + " and got " + str(round(res.best_values['width_x']*2.35, 2)))

    return [res.best_values['height'], res.best_values['center_x'], res.best_values['width_x']], res.chisqr, res.residual, res

def fit2gaussian(x, y, p, fwhm1, fwhm2, useBounds=True):
    """Returns (height1, x1, width1, height2, x2, width2)
    the gaussian parameters of a 2D distribution found by leastsq"""
    params1 = moments(x, y)
    params2 = (params1[0]/10, params1[1]+0.1, params1[2]*3)				# approx. init factors based off culmulative moments, empirically obtained

    params = params1 + params2

    gmod = Model(two_gaussians)
    if useBounds:
        # the following uses an interpolation of the current FWHM for each fibre to predict the next FWHM
        if len(fwhm1) > 2 and len(fwhm2) > 2:
            c1 = np.polyfit(range(len(fwhm1)), fwhm1, 2)
            fwhm1_predict = np.polyval(c1, len(fwhm1))
            fwhm1_predict_lower = fwhm1_predict - (3*np.std(fwhm1))
            fwhm1_predict_upper = fwhm1_predict + (3*np.std(fwhm1))
            c2 = np.polyfit(range(len(fwhm2)), fwhm2, 2)
            fwhm2_predict = np.polyval(c2, len(fwhm2))
            fwhm2_predict_lower = fwhm2_predict - (3*np.std(fwhm2))
            fwhm2_predict_upper = fwhm2_predict + (3*np.std(fwhm2))
        else:
            fwhm1_predict = params[2]
            fwhm1_predict_lower = round(p['minFWHM'], 2)
            fwhm1_predict_upper = round(p['maxFWHM'], 2)
            fwhm2_predict = params[2]
            fwhm2_predict_lower = round(p['minFWHM'], 2)
            fwhm2_predict_upper = round(p['maxFWHM'], 2)

        gmod.set_param_hint('h1', min=p['minPeakVal'],  max=p['maxPeakVal'])
        gmod.set_param_hint('w1', min=fwhm1_predict_lower/2.35, max=fwhm1_predict_upper/2.35)
        gmod.set_param_hint('c1', min=params[1]-0.1, max=params[1]+0.1)			# make sure two gaussians are coincident within 0.1px of the centroid moment

        gmod.set_param_hint('h2', min=p['minPeakVal'], max=p['maxPeakVal'])		# empirically obtained
        gmod.set_param_hint('w2', min=fwhm1_predict_lower/2.35, max=fwhm1_predict_upper/2.35)
        gmod.set_param_hint('c2', min=params[1]-0.1, max=params[1]+0.1)

    p = gmod.make_params(h1=params[0], c1=params[1], w1=params[2], h2=params[3], c2=params[4], w2=params[5]) 
    res = gmod.fit(y, x=x, params=p)

    logging.debug("(__fit2gaussian__) I predicted the first fwhm as " + str(round(fwhm1_predict, 2)) + "+/-" + str(fwhm1_predict_lower) + " and got " + str(round(res.best_values['w1']*2.35, 2)))
    logging.debug("(__fit2gaussian__) I predicted the second fwhm as " + str(round(fwhm2_predict, 2)) + "+/-" + str(fwhm1_predict_lower) + " and got " + str(round(res.best_values['w2']*2.35, 2)))


    return [res.best_values['h1'], res.best_values['c1'], res.best_values['w1'], res.best_values['h2'], res.best_values['c2'], res.best_values['w2']], res.chisqr, res.residual, res
    
if __name__ == "__main__":
    '''
    Find optimal extraction weights for a given continuum file.
    '''

    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General options")
    group1.add_option('--l', action='store', default='INFO', dest='logLevel', help='logging level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')
    group1.add_option('--d', action='store_true', dest='pdebug', help='plot debug mode')
    group1.add_option('--f', action='store', default="./flat.fits", dest='pathToFlatFile', help='path to flat file')
    group1.add_option('--o', action='store', default="weights.dat", type=str, dest='outFile', help='output filepath')
    group1.add_option('--b', action='store', dest='binning', default=50, help='binning factor along dispersion axis (y)')
    group1.add_option('--t', action='store', dest='fitType', default="onegauss", help='fit type (onegauss||twogauss')   
    parser.add_option_group(group1)

    group2 = optparse.OptionGroup(parser, "Peak filtering sensitivity options")
    group2.add_option('--minp', action='store', default=50, type=float, dest='minPeakVal', help='minimum data value for peak moment to be considered valid') 
    group2.add_option('--maxp', action='store', default=2000, type=float, dest='maxPeakVal', help='maximum data value for peak moment to be considered valid')
    group2.add_option('--minf', action='store', default=1.5, type=float, dest='minFWHM', help='minimum peak moment for FWHM to be considered valid') 
    group2.add_option('--maxf', action='store', default=8, type=float, dest='maxFWHM', help='maximum peak moment for FWHM to be considered valid')
    group2.add_option('--maxd', action='store', default=0.1, type=float, dest='maxResidDeviation', help='maximum percentage deviation of average residual from (first) peak value')

    group3 = optparse.OptionGroup(parser, "Fitting options") 
    group3.add_option('--fx', action='store', default=7, type=int, dest='fittingApertureX', help='fitting aperture x (px)') 
    group3.add_option('--or', action='store', default=4, type=int, dest='fitOrder', help='fitting order')
    group3.add_option('--sxw', action='store', default='30,100,1130,1200', type=str, dest='scatteredLightX', help='scattered light window (x1,x2,x3,x4)')

    args = parser.parse_args()
    options, args = parser.parse_args()

    params = {
              'logLevel' : str(options.logLevel), 
              'pdebug' : bool(options.pdebug),
              'pathToFlatFile' : str(options.pathToFlatFile),
              'outFile' : str(options.outFile),
              'binning' : int(options.binning),
              'fitType' : str(options.fitType),
              'minPeakVal' : float(options.minPeakVal),
              'maxPeakVal' : float(options.maxPeakVal),
              'minFWHM' : float(options.minFWHM),
              'maxFWHM' : float(options.maxFWHM),
              'maxResidDeviation' : float(options.maxResidDeviation),
              'fittingApertureX' : int(options.fittingApertureX),
              'order' : int(options.fitOrder),
	      'scatteredLightX' : [int(x) for x in options.scatteredLightX.split(',')]
    }

    logging.basicConfig(format='%(levelname)s: %(message)s', level=getattr(logging, params['logLevel'].upper()))

    ## **********************************
    ## GET ENVIRONMENT VARS AND SET PATHS
    # get L2 environment variables
    L2_bin_dir_path = os.getenv("L2_BIN_DIR")
    # set L2 paths
    L2_bin_frfind  = L2_bin_dir_path + "/frfind"
    L2_bin_frclean = L2_bin_dir_path + "/frclean"

    im_err = errors()
    flatFile = FITSFile(params['pathToFlatFile'], im_err)
    if flatFile.openFITSFile():
        if not flatFile.getHeaders(0):
            im_err._setError(-2)
            im_err.handleError() 
        if not flatFile.getData(0):
            im_err._setError(-3)
            im_err.handleError()          
    else:
        im_err._setError(-1)
        im_err.handleError()

    ## ********
    ## BIN DATA
    ##
    binned_data = []
    for j in range(0, len(flatFile.data)-params['binning'], params['binning']):
        lower = j
        upper = j+params['binning']
        binned_data.append(np.mean(flatFile.data[lower:upper], axis=0))
    pyfits.writeto("binned.fits", binned_data, flatFile.headers)

    flatFile = FITSFile("binned.fits", im_err)
    if flatFile.openFITSFile():
        if not flatFile.getHeaders(0):
            im_err._setError(-2)
            im_err.handleError() 
        if not flatFile.getData(0):
            im_err._setError(-3)
            im_err.handleError()          
    else:
        im_err._setError(-1)
        im_err.handleError()

    ## ************
    ## INPUT CHECKS
    ## some sanity checks
    ##
    if params['fittingApertureX'] % 2 == 0:	# fittingApertureX must be odd
        im_err._setError(-5)
        im_err.handleError()

    ## ***************
    ## PROCESS WITH L2
    ## frfind
    ##
    logging.info("(__main__) executing frfind") 
    proc = subprocess.Popen([L2_bin_frfind, "binned.fits", "3", "3", "10", "30", "2", str(1000/params['binning'])], stdout=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        print(line.rstrip())
    proc.wait()
    proc = subprocess.Popen([L2_bin_frclean, "6", "30"], stdout=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        print(line.rstrip())
    proc.wait()

    ## ************************
    ## PARSE PEAK FINDER OUTPUT
    ##
    ## parse x, y data from frfind_peaks.dat output file into sublisted lists
    ## of format: [row][fibre][pixel], subtracting broad scattered light component.
    ##
    logging.info("(__main__) parsing frfind output") 
    data_x = []
    data_y = []
    this_row_y = None
    with open('frfind_peaks.dat', 'r') as f:
        for line in f:
            if not line.startswith('#') and line != "\n" and line != "-1":
                tmp = line.strip('\t\n').split('\t')
                thisX = int(round(float(tmp[1]))) - 1 # zero indexed
                thisY = int(tmp[2]) - 1               # zero indexed

                if this_row_y == None:
                    this_row_y = thisY
                    this_row_x_values = []
                    this_row_y_values = []
                elif this_row_y != thisY:
                    if len(this_row_x_values) == 144:
                        logging.info("(__main__) found a suitable candidate row at y = " + str(this_row_y)) 
                        data_x.append(this_row_x_values)
                        data_y.append(this_row_y_values)
                    this_row_y = thisY
                    this_row_x_values = []
                    this_row_y_values = []

                # fit broad scattered light component
                sc_x = range(params['scatteredLightX'][0], params['scatteredLightX'][1]) + range(params['scatteredLightX'][2], params['scatteredLightX'][3]) 
                sc_y = flatFile.data[thisY, params['scatteredLightX'][0]:params['scatteredLightX'][1]].tolist() + flatFile.data[thisY, params['scatteredLightX'][2]:params['scatteredLightX'][3]].tolist()
                sc_c = np.polyfit(sc_x, sc_y, 2)

                this_row_fibre_x_data = []
                this_row_fibre_y_data = []
       	        # construct sublist in data containing image data values to fit
                for xi in range(thisX-((params['fittingApertureX']-1)/2), thisX+((params['fittingApertureX']+1)/2)):    # +1 as range is not inclusive of maximum value
                    this_row_fibre_x_data.append(xi)
                    this_row_fibre_y_data.append(flatFile.data[thisY, xi] - np.polyval(sc_c, xi))
                this_row_x_values.append(this_row_fibre_x_data)
                this_row_y_values.append(this_row_fibre_y_data)

    ## FIT FIBRE PSF:
    ## i) remove scattered light contribution
    ## ii) calculate moments of y data
    ## iii) filter moments for FWHM and height
    ## iv) perform fit with boundaries
    ## v) calculate a residual
    ## 
    logging.info("(__main__) fitting fibre PSF with fit type: " + str(params['fitType'])) 
    residuals = []
    fwhm1 = []
    fwhm2 = []
    for idx_r, row in enumerate(data_x):
        this_fibre_residual = []
        this_row_fwhm1 = []
        this_row_fwhm2 = []
        for idx_f, fibre in enumerate(row):
            logging.info("(__main__) processing row " + str(idx_r+1) + " fibre " + str(idx_f+1)) 
            this_row_x_values = []
            this_row_y_values = []
            for idx_p, pixel in enumerate(fibre):
                this_row_x_values.append(data_x[idx_r][idx_f][idx_p])
                this_row_y_values.append(data_y[idx_r][idx_f][idx_p])
 
            # calculate moments
            momParams = moments(np.asarray(this_row_x_values), np.asarray(this_row_y_values))
            momHeight = momParams[0]
            momX      = momParams[1] + thisX
            momFWHMX  = 2.35*momParams[2]

            # filter on moments
            if momHeight < params['minPeakVal'] or momHeight > params['maxPeakVal']:
                logging.info("(__main__) moment height parameter is not within range (" + str(momHeight) + ")")
                this_row_fwhm1.append(None)
                this_row_fwhm2.append(None)
                continue
            elif momFWHMX < params['minFWHM'] or momFWHMX > params['maxFWHM']:
                logging.info("(__main__) moment fwhm parameter is not within range (" + str(momFWHMX) + ")")
                this_row_fwhm1.append(None)
                this_row_fwhm2.append(None)
                continue

            # fit
            if params['fitType'] == "onegauss":
                if idx_r==0:
                    this_fibre_fwhm = np.asarray([None])
                else:
                    this_fibre_fwhm = np.asarray(fwhm1).transpose()[:][idx_f]
                fitParams, fitChiSPDOF, fitResid, fitRes = fitgaussian(np.asarray(this_row_x_values), np.asarray(this_row_y_values), params, this_fibre_fwhm)
                fitHeight1	= fitParams[0]
                fitX1 		= fitParams[1] + thisX
                fitFWHMX1 	= 2.35*fitParams[2]
                fitHeight2	= None
                fitX2 		= None
                fitFWHMX2 	= None

                logging.debug("(__main__) parameters:\t" + str(round(fitParams[0], 2)) + '\t' + str(round(fitParams[1], 2)) + '\t' + str(round(fitParams[2], 2)))

            elif params['fitType'] == "twogauss":
                if idx_r==0:
                    this_fibre_fwhm1 = np.asarray([None])
                    this_fibre_fwhm2 = np.asarray([None])
                else:
                    this_fibre_fwhm1 = np.asarray(fwhm1).transpose()[:][idx_f]
                    this_fibre_fwhm2 = np.asarray(fwhm2).transpose()[:][idx_f]
                fitParams, fitChiSPDOF, fitResid, fitRes = fit2gaussian(np.asarray(this_row_x_values), np.asarray(this_row_y_values), params, this_fibre_fwhm1, this_fibre_fwhm2)
                fitHeight1	= fitParams[0]
                fitX1 		= fitParams[1] + thisX
                fitFWHMX1 	= 2.35*fitParams[2]
                fitHeight2	= fitParams[3]
                fitX2 		= fitParams[4] + thisX
                fitFWHMX2 	= 2.35*fitParams[5]

                logging.debug("(__main__) parameters (1):\t" + str(round(fitParams[0], 2)) + '\t' + str(round(fitParams[1], 2)) + '\t' + str(round(fitParams[2], 2)))
                logging.debug("(__main__) parameters (2):\t" + str(round(fitParams[3], 2)) + '\t' + str(round(fitParams[4], 2)) + '\t' + str(round(fitParams[5], 2)))

            this_percentage_fitRes = pow(np.mean(pow(fitResid, 2)), 0.5)/fitHeight1
            logging.debug("RESIDUAL:\t" + str(this_percentage_fitRes))
            if this_percentage_fitRes > params['maxResidDeviation']:
                logging.info("(__main__) residual deviation is too large (" + str(this_percentage_fitRes) + ")")
                this_row_fwhm1.append(fitFWHMX1)
                this_row_fwhm2.append(fitFWHMX2)
                continue 

            this_row_fwhm1.append(fitFWHMX1)
            this_row_fwhm2.append(fitFWHMX2)

            if params['pdebug']:
                plt.subplot(211)
                plt.plot(this_row_x_values, this_row_y_values, 'ko-')
                plt.plot(this_row_x_values, gaussian(np.asarray(this_row_x_values), *fitParams[0:3]), 'r--')		# first
                if params['fitType'] == "twogauss":
                    plt.plot(this_row_x_values, gaussian(np.asarray(this_row_x_values), *fitParams[3:6]), 'b--')	# second
                    plt.plot(this_row_x_values, two_gaussians(np.asarray(this_row_x_values), *fitParams), 'kx-')	# combined
                plt.subplot(212)
                plt.plot(this_row_x_values, fitResid, 'k-')								# residual
                plt.show()

            this_fibre_residual.append(fitResid)
        residuals.append(this_fibre_residual)
        fwhm1.append(this_row_fwhm1)
        fwhm2.append(this_row_fwhm2)

    n_fibres = idx_f+1
    n_rows = idx_r+1

    with open(params['outFile'], 'w') as out:
        fig = plt.figure()
        out.write("# Polynomial Order:\t" + str(params['order']) + '\n\n')   
        for idx_f in range(0, n_fibres):
            this_fibre_data_fwhm1 = []
            y = []
            for idx_r in range(0, n_rows):
                this_fibre_data_fwhm1.append(fwhm1[idx_r][idx_f])
                y.append(idx_r*params['binning'])
            y = np.asarray(y)
            this_fibre_data_fwhm1 = np.asarray(this_fibre_data_fwhm1)

            idx_fwhm1_finite = np.isfinite(y) & np.isfinite(this_fibre_data_fwhm1)
            this_fibre_data_fwhm1_fit_c = np.polyfit(y[idx_fwhm1_finite], this_fibre_data_fwhm1[idx_fwhm1_finite], params['order']) 
            this_fibre_data_fwhm1_chi_squared = stats.chisquare(this_fibre_data_fwhm1[idx_fwhm1_finite], np.polyval(this_fibre_data_fwhm1_fit_c, y[idx_fwhm1_finite]))[0]

            out.write(str(idx_f+1) + '\t' + '\t'.join(format(x, "1.10e") for x in this_fibre_data_fwhm1_fit_c[::-1]) + '\t' + str(format(this_fibre_data_fwhm1_chi_squared, "1.2f")) + '\n')

            plt.subplot(12,12,idx_f+1)
            plt.plot(y, this_fibre_data_fwhm1)
            plt.plot(y, np.polyval(this_fibre_data_fwhm1_fit_c, y))

    if params['fitType'] == "twogauss":
        with open(params['outFile'], 'w') as out:
            fig = plt.figure()
            out.write("# Polynomial Order:\t" + str(params['order']) + '\n\n')   
            for idx_f in range(0, n_fibres):
                this_fibre_data_fwhm2 = []
                y = []
                for idx_r in range(0, n_rows):
                    this_fibre_data_fwhm2.append(fwhm2[idx_r][idx_f])
                    y.append(idx_r*params['binning'])
                y = np.asarray(y)
                this_fibre_data_fwhm2 = np.asarray(this_fibre_data_fwhm2)

                idx_fwhm2_finite = np.isfinite(y) & np.isfinite(this_fibre_data_fwhm2)
                this_fibre_data_fwhm2_fit_c = np.polyfit(y[idx_fwhm2_finite], this_fibre_data_fwhm2[idx_fwhm2_finite], params['order']) 
                this_fibre_data_fwhm2_chi_squared = stats.chisquare(this_fibre_data_fwhm2[idx_fwhm2_finite], np.polyval(this_fibre_data_fwhm2_fit_c, y[idx_fwhm2_finite]))[0]

                out.write(str(idx_f+1) + '\t' + '\t'.join(format(x, "1.10e") for x in this_fibre_data_fwhm2_fit_c[::-1]) + '\t' + str(format(this_fibre_data_fwhm2_chi_squared, "1.2f")) + '\n')
    
                plt.subplot(12,12,idx_f+1)
                plt.plot(y, this_fibre_data_fwhm2)
                plt.plot(y, np.polyval(this_fibre_data_fwhm2_fit_c, y))

    plt.show()

    os.remove("binned.fits")
    os.remove("frfind_peaks.dat")
    os.remove("frclean_cleaned.dat")
    os.remove("additional_keys")
    os.remove("error_codes")
    flatFile.closeFITSFile()

