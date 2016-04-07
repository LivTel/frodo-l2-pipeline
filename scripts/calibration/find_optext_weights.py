#!/usr/bin/python

from numpy import * 

from pylab import *
import pylab as plt
from scipy import optimize, stats, interpolate
from mpl_toolkits.mplot3d import Axes3D
import subprocess
from FITSFile import *
import time
from errors import *
import optparse
import copy
from lmfit import Model
from scipy.ndimage.filters import median_filter

def gaussian(x, height, center_x, width_x):
    """Returns the value of a gaussian at x with the given parameters"""
    width_x = float(width_x)
    return height*np.exp(-(center_x-x)**2/(2*width_x**2))

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
        if np.count_nonzero(fwhm) > 2:
            fwhm = fwhm[np.nonzero(fwhm)]				# last !None element
            fwhm_predict = fwhm[-1]
            fwhm_predict_lower = fwhm_predict - (3*np.std(fwhm))
            fwhm_predict_upper = fwhm_predict + (3*np.std(fwhm))
        else:
            fwhm_predict = p['minFWHM'] + ((p['maxFWHM'] - p['minFWHM'])/2)
            fwhm_predict_lower = p['minFWHM']
            fwhm_predict_upper = p['maxFWHM']

        gmod.set_param_hint('h1', min=p['minPeakVal'],  max=p['maxPeakVal'])
        gmod.set_param_hint('c1', min=params[1]-0.1, max=params[1]+0.1)
        gmod.set_param_hint('w1', min=fwhm_predict_lower/2.35, max=fwhm_predict_upper/2.35)

    p = gmod.make_params(height=params[0], center_x=params[1], width_x=params[2]) 
    res = gmod.fit(y, x=x, params=p)

    logging.debug("(__fitgaussian__) I predicted " + str(round(fwhm_predict, 2)) + " and got " + str(round(res.best_values['width_x']*2.35, 2)))

    return [res.best_values['height'], res.best_values['center_x'], res.best_values['width_x']], res.chisqr, res.residual, res
    
if __name__ == "__main__":
    '''
    Find optimal extraction weights for a given continuum file.
    '''

    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General options")
    group1.add_option('--l', action='store', default='INFO', dest='logLevel', help='logging level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')
    group1.add_option('--d', action='store_true', dest='pdebug', help='plot debug mode')
    group1.add_option('--f', action='store', default="./flat.fits", dest='pathToFlatFile', help='path to flat file')
    group1.add_option('--of', action='store', default="weights.dat", type=str, dest='outFile', help='output filepath')
    group1.add_option('--s', action='store', dest='stride', default=50, help='stride factor along dispersion axis (y)')
    parser.add_option_group(group1)

    group2 = optparse.OptionGroup(parser, "Peak filtering sensitivity options")
    group2.add_option('--minp', action='store', default=250, type=float, dest='minPeakVal', help='minimum data value for peak moment to be considered valid') 
    group2.add_option('--maxp', action='store', default=2000, type=float, dest='maxPeakVal', help='maximum data value for peak moment to be considered valid')
    group2.add_option('--minf', action='store', default=1.5, type=float, dest='minFWHM', help='minimum peak moment for FWHM to be considered valid') 
    group2.add_option('--maxf', action='store', default=6, type=float, dest='maxFWHM', help='maximum peak moment for FWHM to be considered valid')
    group2.add_option('--maxd', action='store', default=0.1, type=float, dest='maxResidDeviation', help='maximum percentage deviation of average residual from (first) peak value')
    parser.add_option_group(group2)

    group3 = optparse.OptionGroup(parser, "Fitting options") 
    group3.add_option('--a', action='store', default=7, type=int, dest='fittingApertureX', help='fitting aperture x (px)') 
    group3.add_option('--as', action='store', default=7, type=int, dest='fittingScatteredApertureX', help='fitting scattered light aperture x (px)') 
    group3.add_option('--o', action='store', default=4, type=int, dest='fitOrder', help='fitting order')
    parser.add_option_group(group3)

    args = parser.parse_args()
    options, args = parser.parse_args()

    params = {
              'logLevel' : str(options.logLevel), 
              'pdebug' : bool(options.pdebug),
              'pathToFlatFile' : str(options.pathToFlatFile),
              'outFile' : str(options.outFile),
              'stride' : int(options.stride),
              'minPeakVal' : float(options.minPeakVal),
              'maxPeakVal' : float(options.maxPeakVal),
              'minFWHM' : float(options.minFWHM),
              'maxFWHM' : float(options.maxFWHM),
              'maxResidDeviation' : float(options.maxResidDeviation),
              'fittingApertureX' : int(options.fittingApertureX),
              'fittingScatteredApertureX' : int(options.fittingScatteredApertureX),
              'order' : int(options.fitOrder)
    }

    logging.basicConfig(format='%(levelname)s: %(message)s', level=getattr(logging, params['logLevel'].upper()))

    ## **************
    ## INPUT CHECKS *
    ## **************
    ##
    ## a) fittingApertureX must be odd
    if params['fittingApertureX'] % 2 == 0:	
        im_err._setError(-5)
        im_err.handleError()
    ## b) fittingScatteredApertureX must be odd
    if params['fittingApertureX'] % 2 == 0:	
        im_err._setError(-5)
        im_err.handleError()

    ## ************************************
    ## GET ENVIRONMENT VARS AND SET PATHS *
    ## ************************************
    ##
    L2_bin_dir_path = os.getenv("L2_BIN_DIR")
    L2_bin_frfind  = L2_bin_dir_path + "/frfind"
    L2_bin_frclean = L2_bin_dir_path + "/frclean"
    L2_bin_frtrace = L2_bin_dir_path + "/frtrace"

    ## ****************
    ## OPEN FITS FILE *
    ## ****************
    ##
    im_err = errors()
    flatFile_nobin = FITSFile(params['pathToFlatFile'], im_err)
    if flatFile_nobin.openFITSFile():
        if not flatFile_nobin.getHeaders(0):
            im_err._setError(-2)
            im_err.handleError() 
        if not flatFile_nobin.getData(0):
            im_err._setError(-3)
            im_err.handleError()          
    else:
        im_err._setError(-1)
        im_err.handleError()

    ## ***************************
    ## PROCESS FITS FILE WITH L2 *
    ## ***************************
    ## frfind
    ##
    '''logging.info("(__main__) executing frfind") 
    proc = subprocess.Popen([L2_bin_frfind, params['pathToFlatFile'], "3", "3", "10", "30", "2", "1000"], stdout=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        print(line.rstrip())
    proc.wait()

    ## frclean
    ##
    proc = subprocess.Popen([L2_bin_frclean, "6", "30"], stdout=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        print(line.rstrip())
    proc.wait()

    ## frtrace
    ##
    proc = subprocess.Popen([L2_bin_frtrace, "2", "100", "5"], stdout=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        print(line.rstrip())
    proc.wait()'''

    ## *********************
    ## TAKE STRIDE OF DATA *
    ## ********************* 
    ##    
    data_strided_y = []
    data_strided_z = []
    logging.info("(__main__) taking stride of data (" + str(params['stride']) + "px) along dispersion axis") 
    for idx_r in range(0, len(flatFile_nobin.data), params['stride']):
        data_strided_y.append(idx_r)
        data_strided_z.append(flatFile_nobin.data[idx_r])

    ## *****************************
    ## PARSE L2 OUTPUT AND EXTRACT *
    ## *****************************
    ## Parse data from frtrace_traces.dat output file and extract data from FITS file.
    ##
    ## We create 5 arrays, 3 holding the x, y and intensity (z) data for the strided 
    ## data, and a further two holding the x and y centroids.
    ##
    ## Strided data output format: [row][fibre][pixel].
    ## Centroid data output format: [row][fibre].
    ##
    logging.info("(__main__) parsing frtrace output") 
    c_traces = []
    with open('frtrace_traces.dat', 'r') as f:
        for line in f:
            if not line.startswith('#') and line != "\n" and line != "-1":
                tmp = line.strip('\t\n').split('\t')
                tmp = [float(i) for i in tmp[::-1]]
                c_traces.append(tmp[1:-1])

    data_strided_x_aperture = []
    data_strided_y_aperture = []
    data_strided_z_aperture = []
    data_strided_x_centroid = []
    data_strided_y_centroid = []
    for idx_r in data_strided_y:
        this_row_data_strided_x_aperture = []
        this_row_data_strided_y_aperture = []
        this_row_data_strided_z_aperture = []
        this_row_data_strided_x_centroid = []
        this_row_data_strided_y_centroid = []
        for idx_f in range(0,144):
            this_c = c_traces[idx_f]
            this_centroid = np.polyval(this_c, idx_r)-1		# polynomial from file is 1-indexed
            this_row_fibre_data_strided_x_aperture = []
            this_row_fibre_data_strided_y_aperture = []
            this_row_fibre_data_strided_z_aperture = []
            for xi in range(int(round(this_centroid,2))-((params['fittingApertureX']-1)/2), int(round(this_centroid,2))+((params['fittingApertureX']+1)/2)):    # +1 as range is not inclusive of maximum value
                this_row_fibre_data_strided_x_aperture.append(xi)
                this_row_fibre_data_strided_y_aperture.append(idx_r)
                this_row_fibre_data_strided_z_aperture.append(flatFile_nobin.data[idx_r, xi])
            this_row_data_strided_x_aperture.append(this_row_fibre_data_strided_x_aperture)
            this_row_data_strided_y_aperture.append(this_row_fibre_data_strided_y_aperture)
            this_row_data_strided_z_aperture.append(this_row_fibre_data_strided_z_aperture)
            this_row_data_strided_x_centroid.append(this_centroid)
            this_row_data_strided_y_centroid.append(idx_r)  
        data_strided_x_aperture.append(this_row_data_strided_x_aperture)
        data_strided_y_aperture.append(this_row_data_strided_y_aperture)
        data_strided_z_aperture.append(this_row_data_strided_z_aperture)
        data_strided_x_centroid.append(this_row_data_strided_x_centroid)
        data_strided_y_centroid.append(this_row_data_strided_y_centroid)

    ## ************************************
    ## SCATTERED LIGHT CENTROIDS AND DATA *
    ## ************************************
    ##
    ## First we find the interfibre centroids by taking the middlepoint between two adjacent fibre peaks
    ##
    data_strided_x_centroid_scattered = []
    data_strided_y_centroid_scattered = []
    for row_x, row_y in zip(data_strided_x_centroid, data_strided_y_centroid):
        data_strided_x_centroid_scattered.append([i+((j-i)/2) for i, j in zip(row_x[0:-1], row_x[1:])])
        data_strided_y_centroid_scattered.append([i+((j-i)/2) for i, j in zip(row_y[0:-1], row_y[1:])])

    ##
    ## Then, for each centroid, we perform a cubic interpolation 
    ##
    data_strided_z_centroid_scattered = []
    for idx_r, row in enumerate(data_strided_y_centroid_scattered):
        this_row_data_strided_z_centroid_scattered = []
        for idx_c, centroid_scattered in enumerate(row):
            this_centroid_x_scattered = data_strided_x_centroid_scattered[idx_r][idx_c]
            this_row_interfibre_data_strided_x_centroid_scattered_aperture = []
            this_row_interfibre_data_strided_z_centroid_scattered_aperture = []
            for xi in range(int(round(this_centroid_x_scattered,2))-((params['fittingScatteredApertureX']-1)/2), int(round(this_centroid_x_scattered,2))+((params['fittingScatteredApertureX']+1)/2)):
                this_row_interfibre_data_strided_x_centroid_scattered_aperture.append(xi)
                this_row_interfibre_data_strided_z_centroid_scattered_aperture.append(flatFile_nobin.data[idx_r, xi])

            this_min_idx = np.argmin(this_row_interfibre_data_strided_z_centroid_scattered_aperture)	# find minimum index (z)
            this_centroid_range = [this_min_idx-1:this_min_idx+2]
            this_c = np.polyfit(this_row_interfibre_data_strided_x_centroid_scattered_aperture[], this_row_interfibre_data_strided_z_centroid_scattered_aperture[this_centroid_range], 2)	
            print this_c  
            this_row_data_strided_z_centroid_scattered.append(np.min(np.polyval(this_valley_c, this_row_interfibre_data_strided_x_centroid_scattered_aperture[this_min_idx-1:this_min_idx+2])))						# find minimum
        data_strided_z_centroid_scattered.append(this_row_data_strided_z_centroid_scattered)

    ## *************************
    ## SCATTERED LIGHT REMOVAL *
    ## *************************
    ##   
    data_strided_x_centroid_scattered  = np.asarray(data_strided_x_centroid_scattered)
    data_strided_y_centroid_scattered  = np.asarray(data_strided_y_centroid_scattered)
    data_strided_z_centroid_scattered  = np.asarray(data_strided_z_centroid_scattered)
    # smooth along spatial
    z_spatial = []
    for idx_r, row in enumerate(data_strided_y_centroid_scattered):  
        #f = interpolate.UnivariateSpline(data_strided_x_centroid_scattered[idx_r], data_strided_z_centroid_scattered[idx_r], k=3)
        f = np.polyfit(data_strided_x_centroid_scattered[idx_r], data_strided_z_centroid_scattered[idx_r], 1)
        #z_spatial.append(f(data_strided_x_centroid_scattered[idx_r]))
        z_spatial.append(np.polyval(f, data_strided_x_centroid_scattered[idx_r]))
        '''plt.plot(data_strided_x_centroid_scattered[idx_r], data_strided_z_centroid_scattered[idx_r])
        plt.plot(data_strided_x_centroid_scattered[idx_r], np.polyval(f, data_strided_x_centroid_scattered[idx_r]))
        plt.show()'''
    z_spatial = np.asarray(z_spatial)

    # smooth along dispersion
    z_disp = [] 
    for idx_f, fib in enumerate(z_spatial.transpose()):   
        #f = interpolate.UnivariateSpline(data_strided_x_centroid_scattered.transpose()[idx_f], z_spatial.transpose()[idx_f], k=1)
        f = np.polyfit(data_strided_x_centroid_scattered.transpose()[idx_f], z_spatial.transpose()[idx_f], 3)
        #z_disp.append(f(data_strided_x_centroid_scattered.transpose()[idx_f]))
        z_disp.append(np.polyval(f, data_strided_x_centroid_scattered.transpose()[idx_f]))
        '''plt.plot(data_strided_x_centroid_scattered.transpose()[idx_f], z_spatial.transpose()[idx_f])
        plt.plot(data_strided_x_centroid_scattered.transpose()[idx_f], np.polyval(f, data_strided_x_centroid_scattered.transpose()[idx_f]))
        plt.show()'''
    z_disp = np.asarray(z_disp).transpose()

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(data_strided_x_centroid_scattered, data_strided_y_centroid_scattered, data_strided_z_centroid_scattered, cmap="Reds", alpha=0.5)
    ax.plot_surface(data_strided_x_centroid_scattered, data_strided_y_centroid_scattered, z_disp)
    plt.show()
    exit(0)
    
    ## *************
    ## PSF FITTING *
    ## *************
    ## i) calculate moments of y data
    ## ii) filter moments for FWHM and height
    ## iii) perform fit with boundaries
    ## iv) calculate a residual
    ##
    logging.info("(__main__) calculating fibre PSF") 
    residuals = []
    fwhm1 = []
    for idx_r, row in enumerate(data_strided_x_aperture):
        this_fibre_residual = []
        this_row_fwhm1 = []
        for idx_f, fibre in enumerate(row):
            logging.info("(__main__) processing row " + str(idx_r+1) + " fibre " + str(idx_f+1)) 
            this_row_x_values = []
            this_row_z_values = []
            for idx_p, pixel in enumerate(fibre):
                this_row_x_values.append(data_strided_x_aperture[idx_r][idx_f][idx_p])
                this_row_z_values.append(data_strided_z_aperture[idx_r][idx_f][idx_p])
 
            # calculate moments
            momParams = moments(np.asarray(this_row_x_values), np.asarray(this_row_z_values))
            momHeight = momParams[0]
            momX      = momParams[1]
            momFWHMX  = 2.35*momParams[2]
 
            # filter on moments
            if momHeight < params['minPeakVal'] or momHeight > params['maxPeakVal']:
                logging.info("(__main__) moment height parameter is not within range (" + str(momHeight) + ")")
                this_row_fwhm1.append(None)
                continue

            # fit
            if idx_r==0:
                this_fibre_fwhm = np.asarray([])
            else:
                this_fibre_fwhm = np.asarray(fwhm1).transpose()[:][idx_f]
            fitParams, fitChiSPDOF, fitResid, fitRes = fitgaussian(np.asarray(this_row_x_values), np.asarray(this_row_z_values), params, this_fibre_fwhm)
            fitHeight1	= fitParams[0]
            fitX1 	= fitParams[1]
            fitFWHMX1 	= 2.35*fitParams[2]

            logging.debug("(__main__) parameters:\t" + str(round(fitParams[0], 2)) + '\t' + str(round(fitParams[1], 2)) + '\t' + str(round(fitParams[2], 2)))

            this_percentage_fitRes = pow(np.mean(pow(fitResid, 2)), 0.5)/fitHeight1
            logging.debug("(__main__) residual:\t" + str(this_percentage_fitRes))
            if this_percentage_fitRes > params['maxResidDeviation']:
                logging.info("(__main__) residual deviation is too large (" + str(this_percentage_fitRes) + ")")
                this_row_fwhm1.append(None)
                continue 

            this_row_fwhm1.append(fitFWHMX1)

            if params['pdebug']:
                plt.subplot(211)
                plt.plot(this_row_x_values, this_row_y_values, 'ko-')
                plt.plot(this_row_x_values, gaussian(np.asarray(this_row_x_values), *fitParams[0:3]), 'r--')		# first
                plt.subplot(212)
                plt.plot(this_row_x_values, fitResid, 'k-')								# residual
                plt.show()

            this_fibre_residual.append(fitResid)
        residuals.append(this_fibre_residual)
        fwhm1.append(this_row_fwhm1)

    ## ***********
    ## FIT FWHM1 *
    ## ***********
    ##
    logging.info("(__main__) fitting per fibre polynomial to fibre PSF") 
    fit_FWHM = []	# per fibre
    fit_chiS = []	# per fibre
    for idx_f in range(0, 144):
        this_fibre_data_strided_fwhm1 = []
        this_fibre_data_strided_y = []
        for idx_r in range(0, len(data_strided_y_centroid)):
            if fwhm1[idx_r][idx_f] is not None:
                this_fibre_data_strided_fwhm1.append(fwhm1[idx_r][idx_f])
                this_fibre_data_strided_y.append(data_strided_y_centroid[idx_r][idx_f])
        this_fibre_data_strided_y = np.asarray(this_fibre_data_strided_y)
        this_fibre_data_fwhm1 = np.asarray(this_fibre_data_strided_fwhm1)

        this_fibre_data_strided_fwhm1_fit_c = np.polyfit(this_fibre_data_strided_y, this_fibre_data_strided_fwhm1, params['order']) 
        this_fibre_data_strided_fwhm1_chi_squared = stats.chisquare(this_fibre_data_strided_fwhm1, np.polyval(this_fibre_data_strided_fwhm1_fit_c, this_fibre_data_strided_y))[0]

        fit_FWHM.append(this_fibre_data_strided_fwhm1_fit_c)
        fit_chiS.append(this_fibre_data_strided_fwhm1_chi_squared)
    
    ## *********************
    ## OUTPUT FITS TO FILE *
    ## *********************
    ##
    with open(params['outFile'], 'w') as out:
        out.write("# Polynomial Order:\t" + str(params['order']) + '\n\n')   
        for idx_f in range(0, 144):
            this_fibre_data_strided_fwhm1_fit_c = fit_FWHM[idx_f]
            this_fibre_data_strided_fwhm1_chi_squared = fit_chiS[idx_f]
            out.write(str(idx_f+1) + '\t' + '\t'.join(format(x, "1.10e") for x in this_fibre_data_strided_fwhm1_fit_c[::-1]) + '\t' + str(format(this_fibre_data_strided_fwhm1_chi_squared, "1.2f")) + '\n')

    ## *************
    ## TRACE PLOTS *
    ## *************
    ##
    #  a) centroids
    for idx_f, i in enumerate(np.asarray(data_strided_x_centroid).transpose()):
        x = np.asarray(data_strided_x_centroid).transpose()[idx_f]
        y = np.asarray(data_strided_y_centroid).transpose()[idx_f]
        plt.plot(x, y, 'rx-')
    #  b) 5px aperture tramline
        x_m2point5 = x-2.5
        plt.plot(x_m2point5, y, 'r--')
        x_p2point5 = x+2.5
        plt.plot(x_p2point5, y, 'r--')
    #  c) 1 FWHM tramlines
        x_fit_lower_fwhm = x-np.polyval(fit_FWHM[idx_f], y)
        plt.plot(x_fit_lower_fwhm, y, 'b--')
        x_fit_upper_fwhm = x+np.polyval(fit_FWHM[idx_f], y)
        plt.plot(x_fit_upper_fwhm, y, 'b--')
    #  d) scattered light centroids
    for idx_if, i in enumerate(np.asarray(data_strided_x_centroid_scattered).transpose()):
        x = np.asarray(data_strided_x_centroid_scattered).transpose()[idx_if]
        y = np.asarray(data_strided_y_centroid_scattered).transpose()[idx_if]
        plt.plot(x, y, 'g-')
    #  e) image data
    plt.imshow(flatFile_nobin.data, cmap="Greys", aspect='auto', vmax=np.percentile(flatFile_nobin.data, 99.5), vmin=np.percentile(flatFile_nobin.data, 0.5))
    plt.colorbar()
    plt.show()

    #  f) fwhm as a function of dispersion for each fibre
    for idx_f, i in enumerate(np.asarray(data_strided_x_centroid).transpose()):
        x = np.asarray(data_strided_y_centroid).transpose()[idx_f]
        y = np.polyval(fit_FWHM[idx_f], x)
        plt.subplot(12,12,idx_f+1)
        plt.plot(x, y)
    plt.show()


    # cleanup
    '''os.remove("binned.fits")
    os.remove("frfind_peaks.dat")
    os.remove("frclean_cleaned.dat")
    os.remove("frtrace_traces.dat")
    os.remove("additional_keys")
    os.remove("error_codes")
    flatFile_nobin.closeFITSFile()'''

