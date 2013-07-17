#!/usr/bin/python

from numpy import * 
from pylab import *
import pylab
from scipy import optimize
import subprocess
from FITSFile import *
import time
from errors import *
import optparse
from plots import *

def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments"""
    total = data.sum()
    X, Y = indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: ravel(gaussian(*p)(*indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params, maxfev=1000)
    return p

if __name__ == "__main__":
    '''
    Produce focus descriptors for single FRODOSpec arc file.
    '''

    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General options")
    group1.add_option('--loglevel', action='store', default='INFO', dest='logLevel', help='logging level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')
    group1.add_option('--interactive', action='store_true', default=False, dest='interactive', help='show plots?')
    group1.add_option('--arcfile', action='store', default="./r_a_20130314_2_1_1_1.fits", dest='pathToArcFile', help='path to arc file')
    group1.add_option('--outfile', action='store', default="output", type=str, dest='outFile', help='output filepath for fitting results') 
    group1.add_option('--outimage', action='store', default="output", type=str, dest='outImage', help='output plot file') 
    group1.add_option('--dispccdsize', action='store', default=4096, type=int, dest='dispCCDSize', help='size of CCD along dispersion axis (px)') 
    group1.add_option('--spatccdsize', action='store', default=2048, type=int, dest='spatCCDSize', help='size of CCD along spatial axis (px)') 
    parser.add_option_group(group1)

    group2 = optparse.OptionGroup(parser, "Peak fitlering sensitivity options", "[minpeakval], [dertol] are useful options to alter the number of peaks detected. ~4000 is a good number")
    group2.add_option('--dertol', action='store', default=250, type=str, dest='peakfind_derTol', help='derivative tolerance for peakfinding routine') 
    group2.add_option('--maxiter', action='store', default=35, type=int, dest='maxNumIterCycles', help='maximum number of cycles for filtering')     
    group2.add_option('--itersigma', action='store', default=3, type=float, dest='iterSigmaCut', help='sigma FWHM for iterative cut') 
    group2.add_option('--minpeakval', action='store', default=500, type=float, dest='minPeakVal', help='minimum data value for fitted peak to be considered valid') 
    group2.add_option('--maxpeakval', action='store', default=40000, type=float, dest='maxPeakVal', help='maximum data value for fitted peak to be considered valid')
    group2.add_option('--minfwhm', action='store', default=1.5, type=float, dest='minFWHM', help='minimum FWHM for fitted peak to be considered valid') 
    group2.add_option('--maxfwhm', action='store', default=8, type=float, dest='maxFWHM', help='maximum FWHM for fitted peak to be considered valid')
    group2.add_option('--dupdist', action='store', default=1, type=float, dest='duplicateMinDist', help='distance between consecutive peaks where peak is considered duplicate (px)') 
    parser.add_option_group(group2)
    
    group3 = optparse.OptionGroup(parser, "Fitting options")
    group3.add_option('--xstep', action='store', default=50, type=int, dest='xStep', help='number of pixels along spatial axis to skip when fitting best fit line') 
    group3.add_option('--ystep', action='store', default=50, type=int, dest='yStep', help='number of pixels along dispersion axis to skip when fitting best fit line') 
    group3.add_option('--xapfit', action='store', default=7, type=int, dest='fittingApertureX', help='fitting aperture x (px)') 
    group3.add_option('--yapfit', action='store', default=7, type=int, dest='fittingApertureY', help='fitting aperture y (px)')
    parser.add_option_group(group3)

    args = parser.parse_args()
    options, args = parser.parse_args()

    params = {
              'logLevel' : str(options.logLevel),
              'interactive' : bool(options.interactive),
              'peakfind_derTol' : str(options.peakfind_derTol),
              'fittingApertureX' : int(options.fittingApertureX),
              'fittingApertureY' : int(options.fittingApertureY), 
              'duplicateMinDist' : float(options.duplicateMinDist),
              'pathToArcFile' : str(options.pathToArcFile),
              'maxNumIterCycles' : int(options.maxNumIterCycles),
              'iterSigmaCut' : float(options.iterSigmaCut),
              'minPeakVal' : float(options.minPeakVal),
              'maxPeakVal' : float(options.maxPeakVal),
              'minFWHM' : float(options.minFWHM),
              'maxFWHM' : float(options.maxFWHM),
              'xStep' : int(options.xStep),
              'yStep' : int(options.yStep),
              'dispCCDSize' : int(options.dispCCDSize),
              'spatCCDSize' : int(options.spatCCDSize),
              'outFile' : str(options.outFile),
              'outImage' : str(options.outImage)
    }

    ## INPUT CHECKS
    ## some sanity checks
    if params['fittingApertureX'] % 2 == 0:	# fittingApertureX must be odd
        im_err._setError(-5)
        im_err.handleError()

    if params['fittingApertureY'] % 2 == 0:	# fittingApertureY must be odd
        im_err._setError(-6)
        im_err.handleError()

    logging.basicConfig(format='%(levelname)s: %(message)s', level=getattr(logging, params['logLevel'].upper()))

    ## GET ENVIRONMENT VARS AND SET PATHS
    # get L2 environment variables
    L2_bin_dir_path = os.getenv("L2_BIN_DIR")
    # set L2 paths
    L2_bin_frpeakfinder = L2_bin_dir_path + "/frpeakfinder"

    im_err = errors()
    arcFile = FITSFile(params['pathToArcFile'], im_err)
    if arcFile.openFITSFile():
        if not arcFile.getHeaders(0):
            im_err._setError(-2)
            im_err.handleError() 
        if not arcFile.getData(0):
            im_err._setError(-3)
            im_err.handleError()          
    else:
        im_err._setError(-1)
        im_err.handleError()

    ## L2
    ## process with frpeakfinder
    logging.info("(__main__) executing frpeakfinder") 
    proc = subprocess.Popen([L2_bin_frpeakfinder, params['pathToArcFile'], "3", "3", params['peakfind_derTol'], "1", "5", "25", "144"], stdout=subprocess.PIPE)
    for line in iter(proc.stdout.readline,''):
        print(line.rstrip())
    proc.wait()

    ## PARSE AND FILTER MOMENTS ON PEAK FINDER OUTPUT
    ## parse frpeakfinder_peaks output file into sublisted list
    ## of format: [[x1, y1, val1], [x2, y2, val2], ...,[xn, yn, valn]]
    ## fit 2D moments and filter for:
    ## i) height
    ## ii) FWHM
    logging.info("(__main__) parsing and filtering frpeakfinder output") 
    parsedOutput = [] 
    with open('frpeakfinder_peaks.dat', 'r') as f:
        for line in f:
            if not line.startswith('#') and line != "\n" and line != "-1":
                tmp = line.strip('\t\n').split('\t')
                thisX = int(round(float(tmp[1]))) - 1 # zero indexed
                thisY = int(tmp[2]) - 1               # zero indexed

                data = []
       	        # construct sublist in data containing image data values to fit (essentially 2D list "array" containing data values)
                for i in range(thisX-((params['fittingApertureX']-1)/2), thisX+((params['fittingApertureX']+1)/2)):     # +1 as range is not inclusive of maximum value
                    data.append([])
                    for j in range(thisY-((params['fittingApertureY']-1)/2), thisY+((params['fittingApertureY']+1)/2)): # +1 as range is not inclusive of maximum value
                        data[len(data)-1].append(arcFile.data[j, i])

                dataNumpy = np.asarray(data)    # convert nested lists to numpy array
                outParams = moments(dataNumpy) 	# calculate rough estimate of parameters

                fitHeight = outParams[0]
                fitX = outParams[1] + thisX
                fitY = outParams[2] + thisY
                fitFWHMX = 2.35*outParams[3]
                fitFWHMY = 2.35*outParams[4]

                if fitHeight < params['minPeakVal'] or fitHeight > params['maxPeakVal']:
                    continue
                elif fitFWHMX < params['minFWHM'] or fitFWHMX > params['maxFWHM']:
                    continue
                elif fitFWHMY < params['minFWHM'] or fitFWHMY > params['maxFWHM']:
                    continue
 
                # append a sublist, and then append x, y and the value into this sublist
                parsedOutput.append([])
                parsedOutput[len(parsedOutput)-1].append(thisX)
                parsedOutput[len(parsedOutput)-1].append(thisY)
                parsedOutput[len(parsedOutput)-1].append(arcFile.data[thisY, thisX])

    ## FIT PEAKS
    ## fit 2D gaussian to every point in parsedL2Output list
    logging.info("(__main__) found " + str(len(parsedOutput)) + " valid peaks") 
    logging.info("(__main__) fitting 2D gaussian to each entry in parsed peakfinder output list") 

    start = time.clock()

    allPeaks = []
    for line in parsedOutput:
        allPeaks.append([])

        thisX = int(round(float(line[0])))
        thisY = int(line[1])
        thisVal = float(line[2])

        data = []
        # construct sublist in data containing image data values to fit (essentially 2D list "array" containing data values)
        for i in range(thisX-((params['fittingApertureX']-1)/2), thisX+((params['fittingApertureX']+1)/2)):     # +1 as range is not inclusive of maximum value
            data.append([])
            for j in range(thisY-((params['fittingApertureY']-1)/2), thisY+((params['fittingApertureY']+1)/2)): # +1 as range is not inclusive of maximum value
                data[len(data)-1].append(arcFile.data[j, i])

        dataNumpy = np.asarray(data)    	# convert nested lists to numpy array
        outParams = fitgaussian(dataNumpy) 	# this is the bottleneck function call (specifically optimize.leastsq())

        allPeaks[len(allPeaks)-1].append(thisY)	# append row number for aiding plotting
        allPeaks[len(allPeaks)-1].append(outParams[0])
        allPeaks[len(allPeaks)-1].append(outParams[1] + thisX - ((params['fittingApertureX']-1)/2))
        allPeaks[len(allPeaks)-1].append(outParams[2] + thisY - ((params['fittingApertureY']-1)/2))
        allPeaks[len(allPeaks)-1].append(2.35*outParams[3])
        allPeaks[len(allPeaks)-1].append(2.35*outParams[4])

    elapsed = (time.clock() - start)
    logging.info("(__main__) gaussian fitting completed in " + str(elapsed) + "s")  

    ## REMOVE DUPLICATES
    ## as frclean finds centroids per row, duplicate entries in y will require removal
    ## n.b. duplicate peaks should have same parameters, so neither will be particularly favoured 
    logging.info("(__main__) removing duplicate entries") 

    start = time.clock()

    noDupPeaks = []
    thisIterFWHMX = []	# used for ascertaining mean x FWHM
    thisIterFWHMY = []	# used for ascertaining mean y FWHM
    # for each peak
    for el in range(len(allPeaks)):
        elRow = allPeaks[el][0]
        elHeight = allPeaks[el][1]
        elCentX = allPeaks[el][2]
        elCentY = allPeaks[el][3]      
        elFWHMX = allPeaks[el][4]
        elFWHMY = allPeaks[el][5]  
        # consider every other peak
        add = True
        for otherEl in range(el+1, len(allPeaks)):
            otherElCentX = allPeaks[otherEl][2]
            otherElCentY = allPeaks[otherEl][3]
            # and calculate cartesian pixel distance between centroids
            dist = sqrt(pow(elCentX-otherElCentX, 2) + pow(elCentY-otherElCentY, 2))
            if dist < params['duplicateMinDist']: 	# there's a duplicate further down the list
                add = False	
                break

        # add peak if no duplicate was found
        if add:
            noDupPeaks.append([])
            noDupPeaks[len(noDupPeaks)-1].append(elRow)
            noDupPeaks[len(noDupPeaks)-1].append(elHeight)
            noDupPeaks[len(noDupPeaks)-1].append(elCentX)
            noDupPeaks[len(noDupPeaks)-1].append(elCentY)
            noDupPeaks[len(noDupPeaks)-1].append(elFWHMX)
            noDupPeaks[len(noDupPeaks)-1].append(elFWHMY)
		
            thisIterFWHMX.append(elFWHMX)
            thisIterFWHMY.append(elFWHMY)

    elapsed = (time.clock() - start)
    logging.info("(__main__) duplicate entry removal completed in " + str(elapsed) + "s") 
    logging.info("(__main__) mean x FWHM is " + str(round(np.mean(thisIterFWHMX), 2)) + " +/- " + str(round(np.std(thisIterFWHMX), 2)) + "px")  
    logging.info("(__main__) mean y FWHM is " + str(round(np.mean(thisIterFWHMY), 2)) + " +/- " + str(round(np.std(thisIterFWHMY), 2)) + "px")  

    ## FILTER PEAK DATASET TO REMOVE REMAINING OUTLIERS
    # i) perform iterative cut of FWHM in noDupPeaks dataset
    # ii) remove single peak rows
    thisIterPeaks = list(noDupPeaks)
    numIterCycles = 0
    while True:
        numIterCycles += 1
        numBadPeaks = 0

        lastIterPeaks = list(thisIterPeaks)
        lastIterFWHMX = list(thisIterFWHMX)
        lastIterFWHMY = list(thisIterFWHMY)
        del thisIterPeaks[:]
        del thisIterFWHMX[:]
        del thisIterFWHMY[:]

        currentRow = lastIterPeaks[0][0]
        countEl = 0
        for el in lastIterPeaks:     
            elRow = el[0]
            elHeight = el[1]
            elCentX = el[2]
            elCentY = el[3]      
            elFWHMX = el[4]
            elFWHMY = el[5]

            if elFWHMX > np.mean(lastIterFWHMX) + params['iterSigmaCut']*np.std(lastIterFWHMX):		# filter FWHM X
                numBadPeaks += 1 
                continue
            elif elFWHMY > np.mean(lastIterFWHMY) + params['iterSigmaCut']*np.std(lastIterFWHMY):	# filter FWHM Y
                numBadPeaks += 1
                continue
            elif elRow != currentRow:									# filter for single element in row
                currentRow = elRow
                if countEl == 1:
                    countEl = 0
                    continue

            thisIterFWHMX.append(elFWHMX)
            thisIterFWHMY.append(elFWHMY)

            thisIterPeaks.append([])
            thisIterPeaks[len(thisIterPeaks)-1].append(elRow)
            thisIterPeaks[len(thisIterPeaks)-1].append(elHeight)
            thisIterPeaks[len(thisIterPeaks)-1].append(elCentX)
            thisIterPeaks[len(thisIterPeaks)-1].append(elCentY)
            thisIterPeaks[len(thisIterPeaks)-1].append(elFWHMX)
            thisIterPeaks[len(thisIterPeaks)-1].append(elFWHMY)

            countEl += 1

        if numBadPeaks == 0:
            logging.info("(__main__) filter succeeded after " + str(numIterCycles) + " cycles")
            break
        elif params['maxNumIterCycles'] == numIterCycles:
            im_err._setError(1)
            im_err.handleError()
            break

    logging.info("(__main__) mean x FWHM is now " + str(round(np.mean(thisIterFWHMX), 2)) + " +/- " + str(round(np.std(thisIterFWHMX), 2)) + "px")  
    logging.info("(__main__) mean y FWHM is now " + str(round(np.mean(thisIterFWHMY), 2)) + " +/- " + str(round(np.std(thisIterFWHMY), 2)) + "px")
    filteredPeaks = list(thisIterPeaks)

    ## PLOTTING
    logging.info("(__main__) generating plots")   
    plot = plots(im_err, filteredPeaks)

    # set up figure
    plt.figure(num=1, figsize=(7, 5), dpi=150)
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 6}
    matplotlib.rc('font', **font)
    plt.subplots_adjust(hspace=.3)

    # add desired plots
    plt.subplot(321)
    plot.scatterPlot('CENTX', 'CENTY', "Position along spatial axis (px)", "Position along dispersion axis (px)", 'bx')
    plt.subplot(322)
    plot.scatterPlot('FWHMX', 'FWHMY', "FWHM X (px)", "FWHM Y (px)", 'rx')
    plt.subplot(323)
    disp_m1, disp_c1, disp_m2, disp_c2, disp_mAv, disp_cAv = plot.dispPos_v_FWHM(params) # n.b. subscripts 1/2 denote FWHMX and FWHMY respectively
    plt.subplot(324)
    spat_m1, spat_c1, spat_m2, spat_c2, spat_mAv, spat_cAv = plot.spatPos_v_FWHM(params) # n.b. subscripts 1/2 denote FWHMX and FWHMY respectively
    plt.subplot(325)
    plot.FWHM_heatmap(params, 'x')
    plt.subplot(326)
    plot.FWHM_heatmap(params, 'y')

    if params['interactive']:
        plt.show()
    else:
        plt.savefig(params['outImage'])

    logging.info("Average dispersion FWHM best fit results. m=" + str(round(disp_mAv, 10)) + ", c=" + str(round(disp_cAv, 10)))
    logging.info("Average spatial FWHM best fit results. m=" + str(round(spat_mAv, 10)) + ", c=" + str(round(spat_cAv, 10)))

    # write header to file if doesn't exist
    if not os.path.exists(params['outFile']):
        with open(params['outFile'], "w") as f:    
            f.write("# arcpath\tdate\tutstart\tgratid\tdisp_m1\tdisp_c1\tdisp_m2\tdisp_c2\tdisp_mAv\tdisp_cAv\tspat_m1\tspat_c1\tspat_m2\tspat_c2\tspat_mAv\tspat_cAv\tmean_FWHMX\tstd_FWHMX\n")   

    # write output to file (m/c in X/Y/combined for dispersion/spatial axes & combined mean/std FWHM for all X/Y/axes
    with open(params['outFile'], "a") as f:
        f.write(os.path.basename(str(params['pathToArcFile'])) + '\t' + str(arcFile.headers['DATE']) + '\t' + str(arcFile.headers['UTSTART']) + '\t' + str(arcFile.headers['GRATID']) + '\t' + str(disp_m1) + '\t' + str(disp_c1) + '\t' + str(disp_m2) + '\t' + str(disp_c2) + '\t' + str(disp_mAv) + '\t' + str(disp_cAv) + '\t' + str(spat_m1) + '\t' + str(spat_c1) + '\t' + str(spat_m2) + '\t' + str(spat_c2) + '\t' + str(spat_mAv) + '\t' + str(spat_cAv) + '\t' + str(round(np.mean(thisIterFWHMX), 2)) + '\t' + str(round(np.std(thisIterFWHMX), 2)) + '\n')

    arcFile.closeFITSFile()

