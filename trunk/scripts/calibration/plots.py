import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import pylab
import logging
import scipy.interpolate

class plots:
    '''
    A class for plotting calibration data
    '''
    def __init__(self, error, peaks):
	self.error	= error
        self.peaks	= peaks
        self.axesEnum   = {'ROW' 	: 0,
                           'HEIGHT' 	: 1,
                           'CENTX'	: 2,
                           'CENTY'	: 3,
                           'FWHMX'	: 4,
                           'FWHMY'	: 5,
                          }

    def scatterPlot(self, axis1, axis2, xtitle, ytitle, linestyle, hold=False):
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)

        for i, el in enumerate(self.peaks):
            thisElRow		= el[0]
            thisElX		= el[self.axesEnum[axis1]]
            thisElY		= el[self.axesEnum[axis2]]   

            # plot current x, y
            plt.plot(thisElX, thisElY, linestyle)
           
    def dispPos_v_FWHM(self, params):
        plt.xlabel("Position along dispersion axis (px)")
        plt.ylabel("FWHM (px)")

        xAv = []
        y1Av = []
        y2Av = []

        x = []
        y1 = []
        y2 = []

        # sort elements into lists
        for el in self.peaks:
            x.append(el[self.axesEnum['CENTY']])
            y1.append(el[self.axesEnum['FWHMX']])  
            y2.append(el[self.axesEnum['FWHMY']]) 

        for i in range(0, params['dispCCDSize']-params['yStep'], params['yStep']):
            # set limits
            xLo = i
            xHi = i+params['xStep']
            
            # find values between limits
            thisIterX = [x[idx] for idx, el in enumerate(x) if xLo < el <= xHi]
            thisIterY1 = [y1[idx] for idx, el in enumerate(x) if xLo < el <= xHi]
            thisIterY2 = [y2[idx] for idx, el in enumerate(x) if xLo < el <= xHi]

            # check we have a non-empty result
            if thisIterX:
                xAv.append(np.mean(thisIterX))
                y1Av.append(np.mean(thisIterY1))
                y2Av.append(np.mean(thisIterY2))

                plt.errorbar(np.mean(thisIterX), np.mean(thisIterY1), xerr=np.std(thisIterX), yerr=np.std(thisIterY1), color='b', linewidth=2)
                plt.errorbar(np.mean(thisIterX), np.mean(thisIterY2), xerr=np.std(thisIterX), yerr=np.std(thisIterY2), color='r', linewidth=2)
                plt.errorbar([], [], xerr=[], yerr=[], color='y', linewidth=2)
 
        # add best fit lines
        m1, c1 = pylab.polyfit(xAv, y1Av, 1) 
        m2, c2 = pylab.polyfit(xAv, y2Av, 1) 
        mAv, cAv = pylab.polyfit(xAv + xAv, y1Av + y2Av, 1)
        xFit = np.arange(min(xAv), max(xAv))
        y1Fit = (m1*xFit) + c1
        y2Fit = (m2*xFit) + c2
        yAvFit = (mAv*xFit) + cAv

        plt.plot(xFit, y1Fit, 'b') 
        plt.plot(xFit, y2Fit, 'r') 
        plt.plot(xFit, yAvFit, 'y')

        plt.legend(["FWHM_X", "FWHM_Y", "FWHM_AV"], prop={'size':6})

        return m1, c1, m2, c2, mAv, cAv

    def spatPos_v_FWHM(self, params):
        plt.xlabel("Position along spatial axis (px)")
        plt.ylabel("FWHM (px)")

        xAv = []
        y1Av = []
        y2Av = []

        x = []
        y1 = []
        y2 = []

        # sort elements into lists
        for el in self.peaks:
            x.append(el[self.axesEnum['CENTX']])
            y1.append(el[self.axesEnum['FWHMX']])  
            y2.append(el[self.axesEnum['FWHMY']]) 

        for i in range(0, params['spatCCDSize']-params['xStep'], params['xStep']):
            # set limits
            xLo = i
            xHi = i+params['yStep']
            
            # find values between limits
            thisIterX = [x[idx] for idx, el in enumerate(x) if xLo < el <= xHi]
            thisIterY1 = [y1[idx] for idx, el in enumerate(x) if xLo < el <= xHi]
            thisIterY2 = [y2[idx] for idx, el in enumerate(x) if xLo < el <= xHi]

            # check we have a non-empty result
            if thisIterX:
                xAv.append(np.mean(thisIterX))
                y1Av.append(np.mean(thisIterY1))
                y2Av.append(np.mean(thisIterY2))

                plt.errorbar(np.mean(thisIterX), np.mean(thisIterY1), xerr=np.std(thisIterX), yerr=np.std(thisIterY1), color='b', linewidth=2)
                plt.errorbar(np.mean(thisIterX), np.mean(thisIterY2), xerr=np.std(thisIterX), yerr=np.std(thisIterY2), color='r', linewidth=2)
                plt.errorbar([], [], xerr=[], yerr=[], color='y', linewidth=2)
 
        # add best fit lines
        m1, c1 = pylab.polyfit(xAv, y1Av, 1)
        m2, c2 = pylab.polyfit(xAv, y2Av, 1) 
        mAv, cAv = pylab.polyfit(xAv + xAv, y1Av + y2Av, 1)
        xFit = np.arange(min(xAv), max(xAv))
        y1Fit = (m1*xFit) + c1
        y2Fit = (m2*xFit) + c2
        yAvFit = (mAv*xFit) + cAv

        plt.plot(xFit, y1Fit, 'b') 
        plt.plot(xFit, y2Fit, 'r') 
        plt.plot(xFit, yAvFit, 'y')

        plt.legend(["FWHM_X", "FWHM_Y", "FWHM_AV"], prop={'size':6})

        return m1, c1, m2, c2, mAv, cAv

    def FWHM_heatmap(self, params, whichAxisFWHM):
        plt.xlabel("Position along spatial axis (px)")
        plt.ylabel("Position along dispersion axis (px)")

        # arrays for x, y, z1 (FWHM in x), z2 (FWHM in y)
        x = []
        y = []
        z1 = []
        z2 = []

        # prefilter arrays
        X = []
        Y = []
        Z1 = []
        Z2 = []

        # sort elements into lists
        for el in self.peaks:
            X.append(el[self.axesEnum['CENTX']])
            Y.append(el[self.axesEnum['CENTY']])  
            Z1.append(el[self.axesEnum['FWHMX']]) 
            Z2.append(el[self.axesEnum['FWHMY']]) 

        for i in range(0, params['spatCCDSize']-params['xStep'], params['xStep']):
            # set limits
            xLo = i
            xHi = i+params['xStep']
            # filter x (spatial) values between limits
            filterX_X = [X[idx] for idx, el in enumerate(X) if xLo < el <= xHi]
            filterX_Y = [Y[idx] for idx, el in enumerate(X) if xLo < el <= xHi]
            filterX_Z1 = [Z1[idx] for idx, el in enumerate(X) if xLo < el <= xHi]
            filterX_Z2 = [Z2[idx] for idx, el in enumerate(X) if xLo < el <= xHi]
            for j in range(0, params['dispCCDSize']-params['yStep'], params['yStep']):
                # set limits
                yLo = j
                yHi = j+params['yStep']
                # filter y (dispersion) values between limits
                filterXY_X = [filterX_X[idx] for idx, el in enumerate(filterX_Y) if yLo < el <= yHi]
                filterXY_Y = [filterX_Y[idx] for idx, el in enumerate(filterX_Y) if yLo < el <= yHi]
                filterXY_Z1 = [filterX_Z1[idx] for idx, el in enumerate(filterX_Y) if yLo < el <= yHi]
                filterXY_Z2 = [filterX_Z2[idx] for idx, el in enumerate(filterX_Y) if yLo < el <= yHi]
                # check we have a non-empty result
                if filterXY_X:
                    x += filterXY_X
                    y += filterXY_Y
                    z1 += filterXY_Z1
                    z2 += filterXY_Z2

        x = np.asarray(x)
        y = np.asarray(y)
        z1 = np.asarray(z1)
        z2 = np.asarray(z2)

        # which axis FWHM are we plotting?
        if whichAxisFWHM == 'x':
            z = z1
        elif whichAxisFWHM == 'y':
            z = z2

        # set up a regular grid of interpolation points
        xi, yi = np.linspace(x.min(), x.max(), params['xStep']), np.linspace(y.min(), y.max(), params['yStep'])
        xi, yi = np.meshgrid(xi, yi)

        # interpolate
        rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
        zi = rbf(xi, yi)
                
        plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], aspect='auto')
        plt.scatter(x, y, c=z, s=1, marker='.')
        plt.colorbar()
                
