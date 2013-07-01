from errors import *
import pyfits
import os

class FITSFile:
    '''
    A class for handling FITS files.
    '''
    def __init__(self, filepath, error):
        self.filePath		= filepath
	self.error		= error
	self.hduList		= None
	self.headers		= None
	self.data               = None
	self.fileOpen		= False

    def openFITSFile(self):
        '''
        Check a FITS file exists and try to open it.
        '''
	if os.path.exists(self.filePath):
            self.hduList = pyfits.open(self.filePath)
	    self.fileOpen = True
            return True
        else:
            self.error._setError(-1)
            self.error.handleError()
            return False

    def getHeaders(self, HDU):
        '''
        Retrieve headers from a specified HDU.
        '''
        if self.hduList is not None:
            self.headers = self.hduList[HDU].header
            return True
        else:
            self.error._setError(-2)
            self.error.handleError()
            return False

    def getData(self, HDU):
        '''
        Get data from a specified HDU
        '''
        if self.hduList is not None:
            self.data = self.hduList[HDU].data
            return True
        else:
            self.error._setError(-3)
            self.error.handleError()
            return False       

    def closeFITSFile(self):
        '''
        Close FITS file cleanly.
        '''
        if self.hduList is not None:
            self.hduList.close()
	    self.fileOpen = False
            return True
        else:
            self.error._setError(-4)
            self.error.handleError()
            return False
