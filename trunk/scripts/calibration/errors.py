import logging

class errors:
    '''
    A class for handling errors.
    '''
    def __init__(self):
        self._errorCode		= 0
	self._errorCodeDict	= {0:"No errors encountered",
				   1:"(__main__) Iterative cut failed. Too many cycles", 
				   -1:"(openFITSFile) Unable to open FITS file. File not found", 
                                   -2:"(getHeaders) Unable to get headers",
                                   -3:"(getData) Unable to get data",
                                   -4:"(closeFITSFile) Unable to close file. File not open",
                                   -5:"(__main__) fittingApertureX must be odd",
                                   -6:"(__main__) fittingApertureY must be odd"
}

    def _setError(self, newErrorCode):
        '''
        Set internal error code.
        '''
        self._errorCode = newErrorCode
        return True

    def handleError(self):
        '''
        Handle internal error code.
        '''
        errorMsg = self._errorCodeDict.get(self._errorCode)
        if self._errorCode is 0:
            logging.info(errorMsg)
        elif self._errorCode < 0:
            logging.critical(errorMsg)
            exit(0)
        elif self._errorCode > 0:
            logging.warning(errorMsg)
