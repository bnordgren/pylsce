"""Module wraps the rpy2 interface to leverage R's more comprehensive 
implementation of robust fits. This module is limited to a single 
predictor variable and a single response variable."""

import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects.packages import importr

# enable automatic conversion between numpy objects and R
ro.conversion.py2ri = numpy2ri


class RFitResult (object) : 
    """Wraps an R fit object and exposes the results in a manner similar
    to statsmodels RegressionResults"""
    def __init__(self, rfit) : 
        self._base = importr('base')
        self._rfit = rfit
        self.params = np.squeeze(np.array(rfit.rx('coefficients')))
        self.df_resid = rfit.rx('df.residual')[0]
        self.resid = np.array(rfit.rx('residuals'))
        self.nobs = len(self.params)
        self.fittedvalues = rfit.rx('fitted.values')
        self._summary = self._base.summary(self._rfit)
        c = np.squeeze(np.array(self._summary.rx('coefficients')[0]))
        self.bse = c[:,1]
        if c.shape[1] >= 4 : 
            self.pvalues = c[:,3]

    def tvalues(self):  
        c = np.squeeze(np.array(self._summary.rx('coefficients')[0]))
        retval = None
        if c.shape[1] >=3 : 
            retval = c[:,2]
        return retval

    def summary(self) : 
        """Prints a summary of the statistical fit."""
        print self._summary
        

class RFit1D (object) : 
    def __init__(self, y, x) : 
        self.y = y
        self.x = x
        d = { "x" : ro.FloatVector(self.x), 
              "y" : ro.FloatVector(self.y) } 
        self.df = ro.DataFrame(d)
        self.fmla = ro.Formula('y ~ x')
        

    def fit(self) : 
        """Performs the fit"""
        retval = RFitResult(self._doFit())
        return retval


class MMEstimator (RFit1D) : 
    def _doFit(self) : 
        robustbase = importr('robustbase')
        return robustbase.lmrob(self.fmla, self.df)

class TheilSen (RFit1D) : 
    def _doFit(self) : 
        mblm = importr('mblm')
        return mblm.mblm(self.fmla, self.df)


