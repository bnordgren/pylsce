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
        self._rfit = rfit
        self.params = np.array(rfit.rx('coefficients'))

    def summary(self) : 
        """Prints a summary of the statistical fit."""
        base = importr('base')
        print base.summary(self._rfit)
        

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


