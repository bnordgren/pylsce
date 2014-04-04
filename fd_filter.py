"""
Module to support frequency domain filtering of time domain signals. Includes
a filter object, some filter constructors, and some time domain test signal 
generators...
"""
import scipy.signal as s
import numpy as np
import numpy.ma as ma 
import math as m
import aggregator as a
import netCDF4 as nc
import statsmodels.api as sm

class Sampling (object) : 
    """An object to convert samples to period, calculate nyquist,
    calculate normalized angles, etc."""
    def __init__(self, sample_period) : 
        self._sample_period = sample_period
        self._sample_rate = 1./sample_period
        self._nyquist = self._sample_rate / 2.

    def nyquist(self) : 
        return self._nyquist

    def rate(self) : 
        return self._sample_rate

    def period(self) : 
        return self._sample_period

    def normalize(self, freq, rad=False) : 
        if rad : 
            freq = freq / (2 * m.pi)
        return freq / self._nyquist

    def denormalize(self, n_freq, rad=False) : 
        freq = n_freq * self._nyquist
        if rad : 
            freq = 2 * m.pi * freq
        return freq

class PiecemealAverage ( object ) : 
    """Allows you to compute an average of a vector in a piecemeal 
    fashion (adding one vector at a time). When ready, request the 
    averaged vector."""
    def __init__(self, firstvector=None) : 
        if firstvector == None : 
            self._aggregate = None
            self._num = 0
        else : 
            self._aggregate = firstvector.copy()
            self._num = 1
        self._avg = None

    def includeVector(self, vector) : 
        """Adds another vector to the vector average."""
        if self._aggregate == None : 
            self._aggregate = vector.copy()
        else : 
            self._aggregate = self._aggregate + vector
        self._num = self._num + 1
        self._avg = None

    def average(self) : 
        if self._avg == None : 
            self._avg = self._aggregate / self._num
        return self._avg
        
class FrequencyFilter (object) : 

    def __init__(self, b, a, samp) : 
        """Coefficients (y value) and frequencies (x value) of the 
        desired response function. The coefficients should be real
        valued. Frequencies should be sorted in ascending order."""
        self._b = b
        self._a = a
        self._zi = s.lfilter_zi(b,a)
        self._samp = samp
        

    def filter(self, signal) : 
        """Filters the time-domain signal "signal" with the filter represented
        by this object. The sampling period of the signal must be regular
        and must be the same as the period used to design the filter
        parameters. Signal must be an array-like."""
        filteredsig, zo = s.lfilter(self._b, self._a, signal, 
            zi=signal[0]*self._zi)
        return filteredsig

    def filtfilt(self, signal) : 
        """Uses the scipy.signal.filtfilt() function to filter the signal. 
        The sampling period of the signal must be regular and must be the 
        same as the period used to design the filter parameters. Signal
        must be an array-like."""
        return s.filtfilt(self._b, self._a, signal)
        

def get_timefilter() : 
    samp = Sampling(0.25) # sampling period is 0.25 days
    wp = samp.normalize(1./(365*4))     # cutoff freq is time period of 4 years
    ws = samp.normalize(1./(365*2.5))   # time periods of 2.5 years should be rejected
    b, a  = s.iirdesign(wp, ws, 0.5, 6, ftype='butter')
    return FrequencyFilter(b,a,samp)


def process_cruncep(cruncep_file, proc_file, varname='Tair') : 
    # open files
    pfile = a.NetCDFTemplate(cruncep_file, proc_file)
    cfile = nc.Dataset(cruncep_file)

    # initialize
    tstep_len = len(cfile.dimensions['tstep'])
    tstep_dot_len = tstep_len -1
    pfile._ncfile.createDimension('tstep_dot', tstep_dot_len)
    filt = get_timefilter()
    v = cfile.variables[varname]
    fill = v.missing_value
    lf_temp = pfile.create_variable('LF_temps', ('tstep','land'), 
                        'float32', fill, (tstep_len,1))
    lf_temp.units = 'K'
    lf_temp.long_name = 'low frequency temperatures'
    lf_temp.valid_min = v.valid_min
    lf_temp.valid_max = v.valid_max

    lf_temp_dot = pfile.create_variable('LF_temps_dot', ('tstep_dot','land'),
                        'float32', fill, (tstep_dot_len,1) )
    lf_temp_dot.units = 'K/decade'
    lf_temp_dot.long_name = 'time derivative of low frequency temperatures'

    lf_temp_dot_trend = pfile.create_variable('LF_temps_dot_trend',
                        ('land',), 'float32', fill)
    lf_temp_dot_trend.units = 'K/decade'
    lf_temp_dot_trend.long_name = 'average of LF_temps_dot over time period'

    # need two averaging objects
    gt_avg = PiecemealAverage()
    ltdg_avg = PiecemealAverage()
    
    
    # loop over all the land points
    temps = ma.zeros( (v.shape[0],), dtype='float64')
    to_decade = (365.*10.)/0.25
    for i in range(v.shape[1]) : 
        if ( (i%100) == 0) : print i 
        temps[:] = ma.masked_outside(v[:,i],v.valid_min,v.valid_max)
        ftemps = filt.filter(temps)
        lf_temp[:,i] = ftemps

        # include in the "global_temps" average
        gt_avg.includeVector(lf_temp[:,i])

        # calculate the differential
        dtemps = np.diff(ftemps) * to_decade
        lf_temp_dot[:,i] = dtemps

        ltdg_avg.includeVector(dtemps)

        lf_temp_dot_trend[i] = ma.mean(dtemps)

    # request the average vectors and store
    global_temps = gt_avg.average()
    gt = pfile.add_variable(global_temps, 'global_temps', ('tstep',),'float32')
    gt.units = 'K'
    gt.long_name = 'global average of LF_temps'

    lf_temp_dot_global = ltdg_avg.average()
    ltdg = pfile.add_variable(lf_temp_dot_global, 'LF_temps_dot_global',
            ('tstep_dot',), 'float32')
    ltdg.units = 'K/decade'
    ltdg.long_name = 'global average of LF_temps_dot'
    
    # compute the scalar trends by the various methods
    trend = ma.mean(lf_temp_dot_trend[:])
    pfile.LF_dot_trend_global = trend

    trend = ma.mean(lf_temp_dot_global)
    pfile.LF_temps_dot_global_trend=  trend

    # still need to do the line fit to global_temps

    
        
        


