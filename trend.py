import statsmodels.api as sm
import numpy.fft as f
import numpy.ma as ma
import numpy as np
import netCDF4 as nc

class CompressedAxes (object) : 
    """Given a compressed dimension in a netcdf file, this class
    can calculate the compressed (1d) index given the individual (2d) indices,
    or it can calculate the individual indices given the compressed
    index. The strategy here implements the GDT netCDF conventions.
    
    This class does not interpret physical location on the ground.
    """
    def __init__(self, dataset, c_dim, format='F') : 
        """Initializes a CompressedAxes object given a NetCDF dataset
        and the name of the compressed dimension. """
        self._dataset = dataset
        self._c_dim   = c_dim
        self._initCompression()
        self._realIndices = None
        self._format = format
        self._masks = None

    def _initCompression(self) : 
        dims = self._dataset.variables[self._c_dim].compress
        self._dimnames = dims.split()
        
        dimshape = [] 
        for i in range(len(self._dimnames)) : 
            x = len(self._dataset.dimensions[self._dimnames[i]])
            dimshape.append(x)
        self._dimshape = dimshape

        dimfactors = [1] * len(self._dimnames)
        accum = 1
        for i in range (len(dimfactors)) : 
            dimfactors[-(i+1)] = accum
            accum = accum * dimshape[-(i+1)]
        self._dimfactors = dimfactors

    def getCompressedIndex(self, indices) : 
        """Calculates the compressed index given the indices."""
        c_idx = 0
        for i in range(len(indices)) : 
            c_idx = c_idx + (indices[i]*self._dimfactors[i])

        # convert to a "1-based"/Fortran array convention, if necessary
        if self._format == 'F' : 
            c_idx = c_idx + 1

        return c_idx

    def getIndices(self, c_idx) : 
        """Calculates the individual indices given the compressed index"""
        if self._format == 'F' : 
            c_idx = c_idx - 1

        indices = [0] * len(self._dimfactors)
        for i in range(len(self._dimfactors)) : 
            indices[i] = c_idx / self._dimfactors[i]
            c_idx = c_idx - (indices[i] * self._dimfactors[i])
        return tuple(indices)

    def uncompress(self, vector) : 
        """Given a compressed vector, produce an uncompressed 
        2d representation. The vector must be the same length
        as the compressed dimension in the NetCDF file."""
        gi = self.get_grid_indices()
        
        grid = ma.masked_all( self._dimshape, dtype=vector.dtype ) 
        grid[gi] = vector
        
        grid = self.mask_grid(grid)
            
        return grid

    def get_grid_indices(self) : 
        """returns the grid indices which correspond to the compressed vector representation"""
        if self._realIndices == None : 
            num_pts = len(self._dataset.dimensions[self._c_dim])
            c_dim_var = self._dataset.variables[self._c_dim][:]
            self._realIndices = [ [], [] ]
            for i in range(num_pts) : 
                ind = self.getIndices(c_dim_var[i])
                self._realIndices[0].append(ind[0])
                self._realIndices[1].append(ind[1])
        return self._realIndices
        

    def find_grid_indices(self, idx_2d)  :
        """looks up the vector index corresponding to the 2d grid"""
        gi = self.get_grid_indices() 
        stop = len(gi[0])
        i = 0 
        done = False
        while (i < stop) and not done :
            done = (idx_2d[0] == gi[0][i]) and (idx_2d[1] == gi[1][i])
            if not done: 
                i += 1
        if not done : 
                return None
        return i           
        
    def compress(self, grid) : 
        """Given a 2d grid, and the c_dim coordinate variable
        in the netCDF file, create and return a vector representing
        the compressed grid."""
        gi = self.get_grid_indices() 

        v=grid[gi]
        v = self.mask_vec(v)

        return v
        
    def is_masked(self) :
        """returns true if a mask has been stored in this compressor"""
        return self._masks is not None
        
    def mask_grid(self, grid) : 
        """returns the provided grid masked by the stored grid mask"""
        gmask = self.get_grid_mask() 
        if gmask is None : 
            return grid
        grid = ma.array(grid, mask=gmask)
        return grid
        
    def mask_vec(self, vec)  :
        """returns the provided vector masked by the stored vector mask"""
        vmask = self.get_vec_mask()
        if vmask is None : 
            return vec
        vec = ma.array(vec, mask=vmask)
        return vec
        
    def set_grid_mask(self, g_mask) : 
        """sets the vector and grid masks to the provided value
        
        Future compression and uncompression routines will be masked using
        the provided value.
        """
        v_mask = self.compress(g_mask) 
        self._masks = (g_mask, v_mask)
        
    def get_grid_mask(self) : 
        """returns the current grid mask, or None"""
        g_mask = None
        if self._masks is not None : 
            g_mask = self._masks[0]
        return g_mask
        
    def set_vec_mask(self, v_mask)  :
        """sets the vector and grid masks to the provided value
        
        Future compression and uncompression routines will be masked using
        the provided value.
        """
        g_mask = self.uncompress(v_mask)
        self._masks = (g_mask, v_mask)

    def get_vec_mask(self) : 
        """returns the current vector mask, or None"""
        v_mask = None
        if self._masks is not None : 
            v_mask = self._masks[1]
        return v_mask
    
    def remove_mask(self) :
        self._masks = None
        

def compressedAxesFactory(ncfile, dimnames, c_dim, mask=None, bmask=None, format='F') : 
    """Initializes a NetCDF file with the dimensions and coordinate
    variables necessary to support a compressed axis representation of 
    2D grids. Returns the dataset and a CompressedAxes instance."""
    d = nc.Dataset(ncfile, 'w')

    # make the compressed dimension and coordinate variable
    if bmask is None : 
        # cruncep input files have a "mask" variable, where "-1"
        # means "not a good point"
        bmask = np.array([mask[i,j]==(-1) for i in range(mask.shape[0]) for j in range(mask.shape[1]) ])
        bmask = bmask.reshape( (mask.shape[0], mask.shape[1]) )

    # make the uncompressed dimensions.
    d.createDimension(dimnames[0], bmask.shape[0])
    d.createDimension(dimnames[1], bmask.shape[1])

    num_good = np.count_nonzero(bmask==False)
    d.createDimension(c_dim, num_good)
    c_dim_var = d.createVariable(c_dim, 'i4', (c_dim,))
    c_dim_var.compress = '%s %s' % dimnames

    # populate the coordinate variable
    ca = CompressedAxes(d, c_dim,format)
    k = 0 
    for i in range(bmask.shape[0]) : 
        for j in range(bmask.shape[1]) : 
            if not bmask[i,j] : 
                c_dim_var[k] = ca.getCompressedIndex( (i,j))
                k = k + 1

    return d, ca
                
    
        
class ExtendUnlimited (object) : 
    """Opens a series of netcdf files, concatenating variables along
    the specified dimension. The other dimensions must be the same."""
    def __init__(self, fname_pat, sub, cat_dim) : 
        datasets = [ ] 
        for s in sub : 
            datasets.append(nc.Dataset(fname_pat % s))

        starts = [ 0 ]
        lengths = [ len(datasets[0].dimensions[cat_dim]) ]
        for i in range(1,len(datasets)) : 
            starts.append(starts[i-1] + lengths[i-1])
            lengths.append( len(datasets[i].dimensions[cat_dim]) )

        self._datasets = datasets
        self._starts   = starts
        self._lengths  = lengths
        self._cat_dim  = cat_dim

    def getExtIndex(self, ds, idx) : 
        return self._starts[ds] + idx

    def getDsIndex(self, e_idx) : 
        ds=0
        while ((ds < len(self._starts)) and 
              ((self._starts[ds] + self._lengths[ds]) < e_idx)) : 
            ds=ds+1
        return (ds, e_idx - self._starts[ds])

    def calcDsSlices(self, e_slice) : 
        start = e_slice.start
        stop  = e_slice.stop
        if start == None : 
            start = 0 
        if stop == None : 
            stop  = self._starts[-1] + self._lengths[-1]

        startpos = self.getDsIndex(start)
        stoppos  = self.getDsIndex(stop)

        ds_slice = slice(startpos[0], stoppos[0]+1)
        starting_slice = slice(startpos[1],None)
        ending_slice   = slice(stoppos[1])
        return (ds_slice, starting_slice, ending_slice)

    def getExtendedLength(self) : 
        return self._lengths[-1] + self._starts[-1]

    def getShape(self, varname, key, ext_dim) : 
        shape = []
        orig_shape = self._datasets[0].variables[varname].shape
        ext_len = self.getExtendedLength()
        for i in range(len(key)) : 
            if type(key[i]) == int : 
                shape.append(1)
            else : 
                stop = key[i].stop
                start = key[i].start
                if stop == None : 
                    if (i == ext_dim) : 
                        stop = ext_len
                    else:
                        stop = orig_shape[i]
                if start == None : 
                    start = 0
                shape.append(stop-start)

        return shape 

    def getDataKey(self, key, shape) : 
        data_key = []
        for i in range(len(key)) : 
            if type(key[i]) == int : 
                data_key.append(0)
            else : 
                data_key.append(shape[i])
        return data_key


    def getData(self, varname, key, ext_dim) : 
        ext_slice = key[ext_dim]
        if type(ext_slice) == int : 
            ds_coord = self.getDsIndex(ext_slice)
            newkey = list(key)
            newkey[ext_dim] = ds_coord[1]
            data = self._datasets[ds_coord[0]].variables[varname][key]
        else:
            ds_slices = self.calcDsSlices(ext_slice)
            if ds_slices[0].start+1 == ds_slices[0].stop :
                newkey = list(key)
                newkey[ext_dim] = slice(ds_slices[1].start,ds_slices[2].stop)
                data = self._datasets[ds_slices[0].start].variables[varname][newkey]
            else:
                shape = self.getShape(varname, key, ext_dim)
                data = ma.zeros(shape)
                    
                # read from the first ds
                ds_key = list(key)
                ds_key[ext_dim] = ds_slices[1]
                data_key = self.getDataKey(key,shape)
                ext_last = self._lengths[ds_slices[0].start]-ds_slices[1].start
                data_key[ext_dim] = slice(0, ext_last)
                data[data_key] = self._datasets[ds_slices[0].start].variables[varname][ds_key]

                # read from all the middle ds es
                ds_key[ext_dim] = slice(None)
                for i in range(ds_slices[0].start+1, ds_slices[0].stop-1) :
                    print i
                    tmp = ext_last + self._lengths[i]
                    data_key[ext_dim] = slice(ext_last, tmp)
                    data[data_key] = self._datasets[i].variables[varname][ds_key]
                    ext_last = tmp
                
                # read from the last ds
                ds_key[ext_dim] = ds_slices[2]
                tmp = ext_last + ds_slices[2].stop
                data_key[ext_dim] = slice(ext_last, tmp)
                data[data_key] = self._datasets[ds_slices[0].stop-1].variables[varname][ds_key]


        return data


    def getVar(self, varname) : 
        var = self._datasets[0].variables[varname]
        shape = []
        for i in range(len(var.dimensions)) : 
            if var.dimensions[i] == self._cat_dim : 
                shape.append(self._starts[-1] + self._lengths[-1])
                unlim_dim = i
            else : 
                shape.append(var.shape[i])

        ext_var = ExtendedVar(self, varname, unlim_dim)
        ext_var.shape = shape
        ext_var.dimensions = var.dimensions
        return ext_var
        
    def close(self) : 
        for d in self._datasets : 
            d.close()
        
class ExtendedVar(object) : 
    def __init__(self, unlim, varname, unlim_dim) : 
        self._unlim = unlim
        self._varname = varname
        self._unlim_dim = unlim_dim

    def __getitem__(self, key) : 
        return self._unlim.getData(self._varname, key, self._unlim_dim)


class ModelEquation (object) : 
    """Abstract superclass to orchestrate some common operations"""
    def residuals(self, y, x) : 
        """Given a list of independent values in x, calculate a 
        vector of residuals as 'y - y_fit'."""
        return y - self.fitted_values(x)
        
    def fitted_values(self, x) : 
        """Given a list of independent values in x, calculate the model's 
        predicted value at each point."""
        return np.squeeze(np.array([self.evaluate(x_i).real for x_i in x]))


class FourierSeries(object) : 
    """Class represents a fourier series having a specified fundamental
    period and a given number of harmonics. Note, fourier series have 
    plus and minus terms for each harmonic.  This implementation does not 
    include the bias term, so you need to explicitly include it before fitting
    the function to a data series. Hence, one harmonic (just the fundamental 
    frequency) results in two terms, two harmonics result in four, etc. """
    def __init__(self, period, harmonics) : 
        self._period = float(period)
        self._harmonics = harmonics
        self._coef = None

    def no_coefficients(self, x) : 
        """generates individual terms at the observation point x without
        multiplying by coefficients. This is suitable for producing terms
        used to fit the coefficients against data, or for multiplying 
        with the coefficient vector. This vector does not contain
        a constant term."""
        return np.array( [np.exp(1j*2*np.pi*(i/self._period)*x )
            for i in range(-self._harmonics, self._harmonics+1) if i!=0])

    def getNumTerms(self) : 
        return 2*self._harmonics

    def setCoefficients(self, coef) : 
        self._coef = coef

    def getCoefficients(self) : 
        return np.copy(self._coef)

    def evaluate(self, x) : 
        """Calculates the value of this truncated fourier series at the 
        specified observation point."""
        return (self.no_coefficients(x) * self._coef).sum()


class SingleTerm(object) : 
    """Abstract superclass to manage the special case where only a single
    term is managed"""
    def __init__(self, coef=None) : 
        self._coef = coef

    def getNumTerms(self) : 
        return 1

    def setCoefficients(self, coef) : 
        self._coef = coef[0]

    def getCoefficients(self) : 
        return np.array( [self._coef], dtype=complex)

class LinearTerm(SingleTerm) : 
    """Models a linear term (w/o a constant)"""
    def no_coefficients(self, x) : 
        return np.array( [x], dtype = complex)

    def evaluate(self, x): 
        return x * self._coef

class ConstantTerm (SingleTerm) : 
    """Models a constant term"""
    def no_coefficients(self, x) : 
        return np.array( [1], dtype=complex) 
    
    def evaluate(self, x) : 
        return self._coef
        

class AdditiveSeries (ModelEquation) : 
    """Tracks one or more component series which each supply terms to an
    overall series.""" 
    
    def __init__(self, series) : 
        self._series = series
        total_terms = 0 
        series_terms = [ ] 
        for s in series : 
            curterms = s.getNumTerms()
            series_terms.append(curterms)
            total_terms = total_terms + curterms

        self._series_terms = series_terms
        self._total_terms = total_terms
        

    def getNumTerms(self) : 
        return self._total_terms

    def no_coefficients(self, x) : 
        v = np.empty( (self._total_terms,), dtype=complex)
        cur_terms = 0 
        for i in range(len(self._series)) :
            t = self._series_terms[i]
            v[cur_terms:cur_terms+t] = self._series[i].no_coefficients(x)
            cur_terms= cur_terms + t

        return v
            

    def evaluate(self, x) : 
        """Sum up the contributions of each component series."""
        retval = 0 
        for s in self._series : 
            retval = retval + s.evaluate(x)
        return retval 

    def getDesignMatrix(self, time) : 
        """Given a list of times, generates and returns a design
        matrix suitable for use in a least squares regression."""
        numterms = self._total_terms
        numtimes = len(time)

        # construct the exogeneous design matrix
        matrix = np.empty( (numtimes, numterms), dtype=complex ) 
        for o in range(numtimes) : 
            matrix[o,:] = self.no_coefficients(time[o])

        # wrap data with matrix class for lin alg operations
        return np.matrix(matrix, copy=False)

    def fit(self, time, obs, X=None) : 
        """Fits the series. Returns the r squared value. Component 
        series have their coefficients set."""
        assert len(time) == len(obs), "Time and observation arrays must be the same length"
        numterms = self._total_terms
        numtimes = len(time)
        
        if X == None : 
            X = self.getDesignMatrix(time)
        z = np.matrix(np.reshape(obs, (len(obs),1)), copy=False)

        # fit data
        b_hat = (X.H * X).I * X.H * z

        self.setCoefficients(b_hat)

    def setCoefficients(self, coef) : 
        cur_terms = 0 
        for i in range(len(self._series)) :
            t = self._series_terms[i]
            self._series[i].setCoefficients(coef[cur_terms:cur_terms+t])
            cur_terms= cur_terms + t

    def getCoefficients(self) : 
        v = np.empty( (self._total_terms,), dtype=complex)
        cur_terms = 0 
        for i in range(len(self._series)) :
            t = self._series_terms[i]
            v[cur_terms:cur_terms+t] = self._series[i].getCoefficients()
            cur_terms= cur_terms + t

class AbstractTrendFinder (object) : 
    def __init__(self, obs, times, timeslice=slice(None,None,None)):
        self.obs = obs
        self.timeslice = timeslice
        self.times = times[timeslice]

    def fit_all_pixels(self) : 
        """Runs the pixel fit for each pixel in the dataset. Retains the 
        slope, stderr of the slope, and the pval for the slope parameter
        for each pixel. The returned ndarray has three columns: 
        slope, stderr of slope, and pvalue of slope. There is one row
        for each pixel."""

        numobs   = self.obs.shape[1]
        out      = np.empty( (numobs,3), dtype=np.float64) 

        for i in range(numobs) : 
            if (i % 1000 == 0) : 
                print i
            out[i,:] = self.pixelfit(i)

        return out
        
    def getObs(self,i) : 
        return self.obs[self.timeslice,i]



class TrendFinder (AbstractTrendFinder)  :
    """Finds the trend in time series data by first calculating the 
    residual to a fitted periodic function, then fitting a 
    line through the residual. Many of the result parameters are retained."""
    def __init__(self, obs, times, series, timeslice=slice(None,None,None)):
        super(TrendFinder, self).__init__(obs,times,timeslice)
        self.series = series
        self.X1 = series.getDesignMatrix(times)
        self.X2 = sm.add_constant(self.times)

    def pixelfit(self, i) : 
        """Fit the specified pixel"""
        # fit the purely periodic part
        obs_extract = self.getObs(i)
        self.series.fit(self.times, obs_extract, self.X1)

        # fit the linear trend to the residuals.
        residuals = self.series.residuals(obs_extract, self.times)
        results = sm.OLS(residuals, self.X2).fit()

        slope = results.params[1]
        slope_stderr = results.bse[1]
        slope_pval   = results.pvalues[1]
        return (slope, slope_stderr, slope_pval) 
        

class RTrendFinder (AbstractTrendFinder) : 
    """Finds the trend in the time series using one of the more robust
    statistical estimators implemented in R."""
    def __init__(self, obs, times, rclass, timeslice=slice(None,None,None)) :
        super(RTrendFinder, self).__init__(obs,times,timeslice)
        self.rclass = rclass

    def pixelfit(self, i): 
        """Fit the specified pixel using R"""
        obs_extract = self.getObs(i)
        r_fitter = self.rclass(obs_extract, self.times)
        results = r_fitter.fit()
        slope = results.params[1]
        slope_stderr = results.bse[1]
        slope_pval = getattr(results, 'pvalues', None)
        if len(slope_pval) >= 2 : 
            slope_pval = slope_pval[1]
        return (slope, slope_stderr, slope_pval)
    
        
def r_squared(y, f) : 

    assert len(y)==len(f), "Arrays y and f must be same size"
    
    ybar = y.mean()
    sstot = (( y - ybar )**2).sum()
    ssres = (( y - f )**2).sum()
    return 1 - (ssres/sstot)

        
        
        

def fft_grid(vector, outvar, d=0.25) : 
    """performs an fft individually for each pixel, along the time dimension."""
    for i in range(vector.shape[1]) : 
        if i%100 == 0 : 
            print i
        outvar[:,i] = f.fft(vector[:,i])
    

def pixelfit(obs, i, times, series, timeslice=slice(None,None,None)) : 
    """Fits a series object to the time series of data contained within
    a single pixel, optionally extracting the specified time period. 
    Coefficients of the fit have been stored in the 
    series. The rsquared of the fit is returned."""
    # extract the domain over which the fit applies
    t_extract    = times[timeslice]
    obs_extract  = obs[timeslice,i]
    return series.fit(t_extract, obs_extract)

    
def linefit(vector,i,X,timeslice=slice(None,None,None)) : 
    # construct a linefit for the specified pixel across all the years in the dataset. 
    # return the slope and the r^2.
    results = sm.OLS(vector[timeslice,i], X).fit()
    intercept, slope = results.params
    rsq = results.rsquared
    if   -(np.isfinite(rsq)) : 
        rsq = ma.masked
        slope = ma.masked
    return (slope, rsq)

def grid_linefit(grid, timevals=None, timeslice=slice(None,None,None)) :
    """A compressed spatiotemporal grid is provided. A line fit is performed 
    along the time axis for each spatial cell. Two grids are returned,
    each of which is 2d, with the same spatial shape as the input.
    The pixels of one grid contains the slope, the other contains the 
    r squared value of the line fit for that spatial cell.
    A vector of time values may be provided. If not supplied, one 
    will be generated."""
    if timevals == None : 
        timevals = ma.arange(grid.shape[0])
    X = sm.add_constant(timevals, prepend=True)

    outshape = (grid.shape[1],)

    rsq_map = ma.zeros(outshape)
    slope_map = ma.zeros(outshape)

    for i in range(outshape[0]) : 
        if (i%1000) == 0 : 
            print "%d of %d (%f)" % (i, outshape[0], (i*100.0)/outshape[0])
        if ((type(grid) == 'numpy.ma.core.MaskedArray') 
            and grid[0,:].mask[i]) : 
            rsq_map[i] = ma.masked
            slope_map[i] = ma.masked 
        else : 
            m, rsq = linefit(grid,i,X,timeslice)
            rsq_map[i] = rsq
            slope_map[i] = m 

    return (slope_map, rsq_map)


