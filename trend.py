import statsmodels.api as sm
import numpy.ma as ma
import numpy as np

class CompressedAxes (object) : 
    """Given a compressed dimension in a netcdf file, this class
    can calculate the compressed index given the individual indices,
    or it can calculate the individual indices given the compressed
    index. The strategy here implements the GDT netCDF conventions."""
    def __init__(self, dataset, c_dim) : 
        """Initializes a CompressedAxes object given a NetCDF dataset
        and the name of the compressed dimension."""
        self._dataset = dataset
        self._c_dim   = c_dim
        self._initCompression()

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

        return c_idx

    def getIndices(self, c_idx) : 
        """Calculates the individual indices given the compressed index"""
        indices = [0] * len(self._dimfactors)
        for i in range(len(self._dimfactors)) : 
            indices[i] = c_idx / self._dimfactors[i]
            c_idx = c_idx - (indices[i] * self._dimfactors[i])
        return indices
            
        
        
def linefit(grid,i,j,X) : 
    # construct a linefit for the specified pixel across all the years in the dataset. 
    # return the slope and the r^2.
    results = sm.OLS(grid[:,i,j], X).fit()
    intercept, slope = results.params
    rsq = results.rsquared
    if   -(np.isfinite(rsq)) : 
        rsq = ma.masked
        slope = ma.masked
    return (slope, rsq)

def grid_linefit(grid, timevals=None) :
    """A spatiotemporal grid is provided. A line fit is performed 
    along the time axis for each spatial cell. Two grids are returned,
    each of which is 2d, with the same spatial shape as the input.
    The pixels of one grid contains the slope, the other contains the 
    r squared value of the line fit for that spatial cell.
    A vector of time values may be provided. If not supplied, one 
    will be generated."""
    if timevals == None : 
        timevals = ma.arange(grid.shape[0])
    X = sm.add_constant(timevals, prepend=True)

    outshape = (grid.shape[1], grid.shape[2])

    rsq_map = ma.zeros(outshape)
    slope_map = ma.zeros(outshape)

    for i in range(outshape[0]) : 
        for j in range(outshape[1])  :
            if not grid[0,:,:].mask[i,j] : 
                m, rsq = linefit(grid,i,j,X)
                rsq_map[i,j] = rsq
                slope_map[i,j] = m 
            else : 
                rsq_map[i,j] = ma.masked
                slope_map[i,j] = ma.masked 

    return (slope_map, rsq_map)


