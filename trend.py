import statsmodels.api as sm
import numpy.ma as ma
import numpy as np
import netCDF4 as nc

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
        startpos = self.getDsIndex(e_slice.start)
        stoppos  = self.getDsIndex(e_slice.stop)

        ds_slice = slice(startpos[0], stoppos[0]+1)
        starting_slice = slice(startpos[1],None)
        ending_slice   = slice(stoppos[1])
        return (ds_slice, starting_slice, ending_slice)

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
                shape = []
                orig_shape = self._datasets[0].variables[varname].shape
                for i in range(len(key)) : 
                    if type(key[i]) == int : 
                        shape.append(1)
                    else : 
                        stop = key[i].stop
                        start = key[i].start
                        if stop == None : 
                            stop = orig_shape[i]
                        if start == None : 
                            start = 0
                        shape.append(stop-start)
                data = np.zeros(shape)
                    
                # read from the first ds
                ds_key = list(key)
                ds_key[ext_dim] = ds_slices[1]
                data_key = list(key)
                ext_last = self._lengths[ds_slices[0].start]-ds_slices[1].start
                data_key[ext_dim] = slice(0, ext_last)
                data[data_key] = self._datasets[ds_slices[0].start].variables[varname][ds_key]

                # read from all the middle ds es
                ds_key[ext_dim] = slice(None)
                for i in range(ds_slices[0].start+1, ds_slices[0].stop-1) :
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


def linefit(vector,i,X) : 
    # construct a linefit for the specified pixel across all the years in the dataset. 
    # return the slope and the r^2.
    results = sm.OLS(vector[:,i], X).fit()
    intercept, slope = results.params
    rsq = results.rsquared
    if   -(np.isfinite(rsq)) : 
        rsq = ma.masked
        slope = ma.masked
    return (slope, rsq)

def grid_linefit(grid, timevals=None) :
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
        if not grid[0,:].mask[i] : 
            m, rsq = linefit(grid,i,X)
            rsq_map[i] = rsq
            slope_map[i] = m 
        else : 
            rsq_map[i] = ma.masked
            slope_map[i] = ma.masked 

    return (slope_map, rsq_map)


