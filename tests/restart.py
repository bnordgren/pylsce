'''
This file contains some code to support the analysis of restart files.
Restart files contain different variables with different dimension names.
'''

import netCDF4 as nc
import numpy.ma as ma
from aggregator import NetCDFCopier

def compare_summary(first, second, threshold=1.e-6) :
    retval = True
    differences = []
    for pft in range(13) : 
        difference = first[:,pft].sum() - second[:,pft].sum()
        differences.append(difference)
        retval = retval and (abs(difference) <= threshold)
    return (retval, differences)

def load_sechiba(year, variable='veget_max') : 
    sechiba = nc.Dataset('sechiba_restart_%d.nc' % year)
    return sechiba.variables[variable]

def load_stomate(year, variable='PFTpresent') : 
    stomate = nc.Dataset('stomate_restart_%d.nc' % year)
    return stomate.variables[variable]

def sumByPFT(data):
    totals = []
    for pft in range(13) : 
        totals.append(data[:,pft].sum())
    return totals

class RestartFileFixer (NetCDFCopier) : 
    """Tools to "fix" a restart file so that it can be recognized as a 
    geospatial dataset by GDAL. Also corrects the "missing_value" attribute
    such that it is of the same data type as the associated data, and 
    renames it to _FillValue so that command line tools will correctly 
    pick it up."""
    def __init__(self, basename, year) : 
        self._basename = basename
        self._year     = year
        self._openFixed()
        self._initProspects()

    def _openFixed(self) : 
        """Attempts to open the output (fixed) restart file"""
        self.fixedfile = '%s_fixed_%d.nc' % (self._basename, self._year)
        try : 
            self._ncfile = nc.Dataset(self.fixedfile, 'r+')
        except :
            self._ncfile = nc.Dataset(self.fixedfile, "w")
            self._initDomain()

    def _openExemplar(self) : 
        """Returns a handle to the input (not fixed) file."""
        return nc.Dataset('%s_%d.nc' % ( self._basename, self._year ) )

    def _initDomain(self) : 
        """Copies the y dimension into "lat" and the x dimension into "lon"."""
        source = self._openExemplar()
        lon = self._ncfile.createDimension('lon', len(source.dimensions['x']))
        lat = self._ncfile.createDimension('lat', len(source.dimensions['y']))


        # Create the geographic coordinate variables
        latvar = self._ncfile.createVariable('lat', 
                source.variables['nav_lat'].dtype, ('lat',) )
        lonvar = self._ncfile.createVariable('lon', 
                source.variables['nav_lon'].dtype, ('lon',) )

        # populate the variables
        latvar[:] = source.variables['nav_lat'][:,0]
        lonvar[:] = source.variables['nav_lon'][0,:]

        source.close()

    def copyVariable(self, varname) : 
        source = self._openExemplar()

        dims = list(source.variables[varname].dimensions[:])
    
        # substitute lat for y and lon for x
        for i in range(len(dims)) :
            if dims[i] == 'x' :
                dims[i] = 'lon'
            elif dims[i] == 'y' : 
                dims[i] = 'lat'
            

        # Ensure that all the required dimensions are present
        for dim in dims :
            self.copyDimension(dim) 

        # get the missing value
        missing = round(source.variables[varname].missing_value, -13)

        # Create the variable in the fixed file
        the_copy = self._ncfile.createVariable(varname, 
            source.variables[varname].dtype,
            dims, fill_value=missing)

        # copy the data to the fixed variable
        the_copy[:] = ma.masked_equal(source.variables[varname][:], missing)

        source.close()

        return the_copy


