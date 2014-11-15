import netCDF4 as nc 
import numpy as np
import numpy.linalg as la
import numpy.ma as ma


class AffineTransform (object)  :
    def __init__(self, matrix) : 
        self._forward = matrix
        self._backward = la.inv(matrix)

    def forward(self, vector) : 
        return self._forward.dot( vector)

    def xform_fw(self, x, y) : 
        return self.forward(np.array([x, y, 1]))

    def backward(self, vector) : 
        return self._backward.dot( vector)

    def xform_bw(self, i, j) : 
        return self.backward(np.array([i,j,1]))
        

def make_transform(x0, y0, xsize, ysize, i0=0, j0=0) : 
    xt0 = x0 - xsize * i0
    yt0 = y0 - ysize * j0
    matrix = np.array( [[xsize, 0, xt0], [0, ysize, yt0], [0, 0, 1]])
    return AffineTransform(matrix)

class NetCDFCopier (object) :  
    """An abstract superclass to support the copying of variables and
    dimensions between netcdf files. The child class needs to set a 
    self._ncfile property for the output and needs to provide an 
    "_openExemplar()" method which returns a handle to an open 
    netCDF dataset from the source file."""

    def _initProspects(self) : 
        """Grabs a list of variables from one of the input files."""
        f = self._openExemplar()
        self.prospects = f.variables.keys()
        f.close()

    def copyDimension(self, dim) : 
        """If the named dimension does not exist in the output file, 
        it is created and the coordinate variable is copied."""
        if dim in self._ncfile.dimensions : 
            return

        infile = self._openExemplar()
        d = infile.dimensions[dim] 
        self._ncfile.createDimension(dim, len(d))

        if dim in infile.variables :
            dim_copy = self._ncfile.createVariable(dim, 
                infile.variables[dim].dtype,
                (dim,))

            dim_copy[:] = infile.variables[dim][:]

        infile.close()
        
    def createDimension(self, dim, size) : 
        self._ncfile.createDimension(dim,size)

    def close(self) : 
        self._ncfile.close()



class NetCDFTemplate ( NetCDFCopier ) :
    """A class to open a new NetCDF file and copy dimensions as needed
    from a template NetCDF file."""

    def __init__(self, templatefile, newfile) : 
        self._ncfile = nc.Dataset(newfile, 'w')
        self._templatefilename = templatefile
        self._filename = newfile

    def _openExemplar(self) : 
        return nc.Dataset(self._templatefilename)

    def create_variable(self, name, dims, dtype, fill=None, chunk=None) : 
        """Creates a variable in the netcdf file with the given name,
        dimensions and data type."""

        # copy any missing dimensions from the template file
        for dim in dims : 
            self.copyDimension(dim)

        tgt = self._ncfile
        return tgt.createVariable(name, dtype, dims, fill_value=fill, chunksizes=chunk)


    def add_variable(self, data, name, dims, dtype=None) : 
        """Creates a variable in the netcdf file with the given
        name, dimensions and data. If the dimensions do not yet exist 
        in the output file, they are copied from the template. If dtype
        is not specified, it is set to be the same as the data type of 
        the data. This method copies the data into the netCDF file 
        after creating the variable.
        If successful, the new variable is returned."""

        fill_value = getattr(data, 'fill_value', None)
        newvar = self.create_variable(name, dims, dtype, fill_value)
        newvar[:] = data[:]

        return newvar

        

class OrchideeAggregation(NetCDFCopier) : 
    """This class remembers some key characteristics about a source 
    fileset so that aggregate variables may be written to an output 
    file."""
    def __init__(self, basename, cpus, year, fname=None) : 
        self._basename = basename
        self._cpus     = cpus
        self._year     = year
        if fname == None : 
            fname = "%s_%d.nc" % ( basename, year ) 
        self._fname = fname
        self._openAggregate()
        self._initProspects()

    def _openAggregate(self) : 
        """Checks for an existing file. Opens it if possible, creates
        it if not."""
        try : 
            self._ncfile   = nc.Dataset(self._fname, 'a')
        except RuntimeError : 
            self._ncfile   = nc.Dataset(self._fname, 'w')
            self._initDomain()

    def _openExemplar(self) : 
        """Opens one of the input files and returns the handle. You must
        close it. """
        return self.openSource(0)

    def openSource(self, cpu) : 
        """Opens one of the input files and returns the handle. You must
        close it. You provide the cpu number associated with the file. """
        return nc.Dataset("%s_%04d_%d.nc" % (self._basename, cpu, self._year))

    def _initDomain(self) : 
        """Initializes the spatial domain information in the output file."""
        infile = self._openExemplar()
        domain = infile.DOMAIN_size_global
        self._ncfile.createDimension('lat', domain[1])
        self._ncfile.createDimension('lon', domain[0])
        
        dtype = infile.variables['lat'].dtype
        latvar = self._ncfile.createVariable('lat', dtype, ('lat',))
        lonvar = self._ncfile.createVariable('lon', dtype, ('lon',))

        inlat = infile.variables['lat']
        latsz = inlat[1]-  inlat[0]

        inlon = infile.variables['lon']
        lonsz = inlon[1] - inlon[0]

        # indices in the file are 1-based (origin is at (1,1))
        origin_local = infile.DOMAIN_position_first

        origin_lon = inlon[0] - lonsz*(origin_local[0]-1)
        origin_lat = inlat[0] - latsz*(origin_local[1]-1)

        latvar[:] = np.arange(start=origin_lat, 
                        stop=(origin_lat+domain[1]*latsz), 
                        step=latsz, dtype=dtype)
        lonvar[:] = np.arange(start=origin_lon, 
                        stop=(origin_lon+domain[0]*lonsz), 
                        step=lonsz, dtype=dtype)

        infile.close()

    def aggregate_var(self, var) : 
        # Don't aggregate same variable more than once
        if var in self._ncfile.variables : 
            return

        agg, dims, dtype = aggregate_var(self._basename, 
                self._cpus, self._year, var)

        # Ensure that all necessary dimensions exist
        for dim in dims : 
            self.copyDimension(dim)

        # Create the output variable
        aggvar = self._ncfile.createVariable(var, dtype, dims)

        # Store the data
        aggvar[:] = agg[:]

        return aggvar

    def get_var(self, var):
        """Aggregates the named variable (if necessary) and returns it."""
        
        if var not in self._ncfile.variables : 
            self.aggregate_var(var)

        return self._ncfile.variables[var]

    def list_domains(self) : 
        """Collects the local domains of each component file."""
        domains = [] 
        for cpu in range(self._cpus) : 
            src = self.openSource(cpu)
            first = src.DOMAIN_position_first
            last  = src.DOMAIN_position_last

            domains.append( (first, last) ) 
            src.close()

        return domains

    def locate_source(self, lat, lon) : 
        """Locates the source file which contains the specified point."""
        cpu = 0
        found = False
        while (not found) and (cpu < self._cpus) : 
            src = self.openSource(cpu)
            first = src.DOMAIN_position_first
            last  = src.DOMAIN_position_last
            # make sure lat/lon are in ascending order for comparison
            lat_range = [first[1], last[1]]
            lon_range = [first[0], last[0]]
            lat_range.sort()
            lon_range.sort()
            found = ( (lat_range[0] <= lat) and (lat <= lat_range[1]) and
                      (lon_range[0] <= lon) and (lon <= lon_range[1]))
            if not found : 
                cpu+=1
            src.close()

        if not found : 
            cpu = None

        return cpu

        

def aggregate_var(basename, cpus, year, variable) :
    first = True
    for cpu in range(cpus) : 
        fname = "%s_%04d_%d.nc" % (basename, cpu, year)
        nc_file = nc.Dataset(fname) 

        if first : 
            first = False
            dtype = nc_file.variables[variable].dtype
            local_shape = nc_file.variables[variable].shape
            dims  = nc_file.variables[variable].dimensions

            domain = nc_file.DOMAIN_size_global

            global_shape = []
            i_dim=0
            window = [] 
            all_data = slice(None)
            for dim in dims : 
                if dim == 'lat' : 
                    global_shape.append(domain[1])
                    i_lat = i_dim
                    window.append(None)
                elif dim == 'lon' : 
                    global_shape.append(domain[0])
                    i_lon = i_dim
                    window.append(None)
                else :
                    global_shape.append(local_shape[i_dim])
                    window.append(all_data)
                    
                i_dim = i_dim + 1
                
            
            raw_window = list(window)
            raw_window[i_lon] = all_data
            raw_window[i_lat] = slice(1,None)
            agg = ma.masked_equal(np.zeros(global_shape, dtype=dtype),0)

        l_origin = nc_file.DOMAIN_position_first - 1
        l_max    = nc_file.DOMAIN_position_last 

        # Orchidee files "overlap" in the y/lat direction. These 
        # overlapping lines need to be merged with data which has already
        # been loaded into the aggregate array. In order to not
        # make assumptions about how cpu numbering relates to layout
        # on the global domain, or even about the amount of overlap,
        # information in the new array is allowed to overwrite only 
        # those cells which are masked out in the aggregate array.
        raw      = nc_file.variables[variable][:]

        window[i_lat] = slice(l_origin[1], l_max[1])
        window[i_lon] = slice(l_origin[0], l_max[0])
        orig = agg[window]
        ind = np.where(orig.mask)
        orig[ind] = raw[ind]
        agg[window] = orig

        nc_file.close()

    return (agg, dims, dtype)

    
class OrchideeSubdataset (object) : 
    def __init__(self, dataset) : 
        
        self.origin = dataset.DOMAIN_local_start
        self.shape = dataset.DOMAIN_local_size

        # forward transforms go from global index to local index.
        self._xform = make_transform(0,0,1,1,
            i0=self.origin[0]-1, j0=self.origin[1]-1)

        self.dataset = dataset
    

class OrchideeSeries (object) : 
    """This eventually may represent a virtual aggregation of multiple
    files which uses a spatial index to fetch data from the correct 
    partition."""

    def __init__(self, basename, cpus=64, yr=2000) : 
        """Opens a series of netCDF files comprising the ORCHIDEE history files
        for a year. These are taken to be one file per thread. The history
        files are either from the sechiba module or the stomate module, 
        as determined by the basename."""

        for cpu in range(cpus) : 
            fname = "%s_%04d_%d.nc" % (basename, cpu, yr)
            nc_file = nc.Dataset(fname)
            sub_ds = OrchideeSubdataset(nc_file)
            

        self._idx = idx
        self.cpus = cpus
        self.year = yr
        self.basename = basename

        
    
    def bounds(self) : 
        return self._idx.bounds

    def get_variable(self, varname) : 
        bounds = self.bounds()
        
    
