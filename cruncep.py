"""Code to build a single year orchidee-parseable cruncep file from 
the individual-variable cruncep files."""

import trend as t
import netCDF4 as nc
import numpy as np 

#
# Rain and snowfall are missing from the variables below because they 
# require special processing in the partition_precip function. They
# come from the file named "rain", in the variable "Total_Precipitation"
# and they appear in the "Snowf" and "Rainf" variables in the output
# file.
#
file_vars = ['lwdown',
             'qair',
             'press',
             'swdown',
             'tair',
             'uwind',
             'vwind']
orchidee_varnames  = ['LWdown', 
             'Qair',
             'PSurf',
             'SWdown',
             'Tair',
             'Wind_E',
             'Wind_N']
orig_varnames = ['Incoming_Long_Wave_Radiation',
             'Air_Specific_Humidity',
             'Pression',
             'Incoming_Short_Wave_Radiation',
             'Temperature', 
             'U_wind_component',
             'V_wind_component']

dims = ('y', 'x')
c_dim = 'land'

def make_nav_latlon( lat, lon ) : 
    """Constructs the 2D variables nav_lat and nav_lon from the two
    1D variables lat and lon. Note that the ORCHIDEE forcing file 
    requires that we add half a cell to longitude and subtract half 
    a cell from latitude relative to what is present in the cruncep
    files."""
    nav_lat = np.empty( (lat.size, lon.size), dtype = lat.dtype)
    nav_lon = np.empty( (lat.size, lon.size), dtype = lon.dtype)

    lat_res = (np.max(lat[:]) - np.min(lat[:]))/(lat.size-1)
    lon_res = (np.max(lon[:]) - np.min(lon[:]))/(lon.size-1)

    for i in range(lat.size) : 
        nav_lat[i,:] = [ lat[i] - (lat_res/2.) ] * lon.size

    for i in range(lon.size) : 
        nav_lon[:,i] = [ lon[i] + (lon_res/2.) ] * lat.size

    return nav_lat, nav_lon

def init_time(d, y, timelen) : 
    """Creates the "tstep" dimension and the "timestp" and "time" 
    coordinate variables."""
 
    d.createDimension('tstep', None)

    tstep = d.createVariable('timestp', 'i', ('tstep',))
    tstep[:] = np.arange(1, timelen+1, dtype='i')
    tstep.units = 'timesteps since %04d-01-01 00:00:00' % y
    tstep.title = 'Time Steps'
    tstep.long_name = 'Time step axis'
    tstep.time_origin = ' %04d-JAN-01 00:00:00' % y
    tstep.tstep_sec = 21600.

    time = d.createVariable('time', 'f', ('tstep',))
    time[:] = (tstep[:]-1) * 21600.
    time.units = 'seconds since %04d-01-01 00:00:00' % y
    time.title = 'Time'
    time.long_name = 'Time axis'
    time.time_origin = ' %04d-JAN-01 00:00:00' % y
    time.calendar = 'noleap'
    
    
def partition_precip(ofile,ca,y): 
    """Partitions the total precipitation in the "rain" cruncep file
    into snowfall and rainfall. If the surface temperature is less than
    0 deg C, the precip is classified as snowfall, otherwise, rain.
    The output is compressed and stored in the provided output file."""

    fname = 'cruncep_rain_%04d.nc' % y
    d = nc.Dataset(fname)
    len_c_dim = len(ofile.dimensions[c_dim])
    tsteps = len(d.dimensions['time_counter'])

    tair   = ofile.variables['Tair']
    orig_v = d.variables['Total_Precipitation']

    # input precip is mm/6-hr. output is kg/m^2/s
    conversion = 1. / (3600.*6.)


    # create the output variables in the NetCDF file, and 
    # set the basic attributes
    snow = ofile.createVariable('Snowf', orig_v.dtype, ('tstep',c_dim),
            chunksizes=(1,len_c_dim), fill_value=1.e20)
    snow.title = "Snowfall"
    snow.units = "kg/m^2/s"
    snow.missing_value = 1.e20
    rain = ofile.createVariable('Rainf', orig_v.dtype, ('tstep',c_dim),
            chunksizes=(1,len_c_dim), fill_value=1.e20)
    rain.title = "Rainfall"
    rain.units = "kg/m^2/s"
    rain.missing_value = 1.e20

    # loop over the timesteps, partitioning and compressing
    # as we go.
    for i in range(tsteps) : 
        # initialize 0-filled buffers for this timestep
        s = np.zeros( (len_c_dim,) , dtype = orig_v.dtype)
        r = np.zeros( (len_c_dim,) , dtype = orig_v.dtype)
        
        # read in temps from file for all of timestep
        temps = np.squeeze(tair[i,:])

        # read in and compress precip for all of timestep
        precip = ca.compress(orig_v[i,:,:]) * conversion

        # partition and copy data
        for j in range(len_c_dim) : 
            if temps[j] > 273.15 : 
                r[j] = precip[j]
            else: 
                s[j] = precip[j]

        # write timestep to output file
        snow[i,:] = s
        rain[i,:] = r

def cruncep_year(y) : 
    """Creates an ORCHIDEE forcing file for the specified year by
    copying and/or manipulating the variables found in the individual
    cruncep files."""

    o_fname = 'cruncep_halfdeg_%04d.nc' % y
    ofile = None

    #for i in [ 4 ] : 
    for i in range(len(file_vars)) : 
        fname = 'cruncep_%s_%04d.nc' % (file_vars[i], y)
        print fname
        d = nc.Dataset(fname)

        # If first time thru loop, init the output file
        if ofile == None : 
            m = d.variables['mask'][:]
            ofile,ca = t.compressedAxesFactory(o_fname, dims, c_dim, m)
            len_c_dim = len(ofile.dimensions[c_dim])

            lat = d.variables['latitude'][:]
            lon = d.variables['longitude'][:]
            nav_lat = ofile.createVariable('nav_lat', lon.dtype, ('y','x'))
            nav_lon = ofile.createVariable('nav_lon', lat.dtype, ('y','x'))

            nav_lat[:], nav_lon[:] = make_nav_latlon( lat[:], lon[:] )
            
            tsteps = len(d.dimensions['time_counter'])
            init_time(ofile, y, tsteps)


        orig_v = d.variables[orig_varnames[i]]
        v = ofile.createVariable(orchidee_varnames[i], orig_v.dtype, 
                ('tstep',c_dim), chunksizes=(1,len_c_dim), fill_value=1.e20)
        v.title = orig_v.title
        v.units = orig_v.units
        v.missing_value = 1.e20
        
        # compress the grid for each timestep...
        for j in range(tsteps) : 
            if file_vars[i] == "swdown" :
                # fix known units issue with the shortwave radiation
                v[j,:] = ca.compress(np.squeeze(orig_v[j,:,:]/21600.))
            else : 
                v[j,:] = ca.compress(np.squeeze(orig_v[j,:,:]))

        d.close()

    partition_precip(ofile,ca,y)

    ofile.close()
