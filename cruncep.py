"""Code to build a single year orchidee-parseable cruncep file from 
the individual-variable cruncep files."""

import trend as t
import netCDF4 as nc
import numpy as np 

#
# We appear to be missing snowfall. Find out where snowfall is.
#
file_vars = ['lwdown',
             'qair',
             'press',
             'rain',
             'swdown',
             'tair',
             'uwind',
             'vwind']
orchidee_varnames  = ['LWdown', 
             'Qair',
             'Psurf',
             'Rainf',
             'SWdown',
             'Tair',
             'Wind_E',
             'Wind_N']
orig_varnames = ['Incoming_Long_Wave_Radiation',
             'Air_Specific_Humidity',
             'Pression',
             'Total_Precipitation',
             'Incoming_Short_Wave_Radiation',
             'Temperature', 
             'U_wind_component',
             'V_wind_component']

dims = ('y', 'x')
c_dim = 'land'

def make_nav_latlon( lat, lon ) : 
    """Constructs the 2D variables nav_lat and nav_lon from the two
    1D variables lat and lon."""
    nav_lat = np.empty( (lat.size, lon.size), dtype = lat.dtype)
    nav_lon = np.empty( (lat.size, lon.size), dtype = lon.dtype)

    for i in range(lat.size) : 
        nav_lat[i,:] = [ lat[i] ] * lon.size

    for i in range(lon.size) : 
        nav_lon[:,i] = [ lon[i] ] * lat.size

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
    
    

def cruncep_year(y) : 
    
    o_fname = 'cruncep_halfdeg_%04d.nc' % y
    ofile = None

    for i in range(len(file_vars)) : 
        fname = 'cruncep_%s_%04d.nc' % (file_vars[i], y)
        print fname
        d = nc.Dataset(fname)

        # If first time thru loop, init the output file
        if ofile == None : 
            m = d.variables['mask'][:]
            ofile,ca = t.compressedAxesFactory(o_fname, dims, c_dim, m)

            lat = d.variables['latitude'][:]
            lon = d.variables['longitude'][:]
            nav_lat = ofile.createVariable('nav_lat', lon.dtype, ('y','x'))
            nav_lon = ofile.createVariable('nav_lon', lat.dtype, ('y','x'))

            nav_lat[:], nav_lon[:] = make_nav_latlon( lat[:], lon[:] )
            
            tsteps = len(d.dimensions['time_counter'])
            init_time(ofile, y, tsteps)


        orig_v = d.variables[orig_varnames[i]]
        v = ofile.createVariable(orchidee_varnames[i], orig_v.dtype, ('tstep',c_dim))
        v.title = orig_v.title
        v.units = orig_v.units
        
        # compress the grid for each timestep...
        for j in range(tsteps) : 
            v[j,:] = ca.compress(np.squeeze(orig_v[j,:,:]))

        d.close()

    ofile.close()
