"""Code to build a single year orchidee-parseable forcing file from 
the individual-variable wfdei files."""

import trend as t
import netCDF4 as nc
import numpy as np 

#
# We appear to be missing Wind Direction (we have Wind Speed).
# Current thought is that SPITFIRE does not take into account the 
# wind direction, so we just align the wind vector along one of 
# the axes.
#
# I do not at this time have confirmation that the non-SPITFIRE 
# portions of ORCHIDEE also ignore the wind vector direction.
#
file_vars = ['LWdown_WFDEI',
             'Qair_WFDEI',
             'PSurf_WFDEI',
             'Rainf_WFDEI_GPCC',
             'Snowf_WFDEI_GPCC',
             'SWdown_WFDEI',
             'Tair_WFDEI',
             'Wind_WFDEI']
#             'vwind']
orchidee_varnames  = ['LWdown', 
             'Qair',
             'Psurf',
             'Rainf',
             'Snowf',
             'SWdown',
             'Tair',
             'Wind_E']
#             'Wind_N']
orig_varnames = ['LWdown',
             'Qair',
             'PSurf',
             'Rainf',
             'Snowf',
             'SWdown',
             'Tair',
             'Wind']
#             'V_wind_component']

dims = ('y', 'x')
c_dim = 'land'

def make_nav_latlon( lat, lon ) : 
    """Constructs the 2D variables nav_lat and nav_lon from the two
    1D variables lat and lon. Note that the ORCHIDEE forcing file 
    requires that we add half a cell to longitude and subtract half 
    a cell from latitude relative to what is present in the wfdei
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
    tstep.tstep_sec = 10800.

    time = d.createVariable('time', 'f', ('tstep',))
    time[:] = (tstep[:]-1) * 10800.
    time.units = 'seconds since %04d-01-01 00:00:00' % y
    time.title = 'Time'
    time.long_name = 'Time axis'
    time.time_origin = ' %04d-JAN-01 00:00:00' % y
    time.calendar = 'noleap'
    
    
def null_wind_component(d, c_dim, len_c_dim) : 
    nowind = d.createVariable('Wind_N', np.float32, ('tstep', c_dim), 
                fill_value=1.e20 )
    nowind.title = "North component of wind vector (null)"
    nowind.units = "m/s"
    nowind.missing_value = 1.e20
    nullgrid = np.zeros( (len_c_dim,) )
    for t in range(len(d.dimensions['tstep'])) :
        nowind[t,:] = nullgrid

def wfdei_year(y) : 
    
    o_fname = 'wfdei_halfdeg_%04d.nc' % y
    ofile = None

    for i in range(len(file_vars)) : 
        cum_tsteps = 0
        for i_month in range(12) :
            fname = '%s/%s_%04d%02d.nc' % (file_vars[i], file_vars[i], y, i_month+1)
            print fname
            d = nc.Dataset(fname)
            orig_v = d.variables[orig_varnames[i]]

            # If first time thru loop, init the output file
            if ofile == None : 
                m = orig_v[0,:].mask
                ofile,ca = t.compressedAxesFactory(o_fname, dims, c_dim, m)
                len_c_dim = len(ofile.dimensions[c_dim])
    
                lat = d.variables['lat'][:]
                lon = d.variables['lon'][:]
                nav_lat = ofile.createVariable('nav_lat', lon.dtype, ('y','x'))
                nav_lon = ofile.createVariable('nav_lon', lat.dtype, ('y','x'))
    
                nav_lat[:], nav_lon[:] = make_nav_latlon( lat[:], lon[:] )
                
                #tsteps = len(d.dimensions['time_counter'])
                init_time(ofile, y, 365*8)


            if i_month == 0 : 
                v = ofile.createVariable(orchidee_varnames[i], orig_v.dtype, 
                    ('tstep',c_dim), chunksizes=(1,len_c_dim), fill_value=1e20)
                v.title = orig_v.title
                v.units = orig_v.units
                v.missing_value = 1.e20
        
            # compress the grid for each timestep...
            tsteps = len(d.dimensions['tstep'])
            for j in range(tsteps) : 
                v[j+cum_tsteps,:] = ca.compress(np.squeeze(orig_v[j,:,:]))

            cum_tsteps += tsteps

            d.close()

    null_wind_component(ofile, c_dim, len_c_dim)
    ofile.close()
