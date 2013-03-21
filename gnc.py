"""
The name means general nc operators. 
"""
import matplotlib.pyplot as plt
import matplotlib as mat
import g
import numpy as np
import netCDF4 as nc
import pdb as pdb
import mathex
import copy
import pb
import bmap
from inspect import isfunction
import Pdata
import datetime
import pandas as pa
import os
import sys

home_dir = os.path.expanduser('~')
pylab_dir = home_dir+'/'+'python/python_lib'
basedata_dir = pylab_dir + '/base_data'

def append_doc_of(fun):
    def decorator(f):
        f.__doc__ += fun.__doc__
        return f
    return decorator

def find_index_by_region(globHDlon, globHDlat, region_lon, region_lat):
    """
    Find the four indices for the subregion against the whole region that will be formed by merging.

    Parameters:
    -----------
    globHDlon and globHDlat are vectors of lon/lat for the connected file;
    region_lon and region_lat are vectors of lon/lat for the individuel files.

    Warnings:
    ---------
    As this function is designed in case to merge some regional files to a global one (or connected one), in this case the globHDlon/globHDlat could be the lon/lat
        of the whole region that's to be formed by merging, and region_lon/region_lat could be the lon/lat of each subregion that participate in merging. So a strict
        rule has been applied that requires the globHDlon and region_lon (as well as lat) to be in the SAME ascending/descending order.

    Returns:
    --------
    lon_index_min, lon_index_max, lat_index_min, lat_index_max: the four indices indicating the position of the subregion in the bigger grid. Now to have region_lon
        by using lon_index_min and lon_index_max, one should use globHDlon[lon_index_min:lon_index_max+1] due to python's indexing method.
    """
    if not pb.indexable_check_same_order(globHDlon,region_lon):
        raise ValueError("the longtitudes are not in same order!")
    if not pb.indexable_check_same_order(globHDlat,region_lat):
        raise ValueError("the latitudes are not in same order!")

    lon_index = np.nonzero((globHDlon >= region_lon.min()) & (globHDlon <= region_lon.max()))[0]
    lon_index_min = np.min(lon_index)
    lon_index_max = np.max(lon_index)

    lat_index = np.nonzero((globHDlat >= region_lat.min()) & (globHDlat <= region_lat.max()))[0]
    lat_index_min = np.min(lat_index)
    lat_index_max = np.max(lat_index)
    return (lon_index_min, lon_index_max, lat_index_min, lat_index_max)

def find_index_by_vertex(globHDlon, globHDlat, (vlon1,vlon2), (vlat1,vlat2)):
    """
    globHDlon and globHDlat are vectors of lon/lat for the connected file;
    (vlon1,vlon2), (vlat1,vlat2) are (min, max) of lon/lat for the individuel files.
    Return:
        (lon_index_min, lon_index_max, lat_index_min, lat_index_max)
    Test:

    """
    lon_index = np.nonzero((globHDlon >= vlon1) & (globHDlon <= vlon2))[0]
    lon_index_min = np.min(lon_index)
    lon_index_max = np.max(lon_index)

    lat_index = np.nonzero((globHDlat >= vlat1) & (globHDlat <= vlat2))[0]
    lat_index_min = np.min(lat_index)
    lat_index_max = np.max(lat_index)
    return (lon_index_min, lon_index_max, lat_index_min, lat_index_max)

def test_find_index_by_vertex():
    d = pb.pfload(basedata_dir+'/landmask_et_latlon.pf')
    (lon_index_min, lon_index_max, lat_index_min, lat_index_max) = find_index_by_vertex(d.globHDlon, d.globHDlat, (-50,20),(-35,25))
    lon = d.globHDlon[lon_index_min: lon_index_max+1]
    assert lon[0]==-49.75 and lon[-1]==19.75
    lat = d.globHDlat[lat_index_min: lat_index_max+1]
    assert lat[0]==24.75 and lat[-1]==-34.75

def test_find_index_by_region():
    d = pb.pfload('base_data/landmask_et_latlon.pf')
    region_lon = np.arange(-49.75,20,0.5)
    region_lat = np.arange(24.75,-35,-0.5)
    (lon_index_min, lon_index_max, lat_index_min, lat_index_max) = find_index_by_region(d.globHDlon, d.globHDlat, region_lon, region_lat)
    reglon_new = d.globHDlon[lon_index_min: lon_index_max+1]
    np.testing.assert_array_equal(region_lon, reglon_new)
    reglat_new = d.globHDlat[lat_index_min: lat_index_max+1]
    np.testing.assert_array_equal(region_lat, reglat_new)

def _construct_slice_by_dim(numdim,(lon_index_min, lon_index_max, lat_index_min, lat_index_max)):
    """
    construct the slicing which we used to slice the global data in order to fill in regional data.
    """
    if numdim == 4:
        subslice=np.s_[:,:,lat_index_min:lat_index_max+1, lon_index_min:lon_index_max+1]
    elif numdim == 3:
        subslice=np.s_[:,lat_index_min:lat_index_max+1, lon_index_min:lon_index_max+1]
    elif numdim == 2:
        subslice=np.s_[lat_index_min:lat_index_max+1, lon_index_min:lon_index_max+1]
    else:
        raise ValueError("Strange that numdim is 1")
    return subslice


def dic2ncHalfDegree(filename,timevar=('year',1),vardic='none',histtxt='none'):
    """
    Definition:
        dic2nc(filename,timevar=('year',1),vardic='none',histtxt='none')
    Note:
        1. Only use for 0.5dX0.5d resolution
        2. the default dimensions are lat, lon, time_counter
        3. timevar=('year',1), the first element of the tuple is the longname of time variable, the second element is the time length.
        4. vardic=[('ba','burned area','ha',ba),('fc','fire counts','1/day',fc)], the first two in the variable tuple are shortname and longname for varialbe, 
           the third one is unit, and the last one is data. Note the data should be a numpy masked array with dimension of "timelength * 360 * 720".
        5. The function can handle the missing value automatically from the fill_value attribute of the masked array.
    Example:
        >>> a=np.arange(360*720.).reshape(360,720)
        >>> b=np.ma.masked_less(a,100000)
        >>> ff=np.ma.array([b,b])
        >>> ff.shape
          (2, 360, 720)
        >>> dic2nc('test.nc',timevar=('year',2),vardic=[('tv','test variable','no unit',ff)],histtxt='only for test purpose')
    """
    print 'Begin to write to netcdf file',filename
    rootgrp=nc.Dataset(filename,'w',format='NETCDF3_CLASSIC')
    
    #creat dimension
    lat=rootgrp.createDimension('lat',360)
    lon=rootgrp.createDimension('lon',720)
    rootgrp.createDimension('time_counter',None)
    
    #creat dimension variable
    lat=rootgrp.createVariable('lat','f4',('lat',))
    lon=rootgrp.createVariable('lon','f4',('lon',))
    lat.long_name = "latitude"
    lon.long_name = "longitude"
    lon.units="degrees_east"
    lat.units="degrees_north"
    y=np.arange(-179.75,180,0.5)
    x=np.arange(89.75,-90,-0.5)
    lat[:]=x
    lon[:]=y
    
    time=rootgrp.createVariable('time','i4',('time_counter',))
    time[:]=np.arange(1,timevar[1]+1)
    time.long_name = timevar[0]
    
    
    #creat ordinary variable
    if vardic=='none':
        print 'no data is provided!'
    else:
        for vardata in vardic: 
            ba=rootgrp.createVariable(vardata[0],'f4',('time_counter','lat','lon',))
            ba.long_name=vardata[1]
            ba.units=vardata[2]
            if vardata[3].ndim==2:
                ba[0,:,:]=vardata[3]
            else:
                ba[:]=vardata[3]
            if np.ma.isMA(vardata[3]):
                ba.missing_value=vardata[3].fill_value
            else:
                ba.missing_value=float(-9999.)

            print 'Variable  --',vardata[0],'--  is fed into the file'
    rootgrp.history=histtxt
    rootgrp.close()


def dic_ndarray_to_nc_HalfDegree(filename,dic_of_ndarray,unit='unitless',histtxt='No other history provided'):
    """
    write a dictionary of ndarrays of HalfDegree to a netcdf file.

    Parameters:
    ----------
    dic_of_ndarray: a dictionary of ndarray, the dim of ndarray could only be 2-dim or 3-dim.

    Notes:
    ------
    1. The keys will be used as variable names in the nc file.

    """
    if mathex.ndarray_list_check_unique_shape(dic_of_ndarray.values()):
        array_shape = np.shape(dic_of_ndarray.values()[0])
    else:
        raise ValueError("Please check: not all the ndarray in the dictionary share the same dimension")

    if len(array_shape) <2 or len(array_shape) >3:
        raise ValueError("This function currently only accepts the ndarray with 2 or 3 dimensions, the given dimension is {0}".format(array_shape))
    else:
        if len(array_shape) == 2:
            timevar = ('year',1)
        else:
            timevar = ('year',array_shape[0])
        vardic=[]
        for key,array in dic_of_ndarray.items():
            vardic.append((key,key,unit,array))
    histtxt = 'file created at ' + str(datetime.datetime.today()) + '\n' + histtxt
    dic2ncHalfDegree(filename,timevar=timevar,vardic=vardic,histtxt=histtxt)

def txt2nc_HalfDegree(filename,name_list=None,varname_list=None,name_keyword=False,common_prefix_in_name='',name_surfix='.txt',unit='unitless',histtxt='No other history provided',land_mask=True):
    """
    Write one or a group of txt files directly to a netcdf file.

    Parameters:
    -----------
    name_list: the name_list allows a flexible way to write txt to netcdf files.
        1. when name_list is a string, write only from one txt file, varname_list should also be a string specifying the varialbe name that's in the new nc file.
        2. when name_list is a string:
            2.1 when name_keyword is True, the members in the name_list will be used as nameid together with common_prefix_in_name and name_surfix to construct the complete
                file names, the nameid will be automacatilly used as variable names. In this case varname_list will be overwritten if it's not None
            2.2 when name_keyword is False, the members in name_list should indicate full path of txt files, and strings in varname_list will be used as the varialble names
                for nc files.
    name_keyword: see above.
    land_mask: boolean variable. if True the input arrays will be applied with land mask.

    Notes
    -----
    1. the txt files are read by np.genfromtxt, so literally all file types that could be treated with np.genfromtxt could be used as input files.

    See also
    --------
    dic_ndarray_to_nc_HalfDegree
    dic2ncHalfDegree
    """
    land = pb.pfload(basedata_dir+'/landmask_et_latlon.pf')

    #prepare name_list when name_keyword is True
    if name_keyword and isinstance(name_list,list):
        varname_list = name_list[:]
        name_list = [common_prefix_in_name + nameid + name_surfix for nameid in name_list]

    dic_of_ndarray={}
    #we want to write from a list of files
    if isinstance(name_list,list):
        for nameid,varname in zip(name_list,varname_list):
            array = np.genfromtxt(nameid)
            if land_mask:
                array = np.ma.masked_array(array, mask=land.globlandmaskHD)
            dic_of_ndarray[varname] = array
    #we want to write from only one file
    elif isinstance(name_list,str):
        array = np.genfromtxt(name_list)
        if land_mask:
            array = np.ma.masked_array(array, mask=land.globlandmaskHD)
        dic_of_ndarray[varname_list] = array
    else:
        raise ValueError("name_list could only be string or list!")

    dic_ndarray_to_nc_HalfDegree(filename,dic_of_ndarray,unit=unit,histtxt=histtxt)

def dic2nc(filename,latvar=1,lonvar=1,timevar=('year',1),vardic='none',histtxt='none'):
    """
    Definition:
       dic2nc(filename,latvar=1,lonvar=1,timevar=('year',1),vardic='none',histtxt='none') 
    Note:
        1. latvar, lonvar receive lat and lon variables, eg. latvar=np.arange(89.75,-90,-0.5), lonvar=np.arange(-179.75,180,0.5)
        2. the default dimensions are (time_counter, lat, lon)
        3. timevar=('year',1), the first element of the tuple is the longname of time variable, the second element is the time length.
        4. vardic=[('ba','burned area','ha',ba),('fc','fire counts','1/day',fc)], the first two in the variable tuple are shortname and longname for varialbe, 
           the third one is unit, and the last one is data. Note the data should be a numpy masked array with dimension of "timelength * 360 * 720".
        5. The function can handle the missing value automatically from the fill_value attriabute of the masked array.
    Example:
        >>> a=np.arange(360*720.).reshape(360,720)
        >>> b=np.ma.masked_less(a,100000)
        >>> ff=np.ma.array([b,b])
        >>> ff.shape
          (2, 360, 720)
        >>> dic2nc('test.nc',latvar=np.arange(89.75,-90,-0.5),lonvar=np.arange(-179.75,180,0.5),timevar=('year',2),vardic=[('tv','test variable','no unit',ff)],histtxt='only for test purpose')
    """
    print 'Begin to write to netcdf file',filename
    rootgrp=nc.Dataset(filename,'w',format='NETCDF3_CLASSIC')
    
    #creat dimension
    lat=rootgrp.createDimension('lat',len(latvar))
    lon=rootgrp.createDimension('lon',len(lonvar))
    rootgrp.createDimension('time_counter',None)
    
    #creat dimension variable
    lat=rootgrp.createVariable('lat','f4',('lat',))
    lon=rootgrp.createVariable('lon','f4',('lon',))
    lat.long_name = "latitude"
    lon.long_name = "longitude"
    lon.units="degrees_east"
    lat.units="degrees_north"
    lat[:]=latvar
    lon[:]=lonvar
    
    time=rootgrp.createVariable('time','i4',('time_counter',))
    time[:]=np.arange(1,timevar[1]+1)
    time.long_name = timevar[0]
    
    
    #creat ordinary variable
    if vardic=='none':
        print 'no data is provided!'
    else:
        for vardata in vardic: 
            ba=rootgrp.createVariable(vardata[0],'f4',('time_counter','lat','lon',))
            ba.long_name=vardata[1]
            ba.units=vardata[2]
            if vardata[3].ndim==2:
                ba[0,:,:]=vardata[3]
            else:
                ba[:]=vardata[3]
            if np.ma.isMA(vardata[3]):
                ba.missing_value=vardata[3].fill_value
            else:
                ba.missing_value=1.e+20

            print 'Variable  --',vardata[0],'--  is fed into the file'
    rootgrp.history=histtxt
    rootgrp.close()


class NcWrite(object):
    """
    NcWrite object allows for flexible writing to NetCDF file.
    """

    def __init__(self,filename):
        print 'Begin to write to netcdf file',filename
        self.filename=filename
        self.rootgrp=nc.Dataset(filename,'w',format='NETCDF3_CLASSIC')
        self._diminfo_keys = ['dim_name', 'dimvar_name', 'dimvar_longname', 'dimvar_dtype', 'dimvar_value', 'dimvar_unit', 'unlimited']
        self.default_record_dim_name = 'time_counter'
        self.default_record_dim_var_name = 'time'
        self.dimensions = {}
        self.diminfo_lat = dict(zip(['dim_name', 'dimvar_name', 'dimvar_longname'],['lat', 'lat', 'latitude']))
        self.latdim_name = self.diminfo_lat['dim_name']
        self.latvar_name = self.diminfo_lat['dimvar_name']
        self.diminfo_lon = dict(zip(['dim_name', 'dimvar_name', 'dimvar_longname'],['lon', 'lon', 'longitude']))
        self.londim_name = self.diminfo_lon['dim_name']
        self.lonvar_name = self.diminfo_lon['dimvar_name']
        self.diminfo_time_year = dict(zip(['dim_name', 'dimvar_name', 'dimvar_longname', 'dimvar_dtype', 'dimvar_unit'],['time_counter','time', 'year', 'i4', 'year']))
        self.diminfo_time_month = dict(zip(['dim_name', 'dimvar_name', 'dimvar_longname', 'dimvar_dtype', 'dimvar_unit'],['time_counter','time', 'month', 'i4', 'month'])                                      )
        self.diminfo_time_day = dict(zip(['dim_name', 'dimvar_name', 'dimvar_longname', 'dimvar_dtype', 'dimvar_unit'],['time_counter','time', 'day', 'i4', 'day'])                                      )
        self.timedim_name = self.diminfo_time_year['dim_name']
        self.timevar_name = self.diminfo_time_year['dimvar_name']
        self.rootgrp.history=''

    @staticmethod
    def _replace_none_by_given(orinput,default):
        if orinput == None:
            return default
        else:
            return orinput

    def add_diminfo_lat(self,dim_name=None,dimvar_name=None,dimvar_longname=None):
        self.diminfo_lat = dict(zip(['dim_name', 'dimvar_name', 'dimvar_longname'],[dim_name,dimvar_name,dimvar_longname]))
        self.latdim_name = self.diminfo_lat['dim_name']
        self.latvar_name = self.diminfo_lat['dimvar_name']

    def add_diminfo_lon(self,dim_name=None,dimvar_name=None,dimvar_longname=None):
        self.diminfo_lon = dict(zip(['dim_name', 'dimvar_name', 'dimvar_longname'],[dim_name,dimvar_name,dimvar_longname]))
        self.londim_name = self.diminfo_lon['dim_name']
        self.lonvar_name = self.diminfo_lon['dimvar_name']

    def add_dim(self, diminfo_value = None,**attr_kwargs):
        """
        Write dimensions and dimensional variables to NetCDF file.
        diminfo could be supplied with a list: [dim_name, dimvar_name, dimvar_longname, dimvar_dtype, dimvar_value, dimvar_unit, unlimited], with unlimited as True or False.

        Returns:
        -------
        the dimension that's created.
        """

        assert len(diminfo_value) == len(self._diminfo_keys)
        diminfo = dict(zip(self._diminfo_keys, diminfo_value))

        #creat dimension and dimension variable
        if diminfo['unlimited'] == True:
            self.rootgrp.createDimension(diminfo['dim_name'], None)
        else:
            self.rootgrp.createDimension(diminfo['dim_name'], len(diminfo['dimvar_value']))

        ncdimvar = self.rootgrp.createVariable(diminfo['dimvar_name'], diminfo['dimvar_dtype'], (diminfo['dim_name'],))
        ncdimvar[:] = diminfo['dimvar_value']
        ncdimvar.long_name = diminfo['dimvar_longname']
        ncdimvar.units = diminfo['dimvar_unit']
        #set variable attributes by hand
        for key,value in attr_kwargs.items():
            ncdimvar.setncattr(key, value)
        self.dimensions[diminfo['dim_name']] = self.rootgrp.dimensions[diminfo['dim_name']]


    def add_dim_lat(self, latvar = None, **attr_kwargs):
        """
        A shortcut function to write latitude dimension, cf. add_dim for
            more information. default latvar is global half degree.
        """
        latvar = self._replace_none_by_given(latvar,np.arange(89.75,-90,-0.5))
        lat_diminfo_value = [self.diminfo_lat['dim_name'],
                             self.diminfo_lat['dimvar_name'],
                             self.diminfo_lat['dimvar_longname'],
                             'f4', latvar, 'degrees_north', False]
        self.add_dim(lat_diminfo_value, **attr_kwargs)
        self.latvar = latvar


    def add_dim_lon(self, lonvar = None, **attr_kwargs):
        """
        A shortcut function to write longitude dimension, cf. add_dim
            for more information. default lonvar is global half degree.
        """
        lonvar = self._replace_none_by_given(lonvar,np.arange(-179.75,180,0.5))
        lon_diminfo_value = [self.diminfo_lon['dim_name'],
                             self.diminfo_lon['dimvar_name'],
                             self.diminfo_lon['dimvar_longname'],
                             'f4', lonvar, 'degrees_east', False]
        self.add_dim(lon_diminfo_value, **attr_kwargs)
        self.lonvar = lonvar

    def add_dim_pft(self,**attr_kwargs):
        """
        A shortcut function to write PFT dimension, cf. add_dim
            for more information.
        """
        pft_diminfo_value = ['PFT', 'PFT', 'Plant functional type',
                              'i4', np.arange(1,14), '1', False]
        self.add_dim(pft_diminfo_value, **attr_kwargs)

    def add_dim_time(self, timevar = None, timestep='year', **attr_kwargs):
        """
        A shortcut fuction to write time dimension
            [default dimension name: time_counter;
             default dimension variable name: time],
            cf. add_dim for more information.
        """
        timevar = self._replace_none_by_given(timevar,np.arange(1,2))
        if timestep == 'year':
            diminfo_time = self.diminfo_time_year
        elif timestep == 'month':
            diminfo_time = self.diminfo_time_month
        elif timestep == 'day':
            diminfo_time = self.diminfo_time_day
        else:
            raise ValueError("timestep {0} is not recognized!".format(timestep))

        time_diminfo_value = [diminfo_time['dim_name'],
                              diminfo_time['dimvar_name'],
                              diminfo_time['dimvar_longname'],
                              diminfo_time['dimvar_dtype'],
                              timevar, diminfo_time['dimvar_unit'], True]
        self.add_dim(time_diminfo_value,**attr_kwargs)
        self.timedim_name = diminfo_time['dim_name']
        self.timevar_name = diminfo_time['dimvar_name']
        self.timevar = timevar


    #creat ordinary variable
    def add_var(self, varinfo_value = None, attr_copy_from = None, **attr_kwargs):
        """
        Purpose: Add variables to NetCDF file.
        Arguments:
            1. varinfo_value = [varname, dim_tuple, dtype, varvalue], eg.['ba', ('time_counter', 'lat', 'lon', ), 'f4', ba_data]
            2. use attr_copy_from if the variable attributes are copied from another netCDF4.Variable object.
            3. set variable attributes by using attr_kwargs.
        Notes:
        ------
        1. The function can handle the missing value automatically from the fill_value attriabute of the masked array.
        2. It can handle the case of time axis lenght ==1
        """
        varinfo_keys=['varname', 'dim_tuple', 'dtype', 'varvalue']
        varinfo = dict(zip(varinfo_keys, varinfo_value))
        var = self.rootgrp.createVariable(varinfo['varname'], varinfo['dtype'], varinfo['dim_tuple'])
        vardata = varinfo['varvalue']

        #handle the case when len(time) == 1
        if var.ndim - vardata.ndim == 1 and var.shape[0] == 1:
            var[:] = vardata[np.newaxis,...]
        else:
            var[:] = vardata

        if np.ma.isMA(vardata):
            var.missing_value = vardata.fill_value
        else:
            var.missing_value = 1.e+20
        #var.setncattr('_FillValue',var.missing_value)

        #copy the variable attributes from another netCDF4.Variable object.
        if isinstance(attr_copy_from, nc.Variable):
            for attr_name in attr_copy_from.ncattrs():
                var.setncattr(attr_name, attr_copy_from.getncattr(attr_name))

        #set variable attributes by hand
        for key,value in attr_kwargs.items():
            var.setncattr(key, value)
        print 'Variable  --',varinfo['varname'],'--  is fed into the file'

    def add_2dim_lat_lon(self, lat_value=None, lon_value=None):
        """
        This is shortcut function for adding 3dim_time_lat_lon; all 3 dimensions will use default dim name.
        """
        self.add_dim_lat(latvar=lat_value)
        self.add_dim_lon(lonvar=lon_value)

    def add_3dim_time_lat_lon(self, time_length=None, lat_value=None, lon_value=None):
        """
        This is shortcut function for adding 3dim_time_lat_lon; all 3 dimensions will use default dim name.
        """
        self.add_dim_lat(latvar=lat_value)
        self.add_dim_lon(lonvar=lon_value)
        time_length = self._replace_none_by_given(time_length,1)
        self.add_dim_time(np.arange(time_length)+1)

    def add_4dim_time_pft_lat_lon(self, time_length=None, lat_value=None, lon_value=None):
        """
        This is shortcut function for adding 3dim_time_lat_lon; all 3 dimensions will use default dim name.
        """
        self.add_dim_lat(latvar=lat_value)
        self.add_dim_lon(lonvar=lon_value)
        time_length = self._replace_none_by_given(time_length,1)
        self.add_dim_time(np.arange(time_length)+1)
        self.add_dim_pft()

    def add_var_3dim_time_lat_lon(self, varname, data, attr_copy_from=None,
                                  **attr_kwargs):
        self.add_var((varname,
                     (self.timedim_name, self.latdim_name,self.londim_name, ),
                     'f4', data), attr_copy_from=attr_copy_from, **attr_kwargs)

    def add_var_2dim_lat_lon(self, varname, data, attr_copy_from=None, **attr_kwargs):
        self.add_var((varname, (self.latdim_name, self.londim_name, ), 'f4', data), attr_copy_from=attr_copy_from, **attr_kwargs)

    def add_var_4dim_time_pft_lat_lon(self, varname, data, attr_copy_from=None, **attr_kwargs):
        self.add_var((varname, (self.timedim_name, 'PFT', self.latdim_name, self.londim_name, ), 'f4', data), attr_copy_from=attr_copy_from, **attr_kwargs)

    def add_history_attr(self,histtxt='No history provided'):
        self.rootgrp.history=histtxt

    def add_global_attributes(self,attr_dic):
        self.rootgrp.setncatts(attr_dic)

    def close(self):
        self.rootgrp.history = self.rootgrp.history + \
                    '\nfile created at ' + str(datetime.datetime.today())
        self.rootgrp.close()

    def _construct_data_by_dim(self,numdim):
        """
        construct the global ndarray that we need and fill in with the data from each region.
        """
        dimlen_dic = pb.Dic_Apply_Func(len,self.dimensions)
        if numdim == 4:
            glob_data = np.ma.zeros((dimlen_dic[self.timedim_name],dimlen_dic['PFT'],dimlen_dic[self.latdim_name],dimlen_dic[self.londim_name]))
            glob_data.mask = True
        elif numdim == 3:
            glob_data = np.ma.zeros((dimlen_dic[self.timedim_name],dimlen_dic[self.latdim_name],dimlen_dic[self.londim_name]))
            glob_data.mask = True
        elif numdim == 2:
            glob_data = np.ma.zeros((dimlen_dic[self.latdim_name],dimlen_dic[self.londim_name]))
            glob_data.mask = True
        else:
            raise ValueError("Strange that numdim is 1")
        return glob_data


    def add_var_smart_ndim(self,varname,numdim,data,pftdim=False,**attr_kwargs):
        '''
        Select the add_var_* method in a smart way by knowing the ndim.

        Parameters:
        -----------
        ndim: the number of dimensions for the varname
        '''
        if numdim == 4:
            if data.ndim == 4:
                self.add_var_4dim_time_pft_lat_lon(varname, data, **attr_kwargs)
            elif data.ndim == 3:
                if pftdim == False:
                    self.add_var_3dim_time_lat_lon(varname, data, **attr_kwargs)
                else:
                    self.add_var_4dim_time_pft_lat_lon(varname, data, **attr_kwargs)
            elif data.ndim == 2:
                self.add_var_3dim_time_lat_lon(varname, data, **attr_kwargs)
            else:
                raise ValueError("data ndim < 2!")
        elif numdim == 3:
            self.add_var_3dim_time_lat_lon(varname, data, **attr_kwargs)
        elif numdim == 2:
            self.add_var_2dim_lat_lon(varname, data, **attr_kwargs)
        else:
            raise ValueError("Strange that numdim is 1")



    def _add_var_by_dim(self,varname,numdim,data):
        """
        Used only by add_var_from_file_list
        """
        glob_data = self._construct_data_by_dim(numdim)
        for subdata in data:
            lon_index_min, lon_index_max, lat_index_min, lat_index_max = find_index_by_region(self.lonvar, self.latvar, subdata.lonvar, subdata.latvar)
            subslice = _construct_slice_by_dim(numdim,(lon_index_min, lon_index_max, lat_index_min, lat_index_max))
            glob_data[subslice] = subdata.d1.__dict__[varname] #note here the mask of glob_data will be changed automatically.
            print "data fed from file --{0}--".format(subdata.filename)

        if numdim == 4:
            self.add_var_4dim_time_pft_lat_lon(varname, glob_data, attr_copy_from=subdata.d0.__dict__[varname])
        elif numdim == 3:
            self.add_var_3dim_time_lat_lon(varname, glob_data, attr_copy_from=subdata.d0.__dict__[varname])
        elif numdim == 2:
            self.add_var_2dim_lat_lon(varname, glob_data, attr_copy_from=subdata.d0.__dict__[varname])
        else:
            raise ValueError("Strange that numdim is 1")

    def add_var_from_file_list(self,input_file_list,varlist,Ncdata_latlon_dim_name=None):
        """
        Mainly used for merging nc files spatially

        Parameters:
        -----------
        input_file_list: the NetCDF input file list for merging.
        varlist: variable name list that will appear in merged nc file. Note each input file have have the specified variable and the dimension accros all 
            input files must be the same.
        Ncdata_latlon_dim_name: the selective varialbe that's used in Ncdata when open a nc file.

        """
        data = [Ncdata(filename,latlon_dim_name=Ncdata_latlon_dim_name) for filename in input_file_list]
        subdata_first = data[0]
        for varname in varlist:
            if not pb.object_list_check_unique_attribute([subdata.d0.__dict__[varname] for subdata in data],'dimensions'):
                raise ValueError("The variable {0} in all input files does not have the same dimension!".format(varname))
            else:
                numdim = len(subdata_first.d0.__dict__[varname].dimensions)
                self._add_var_by_dim(varname,numdim,data)

        glob_attr_dic = subdata_first.global_attributes
        #pdb.set_trace()
        #add time stamp for this operation
        if 'history' in glob_attr_dic:
            glob_attr_dic['history'] = 'file created at ' + str(datetime.datetime.today()) + '\n' + \
                                       'by merging files:' + '\n'+ \
                                       '--'+('--'+'\n'+'--').join(input_file_list)+'--'+'\n' + \
                                        str(glob_attr_dic['history'])
        else:
            glob_attr_dic['history'] = 'file created at ' + str(datetime.datetime.today()) + '\n' + \
                                       'by merging files:' + '\n'+ \
                                       '--'+('--'+'\n'+'--').join(input_file_list)+'--'

        #write the global attributes
        if nc.__version__ == '1.0.1':
            self.add_global_attributes(glob_attr_dic)
        elif nc.__version__ == '0.9.7':
            for key,value in glob_attr_dic.items():
                self.rootgrp.setncattr(key,value)

def nc_spatial_concat_ncfiles(outfile,input_file_list,timestep=None,time_length=None,
                   varlist=None,latvar=None,lonvar=None,
                   latinfo=None,loninfo=None,pft=False,
                   Ncdata_latlon_dim_name=None):
    """
    A shortcut for merging spatially a list of files. For detailed control
        of dimensions and varialbes, use NcWrite first, followed by
        add_var_from_file_list.

    Parameters:
    -----------
    timestep: 'year' or 'month', default is 'year'
    time_length: np.arange(1,time_length+1) will be the time variable value.
    varlist: the variable list that's to be retained in merged file.
        default inclules all variables except the dimension variable.
    latvar,lonvar: the lat/lon for megered data. default is 0.5 degree
        resolution with global coverage.
    latinfo,loninfo: tuple containing ('lat/lon_dim_name','lat/lon_var_name',
        'lat/lon_var_longname'); default for lat is ('lat','lat','latitude')
        and for lon is ('lon','lon','longitude').
    pft: if pft==True, then PFT dimensions from ORCHIDEE will be added.
    Ncdata_latlon_dim_name: the specified lat/lon dimension names that are
        used when calling Ncdata.

    see also
    --------
    gnc.Ncdata

    Test
    ----
    nc_subgrid_csv and nc_merge_files are tested against each other
        in the gnc_test.py.
    """
    ncfile = NcWrite(outfile)
    ncfile.add_diminfo_lon('lon','lon','longitude')
    ncfile.add_diminfo_lat('lat','lat','latitude')
    ncfile.add_2dim_lat_lon(np.arange(89.75,-90,-0.5),
                            np.arange(-179.75,180,0.5))
    if pft == True:
        ncfile.add_dim_pft()
    if timestep == None:
        timestep = 'year'
    ncfile.add_dim_time(np.arange(1,1+time_length),timestep=timestep)
    subdata_first = Ncdata(input_file_list[0],
                           latlon_dim_name=Ncdata_latlon_dim_name)
    if varlist == None:
        varlist = pb.StringListAnotB(subdata_first.list_var(),
                                     subdata_first.dimvar_name_list)
    ncfile.add_var_from_file_list(input_file_list,varlist,
                         Ncdata_latlon_dim_name=Ncdata_latlon_dim_name)
    ncfile.close()

def _set_default_ncfile_for_write(ncfile,**kwargs):
    '''
    set default lat,lon,time,pft... ect for writing to nc file.

    kwargs:
    -------
    timestep: 'year' or 'month', default is 'year'
    time_length: np.arange(1,time_length+1) will be the time variable value,
        default value is 1.
    latvar,lonvar: the lat/lon for megered data. default is 0.5 degree
        resolution with global coverage.
    latinfo,loninfo: tuple containing ('lat/lon_dim_name','lat/lon_var_name',
        'lat/lon_var_longname'); default for lat is ('lat','lat','latitude')
        and for lon is ('lon','lon','longitude').
    pft: if pft==True, then PFT dimensions from ORCHIDEE will be added.
    '''

    loninfo = kwargs.get('loninfo',('lon','lon','longitude'))
    latinfo = kwargs.get('latinfo',('lat','lat','latitude'))
    latvar = kwargs.get('latvar',np.arange(89.75,-90,-0.5))
    lonvar = kwargs.get('lonvar',np.arange(-179.75,180,0.5))
    timestep = kwargs.get('timestep','year')
    time_length = kwargs.get('time_length',1)

    ncfile.add_diminfo_lon(*loninfo)
    ncfile.add_diminfo_lat(*latinfo)
    ncfile.add_2dim_lat_lon(latvar,lonvar)
    ncfile.add_dim_time(np.arange(1,1+time_length),timestep=timestep)
    if kwargs.get('pft',False):
        ncfile.add_dim_pft()



def ncm2y(infile,outfile,mode=None,varlist='all'):
    """
    Purpose: Transform monthly data to yearly mean or sum. Note the input file format is very strict.
    Note:
        1. the input netCDF file must have only 3 (2or4 will fail) dimensions and only 3 dimension varialbes with the corresponding axis as the only dimension.
        2. the input netCDF file other variables must have only 3 dimesions with the first dimension as the unlimited time dimension.
        3. mode = 'sum' or 'mean', use mathex.m2ymean or mathex.m2ysum to do the transformation.
        4. the function copy the dimension name, variable name and data for the other 2 dimensions other than the time dimension direclty from the input file to output file.
        5. by default, if varlist is 'all', all the variables in input file will be calculated.

        6.This function has been tested against NCO ncra and is working correctly.
    Example:
    >>> import gnc
    >>> gnc.ncm2y('testdata/cru1999prm.nc','testdata/cru1999prm.yearsum.nc',mode='sum')
    """
    f1=nc.Dataset(infile,'r')
    if len(f1.dimensions) >3:
        print 'the input file has more than 3 dimensions!'

    fdimdic=f1.dimensions
    fdimname=[]
    for i in fdimdic.keys():
        fdimname.append(str(i))
    for i in fdimname:
        tempdimvar=f1.dimensions[i]
        if tempdimvar.isunlimited():
            globe_unlimdimname=i
    fdimname.remove(globe_unlimdimname)
    globe_limdimname1=fdimname[0]
    globe_limdimname2=fdimname[1]

    for name in f1.variables.keys():
        tempvar_full=f1.variables[str(name)]
        if len(tempvar_full.dimensions)==1:
            if str(tempvar_full.dimensions[0])==globe_unlimdimname:
                unlimdimvar_name=str(name)
            elif str(tempvar_full.dimensions[0])==globe_limdimname1:
                limdimvar1_name=str(name)
            elif str(tempvar_full.dimensions[0])==globe_limdimname2:
                limdimvar2_name=str(name)

    limdimvar1_full=f1.variables[limdimvar1_name]
    limdimvar2_full=f1.variables[limdimvar2_name]
    unlimdimvar_full=f1.variables[unlimdimvar_name]
    
    print 'Begin to write to netcdf file ',outfile
    rootgrp=nc.Dataset(outfile,'w',format='NETCDF3_CLASSIC')
    
    #creat dimension
    lat=rootgrp.createDimension(globe_limdimname1,len(limdimvar1_full))
    lon=rootgrp.createDimension(globe_limdimname2,len(limdimvar2_full))
    rootgrp.createDimension(globe_unlimdimname,None)
    print 'Dimensions ',globe_limdimname1,globe_limdimname2,globe_unlimdimname,' created'
    
    #creat dimension variable
    lat=rootgrp.createVariable(str(limdimvar1_full._name),limdimvar1_full.dtype,(globe_limdimname1,))
    lon=rootgrp.createVariable(str(limdimvar2_full._name),limdimvar2_full.dtype,(globe_limdimname2,))
    if hasattr(limdimvar1_full,'long_name') and hasattr(limdimvar2_full,'long_name'):
        lat.long_name = limdimvar1_full.long_name.encode()
        lon.long_name = limdimvar2_full.long_name.encode()
    lat.units=limdimvar1_full.units.encode()
    lon.units=limdimvar2_full.units.encode()
    lat[:]=limdimvar1_full[:].copy()
    lon[:]=limdimvar2_full[:].copy()
    time=rootgrp.createVariable(unlimdimvar_full._name,unlimdimvar_full.dtype,(globe_unlimdimname,))
    time[:]=np.arange(1,len(unlimdimvar_full[:])/12+1)
    time.long_name = 'time'
    print 'Dimension variables ','--'+str(limdimvar1_full._name)+'--','--'+str(limdimvar2_full._name)+'--','--time-- created'

    varlist_all=[name.encode() for name in f1.variables.keys()]
    varlist_all.remove(str(limdimvar1_full._name))
    varlist_all.remove(str(limdimvar2_full._name))
    varlist_all.remove(str(unlimdimvar_full._name))
    if varlist=='all':
        varlist2=varlist_all
    else:
        varlist2=varlist

    #creat ordinary variable
    for varname in varlist2:
        var_full=f1.variables[varname]
        if str(var_full.dimensions[0])!=globe_unlimdimname:
            print 'the time dimension is not the first dimension for varialbe --',varname,'--!'
        var_value=copy.deepcopy(var_full[:])
        if mode==None:
            raise ValueError('please specify one mode!')
        elif mode=='mean':
            vardata=mathex.m2ymean(var_value)
        elif mode=='sum':
            vardata=mathex.m2ysum(var_value)
        else:
            raise ValueError('only mean or sum can be provided!')
        ba=rootgrp.createVariable(str(var_full._name),var_full.dtype,(str(var_full.dimensions[0]),str(var_full.dimensions[1]),str(var_full.dimensions[2]),))
        if hasattr(var_full,'long_name'):
            ba.long_name=var_full.long_name
        try:
            ba.units=var_full.units
        except AttributeError:
            pass
        if vardata.ndim==2:
            ba[0,:,:]=vardata
        else:
            ba[:]=vardata
        if np.ma.isMA(vardata):
            ba.missing_value=vardata.fill_value
        else:
            ba.missing_value=float(-9999)
        print 'Variable  --',str(var_full._name),'--  is fed into the file ','--',outfile,'--'

    if hasattr(f1,'history'):
        rootgrp.history='year value by applying method month-to-year '+mode+' from file '+'--'+infile+'--\n'+f1.history.encode()
    else:
        rootgrp.history='year value by applying method month-to-year '+mode+' from file '+'--'+infile+'--'
    f1.close()
    rootgrp.close()



def ncmask(infile,outfile,llt=('lat_name','lon_name','rec_name'),mask=None):
    """
    This is only an easy application of pb.ncread and gnc.dic2nc function, can only handle varialbs with 3 dimension.
    Note:
         1. The variable name for latitude,longitude and time must be explicitly stated.
         2. latitude and longitude must be 1D array but NOT meshgrid
    """
    if llt==('lat_name','lon_name','rec_name'):
        raise ValueError('lat,lon,record variable name must be explicitly stated.')
    else:
        d0,d1=pb.ncreadg(infile)
        varlist=d1.__dict__.keys()
        for dimvar in llt:
            varlist.remove(dimvar)
        #creat output varlist
        outvarlist=[]
        for var in varlist:
            data=mathex.ndarray_apply_mask(d1.__dict__[var],mask)
            varfull=d0.__dict__[var]
            varfull_name=varfull._name
            if hasattr(varfull,'long_name'):
                varfull_longname=varfull.long_name
            else:
                varfull_longname=varfull._name
            if hasattr(varfull,'units'):
                varfull_units=varfull.units
            else:
                varfull_units='no units'
            vartuple=(varfull_name,varfull_longname,varfull_units,data)
            outvarlist.append(vartuple)
        #get length of record variable
        try:
            reclen=len(d1.__dict__[llt[2]])
        except TypeError:
            reclen=1
        dic2nc(outfile,latvar=d1.__dict__[llt[0]],lonvar=d1.__dict__[llt[1]],timevar=(d0.__dict__[llt[2]]._name,reclen),\
        vardic=outvarlist,histtxt='created from file '+infile)

def ncfilemap(infile,latvarname=None,lonvarname=None,mapvarname=None,mapdim=None,agremode=None,pyfunc=None,mask=None,unit=None,title=None,\
              projection='cyl',mapbound='all',gridstep=(30,30),shift=False,cmap=None,map_threshold=None,colorbarlabel=None,levels=None,\
              data_transform=False,ax=None,\
              colorbardic=None):
    """
    OneLineInfo: This is a simple wrap of ncdatamap
    Purpose: plot map directly from a nc file.
    Arguments: cf. ncdatamap arguments
    Returns: 
        if ax==None:
            return d0,d1,fig,axt,m,cbar,mapvar
        else:
            return d0,d1,m,cbar,mapvar
        ---
        d0,d1 --> d0,d1=pb.ncreadg(infile)
        mapvar --> final data used for plotting the map.
    See also bmap.contourfmap2
    Example:
        see testdata/gnc_ncmap_eg.py for examples.
    """
    d0,d1=pb.ncreadg(infile)
    if ax==None:
        fig,axt,m,cbar,mapvar=ncdatamap(d0,d1,latvarname=latvarname,lonvarname=lonvarname,mapvarname=mapvarname,mapdim=mapdim,agremode=agremode,pyfunc=pyfunc,\
                                        mask=mask,unit=unit,title=title,\
                                        projection=projection,mapbound=mapbound,gridstep=gridstep,shift=shift,cmap=cmap,map_threshold=map_threshold,\
                                        colorbarlabel=colorbarlabel,levels=levels,data_transform=data_transform,ax=ax,colorbardic=colorbardic)
        return d0,d1,fig,axt,m,cbar,mapvar
    else:
        m,cbar,mapvar=ncdatamap(d0,d1,latvarname=latvarname,lonvarname=lonvarname,mapvarname=mapvarname,mapdim=mapdim,agremode=agremode,pyfunc=pyfunc,\
                                mask=mask,unit=unit,projection=projection,mapbound=mapbound,gridstep=gridstep,shift=shift,cmap=cmap,map_threshold=map_threshold,\
                                colorbarlabel=colorbarlabel,levels=levels,data_transform=data_transform,ax=ax,colorbardic=colorbardic)
        return d0,d1,m,cbar,mapvar


def ncdatamap(d0,d1,latvarname=None,lonvarname=None,mapvarname=None,forcedata=None,mapdim=None,agremode=None,pyfunc=None,mask=None,mask_value=None,unit=None,title=None,\
             projection='cyl',mapbound='all',gridstep=(30,30),shift=False,cmap=None,map_threshold=None,colorbarlabel=None,levels=None,data_transform=False,ax=None,\
             colorbardic={}):
    """
    Purpose: plot map directly from d0,d1=pb.ncreadg('file.nc') object.
    Arguments:
        d0,d1 --> generated by d0,d1=pb.ncreadg(infile)
        latvarname & lonvarname --> Dimension variables in the file; by default the function tries to use ('lat','lon'),or ('latitude','longtitude') pair.
        mapvarname --> varialbe name for mapping. variable dimension should be strictly as (timedim,lat,lon); 4dim varialbe are not handled currently.
        forcedata --> 
            1. when forcedata is ndarray, mapvarname will not be used to retrieve data from d1 but only for displaying variable name purpose. Instead,
               forcedata will be used as data for plotting. This can be used in case data is very big to increase speed. 
            2. In this case, mapvarname will still be used to retrieve variable units from d0.
        unit --> None for use default units inside file; 'unit' to overwrite default; False for no units display.
        title --> None for use default var lont_name inside file; 'I am title' to overwrite default; False for no title display.
        mapdim --> when variable has 3 dim,eg. set mapdim=2 to plot only for the 3rd (mapdim 0-based) data in the time dimension.
        agremode --> 'sum' or 'mean', m2y transformation apply to data prior to mapping. None for ploting original (monthly) data. When agremode!=None, 
            mapdim means the selected year. 'fsum' or 'fmean' apply to the first axis of data.
        pyfunc --> 
            1. a function that will be applied on data. if it's only a number, then data will be multiplied by given number; otherwise one can use pyfunc = 
               lambda x : x*2+3 to define a function. 
            2. pyfunc is used after applying the agremode operation.
            3. Note when using pyfunc, map_threshold will be applied to data after pyfunc transformation as aim of map_threshold is mainly for mapping prupose.
        mask --> a 2dim boolean array
        levels --> levels are used in a way after pyfunc application.
        ax --> axex for which to plot.
        other args --> the same as in function bmap.contourfmap2.
    Returns: 
        if ax==None:
            return fig,ax,m,cbar,mapvar
        else:
            return m,cbar,mapvar
        mapvar --> final data used for plotting the map.
    See also bmap.contourfmap2
    Example:
        see testdata/gnc_ncmap_eg.py for examples. (examples last tested 24/May/2012)
    """
    if latvarname==None and lonvarname==None:
        if 'lat' in d0.__dict__.keys() and 'lon' in d0.__dict__.keys():
            latvar=d0.__dict__['lat'][:]
            lonvar=d0.__dict__['lon'][:]
        elif 'latitude' in d0.__dict__.keys() and 'longitude' in d0.__dict__.keys():
            latvar=d0.__dict__['latitude']
            lonvar=d0.__dict__['longitude']
        else:
            raise ValueError('Default lat and lon names are not found in the data, please speicify them and add into default list')
    else:
        latvar=d0.__dict__[latvarname][:]
        lonvar=d0.__dict__[lonvarname][:]
    #extract 2dim lat,lon to 1dim.
    if latvar.ndim==2:
        latvar=latvar[:,0]
    if lonvar.ndim==2:
        lonvar=lonvar[0,:]

    #read mapvar data and prepare for mapping
    if mapvarname==None:
        if forcedata == None:
            raise ValueError("mapvarname and forcedata both as None!")
        else:
            mapvar = forcedata
        #raise ValueError('please provide a varialbe name or data')
    else:
        if isinstance(mapvarname,str):
            if forcedata==None:
                mapvar=d1.__dict__[mapvarname]
            else:
                mapvar=forcedata
        else:
            raise ValueError('mapvarname must be string type')
    if mapvar.ndim==2 and (mapdim!=None or agremode!=None):
        raise ValueError('{0} has only 2 valid dimension but mapdim or agremode is not None'.format(mapvarname))
    elif mapvar.ndim==3:
        if agremode==None and mapdim==None:
            raise ValueError('{0} has 3 valid dimension, cannot leave both mapdim and agremode as None'.format(mapvarname))
        elif agremode!=None:
            #make sum or mean of monthly data
            if agremode=='sum':
                mapvar=mathex.m2ysum(mapvar)
            elif agremode=='mean':
                mapvar=mathex.m2ymean(mapvar)
            elif agremode=='fsum':
                mapvar=np.ma.sum(mapvar,axis=0)
            elif agremode=='fmean':
                mapvar=np.ma.mean(mapvar,axis=0)

            #extract data for only specified dimension
            if mapdim!=None:
                if mapvar.ndim==2: #the original mapvar has 3 dim with first dim of 12.
                    raise ValueError ('{0} has 3 valid dimension with firest dim size as 12, cannot specify agremode and mapdim simultaneously'.format(mapvarname))
                else: #mapvar.ndim==3
                    mapvar=mapvar[mapdim,:,:]
            #error handling when mapdim==None
            else:
                if mapvar.ndim==2: #the original mapvar has 3 dim with first dim of 12.
                    pass
                else: #mapvar.ndim==3
                    raise ValueError('the dimension after data transformation is {0}, must specify mapdim to plot'.format(mapvar.shape))
        #plot specified dimension for data without transformation
        else:
            mapvar=mapvar[mapdim,:,:]
    elif mapvar.ndim==4:
        raise ValueError('cannot handle 4 dimensional data')
    #apply the defined function on data
    if pyfunc!=None:
        if isfunction(pyfunc):
            mapvar=pyfunc(mapvar)
        elif isinstance(pyfunc,list) and isfunction(pyfunc[0]):
            for subpyfunc in pyfunc:
                mapvar=subpyfunc(mapvar)
        else:
            mapvar=mapvar*pyfunc

    #apply mask
    if mask!=None:
        mapvar=mathex.ndarray_apply_mask(mapvar,mask)

    if mask_value != None:
        mapvar=np.ma.masked_equal(mapvar,mask_value)

    #finally, if lat is provide in increasing sequence, flip over the data.
    if latvar[0]<latvar[-1]:
        latvar=latvar[::-1]
        mapvar=np.flipud(mapvar)
    #make the mapping
    if ax==None:
        fig,axt,m,cbar=bmap.contourfmap2(latvar,lonvar,mapvar,projection=projection,mapbound=mapbound,gridstep=gridstep,shift=shift,cmap=cmap,\
                                        map_threshold=map_threshold,\
                                        colorbarlabel=colorbarlabel,levels=levels,data_transform=data_transform,ax=None,colorbardic=colorbardic)
    else:
        axt=ax
        m,cbar=bmap.contourfmap2(latvar,lonvar,mapvar,projection=projection,mapbound=mapbound,gridstep=gridstep,shift=shift,cmap=cmap,\
                                        map_threshold=map_threshold,\
                                        colorbarlabel=colorbarlabel,levels=levels,data_transform=data_transform,ax=ax,colorbardic=colorbardic)

    mapvar_full=d0.__dict__[mapvarname]
    #fuction to retreive title or unit from either inside data or
    #by external forcing
    def retrieve_external_default(external_var,attribute):
        if external_var!=None:
            outvar=external_var
        elif hasattr(mapvar_full,attribute):
            if external_var==False:
                outvar=None
            else:
                outvar=mapvar_full.getncattr(attribute)
        else:
            outvar=None
        return outvar
    #retrieve title or unit
    map_unit=retrieve_external_default(unit,'units')
    map_title=retrieve_external_default(title,'long_name')
    try:
        agre_title_complement='[yearly '+agremode+']'
    except:
        agre_title_complement=None

    #function handling ax title
    def set_title_unit(ax,title=None,unit=None,agre_title_complement=None):
        try:
            title_unit=title+('\n'+unit)
        except TypeError:
            try:
                title_unit=title
                title_agre=agre_title_complement+' '+title
            except TypeError:
                pass
        finally:
            try:
                title_full=agre_title_complement+' '+title_unit
            except TypeError:
                title_full=title_unit
        if title_full!=None:
            ax.set_title(title_full)
        else:
            pass
    set_title_unit(axt,map_title,map_unit,agre_title_complement)
    if ax==None:
        return fig,axt,m,cbar,mapvar
    else:
        return m,cbar,mapvar

def ncreadg(filename):
    """
    Purpose: read a .nc file using netCDF4 package and read the data into netCDF4.Variable object
    Definition: ncread(filename)
    Arguments:
        file--> file name
    Return: return a list of ncdata object; the 1st one contains original nctCDF4 objects and the 2nd one contains data with duplicate dimensions removed.
    Note:
        1. This is for general purpose

    Example:
        >>> data=g.ncread('cru1901.nc')
        >>> data
          <g.ncdata object at 0x2b1e4d0>
        >>> data.t2m
          <netCDF4.Variable object at 0x2b20c50>

    """
    f1=nc.Dataset(filename,mode="r")
    datanew=ncdata()
    datanew2=ncdata()
    for var in f1.variables.keys():
        var=str(var)
        datanew.__dict__[var]=f1.variables[var]
        datanew2.__dict__[var]=Remove_dupdim(f1.variables[var][:])
    return [datanew,datanew2]
    f1.close()


def ncfile_dim_check(infile,outfile,varlist='all',reginterval=None):
    """
    """
    f1=nc.Dataset(infile,'r')
    if len(f1.dimensions) >3:
        raise ValueError('the input file has more than 3 dimensions!')

    #identity dimension varialbe and unlimited variable
    fdimdic=f1.dimensions
    fdimname=[]
    for i in fdimdic.keys():
        fdimname.append(str(i))
    for i in fdimname:
        tempdimvar=f1.dimensions[i]
        if tempdimvar.isunlimited():
            globe_unlimdimname=i
    fdimname.remove(globe_unlimdimname)
    globe_limdimname1=fdimname[0]
    globe_limdimname2=fdimname[1]

    for name in f1.variables.keys():
        tempvar_full=f1.variables[str(name)]
        if len(tempvar_full.dimensions)==1:
            if str(tempvar_full.dimensions[0])==globe_unlimdimname:
                unlimdimvar_name=str(name)
            elif str(tempvar_full.dimensions[0])==globe_limdimname1:
                limdimvar1_name=str(name)
            elif str(tempvar_full.dimensions[0])==globe_limdimname2:
                limdimvar2_name=str(name)

    limdimvar1_full=f1.variables[limdimvar1_name]
    limdimvar2_full=f1.variables[limdimvar2_name]
    unlimdimvar_full=f1.variables[unlimdimvar_name]

class Ncdata(object):
    """
    NCdata is object facilitating maping and exploring nc file.
    """
    _timedim_candidate_list=['time_counter','time','tstep']
    _timevar_name_can_list=['time_counter','time']
    _default_pftdim='PFT' #not used

    def __init__(self,filename,latlon_dim_name=None):
        """
        Parameters:
        -----------
        latlon_dim_name: a tuple giving the names of lat/lon dimension names of the file (eg. ('latdimname','londimname')), if it's None, they will be guessed from the
            default candidate list.
        """
        self.filename=filename
        ncf=nc.Dataset(filename)
        #get the unlimited dim name
        timedim_detect = False
        for timedim_candidate in Ncdata._timedim_candidate_list:
            if timedim_candidate in ncf.dimensions.keys():
                tempdim=ncf.dimensions[timedim_candidate]
                if tempdim.isunlimited():
                    self.unlimited_dimname=timedim_candidate
                    self.unlimited_dimlen=len(tempdim)
                    timedim_detect = True
                    break
                else:
                    raise ValueError("the dimension {0} exists but is not unlimited dimension".format(timedim_candidate))
            else:
                pass
        if not timedim_detect:
            print "unlimited_dimname not found, please make sure the file has no unlimited dimension"
            self.unlimited_dimname = None

            #raise ValueError("please expand the time dimension name candidates list")

        self.dimensions = ncf.dimensions
        #get the global attributes.
        self.global_attributes = dict(ncf.__dict__)
        ncf.close()

        #guess the lat/lon dimension names
        if latlon_dim_name == None:
            if 'lat' in self.dimensions and 'lon' in self.dimensions:
                self.latdim_name = 'lat'
                self.londim_name = 'lon'
            elif 'latdim' in self.dimensions and 'londim' in self.dimensions:
                self.latdim_name = 'latdim'
                self.londim_name = 'londim'
            elif 'latitude' in self.dimensions and 'longitude' in self.dimensions:
                self.latdim_name = 'latitude'
                self.londim_name = 'longitude'
            elif 'x' in self.dimensions and 'y' in self.dimensions:
                self.latdim_name = 'y'
                self.londim_name = 'x'
            else:
                raise ValueError("lat/lon dimension names could not be guessed, please either provide the latlon_dim_name or expand default guess list")
        else:
            self.latdim_name = latlon_dim_name[0]
            self.londim_name = latlon_dim_name[1]

        #build the dimension names list
        self.dim_name_list = [self.unlimited_dimname, self.latdim_name, self.londim_name]
        if 'PFT' in self.dimensions:
            self.dim_name_list = [self.unlimited_dimname, 'PFT', self.latdim_name, self.londim_name]


        #read in data
        self.d0,self.d1=pb.ncreadg(filename)

        #set default latvar,lonvar name and retrieve lonvar/latvar values.
        if 'lat' in self.d1.__dict__.keys() and 'lon' in self.d1.__dict__.keys():
            self.latvar_name='lat'
            self.lonvar_name='lon'
        elif 'latitude' in self.d1.__dict__.keys() and 'longitude' in self.d1.__dict__.keys():
            self.latvar_name='latitude'
            self.lonvar_name='longitude'
        elif 'nav_lat' in self.d1.__dict__.keys() and 'nav_lon' in self.d1.__dict__.keys():
            self.latvar_name='nav_lat'
            self.lonvar_name='nav_lon'
        elif 'latvar' in self.d1.__dict__.keys() and 'lonvar' in self.d1.__dict__.keys():
            self.latvar_name='latvar'
            self.lonvar_name='lonvar'
        else:
            raise ValueError('Default lat and lon names are not found in the data, please specify by add_latvar_lonvar_name, or update default list')

        if np.ndim(self.d1.__dict__[self.latvar_name]) in [0,1]:
            self.latvar = self.d1.__dict__[self.latvar_name]
        elif np.ndim(self.d1.__dict__[self.latvar_name]) == 2:
            self.latvar = self.d1.__dict__[self.latvar_name][:,0]
        else:
            raise ValueError("the lat variable ndim is {0}".format(np.ndim(self.d1.__dict__[self.latvar_name])))

        if np.ndim(self.d1.__dict__[self.lonvar_name]) in [0,1]:
            self.lonvar = self.d1.__dict__[self.lonvar_name]
        elif np.ndim(self.d1.__dict__[self.lonvar_name]) == 2:
            self.lonvar = self.d1.__dict__[self.lonvar_name][0,:]
        else:
            raise ValueError("the lon variable ndim is {0}".format(np.ndim(self.d1.__dict__[self.lonvar_name])))

        self.geo_limit={'lat':(np.min(self.latvar),np.max(self.latvar)), 'lon':(np.min(self.lonvar),np.max(self.lonvar))}

        #This is for compatiblity with previous version.
        self.lat_name = self.latvar_name
        self.lon_name = self.lonvar_name
        self.lat = self.latvar
        self.lon = self.lonvar

        #get time varialbe
        timevar_detect = False
        for timevar_name_can in Ncdata._timevar_name_can_list:
            if timevar_name_can in self.d1.__dict__.keys():
                if self.d0.__dict__[timevar_name_can].dimensions[0] == self.unlimited_dimname:
                    self.timevar_name = timevar_name_can
                    timevar_detect = True
                    break
                else:
                    raise ValueError("candidate timevar {0} detected but its dimensions is not {1}".format(timevar_name_can,self.unlimited_dimname))
            else:
                pass
        if not timevar_detect:
            print "timevar not found. please make sure there is not time var in the file"
            #raise ValueError("timevar not found, please expand the default candidate list.")
        if timevar_detect:
            self.timevar = self.d1.__dict__[self.timevar_name]
        else:
            self.timevar = None
            self.timevar_name = None

        #build the dimension varialbe list
        self.dimvar_name_list = [self.lonvar_name, self.latvar_name, self.timevar_name]
        self.varlist = pb.StringListAnotB(self.d1.__dict__.keys(),self.dimvar_name_list)

        #get PFT variable
        if 'PFT' in self.d1.__dict__.keys():
            self.pftvar_name = 'PFT'
            self.pftvar = self.d1.__dict__[self.pftvar_name]
            self.dimvar_name_list.append(self.pftvar_name)

    def add_latvar_lonvar_name(self,latvar_name,lonvar_name):
        self.latvar_name=latvar_name
        self.lonvar_name=lonvar_name

    def get_pftsum(self,varlist=None,veget_npindex=np.s_[:]):
        d0=self.d0
        d2=g.ncdata()
        d1=self.d1
        print "*******PFT VEGET_MAX weighted sum begin******"
        for var in d1.__dict__.keys():
            if var!='VEGET_MAX':
                #4-dim variable before squeeze
                if d0.__dict__[var].ndim==4:
                    if 'PFT' not in d0.__dict__[var].dimensions:
                        print "var '{0}' has 4 dimensions but PFT is not one of the dimensions".format(var)
                    #'PFT' is one dimension
                    else:
                        #if 'PFT' is the second dimension
                        if d0.__dict__[var].dimensions[1]=='PFT':
                            #first/time dim lenght is not 1
                            if d1.__dict__[var].ndim == 4:
                                #veget_max=np.ma.masked_equal(d1.__dict__['VEGET_MAX'],0.)
                                #temp=pb.MaskArrayByNan(d1.__dict__[var])*veget_max
                                veget_max=d1.VEGET_MAX[veget_npindex]
                                vardata = d1.__dict__[var][veget_npindex]
                                temp=vardata*veget_max
                                temppftsum=np.ma.sum(temp,axis=1)
                                d2.__dict__[var]=temppftsum
                                print '{0} treated'.format(var)
                            #first/time dim length is 1
                            elif d1.__dict__[var].ndim == 3:
                                #veget_max=np.ma.masked_equal(d1.__dict__['VEGET_MAX'],0.)
                                #temp=pb.MaskArrayByNan(d1.__dict__[var])*veget_max
                                veget_max=d1.VEGET_MAX[veget_npindex]
                                vardata = d1.__dict__[var][veget_npindex]
                                temp=vardata*veget_max
                                temppftsum=np.ma.sum(temp,axis=0)
                                d2.__dict__[var]=temppftsum
                                print '{0} treated'.format(var)
                            else:
                                print "original 4-dim var '{0}' has less than 3 dimensions after being squeezed".format(var)

                        else:
                            print "var '{0}' has 4 dimension but the second dimension is not PFT".format(var)
                #3-dim variable before squeeze
                elif d0.__dict__[var].ndim==3:
                    #if PFT is not one of the remaining dimensions, probably this variable has already been PFT veget weigthed summed.
                    if 'PFT' not in d0.__dict__[var].dimensions:
                        d2.__dict__[var]=d1.__dict__[var]
                    else:
                        print "original var '{0}' has 3 dim but PFT is still one dimension".format(var)
                else:
                    pass

        try:
            d2.__dict__['Areas']=d1.__dict__['Areas']
        except KeyError:
            print "Warning! Areas not in the history file!"
        self.pftsum=d2
        print "*******END of PFT sum******"

    def get_spa(self,pftop=True):
        """
        Get the PFT weighed spatial mean and sum, only non-masked values will be considered.
        """
        if not hasattr(self,'pftsum'):
            if pftop==True:
                self.get_pftsum()
                d2=self.pftsum
            else:
                d2=self.d1
        else:
            d2=self.pftsum

        d3=g.ncdata()  #sum
        d4=g.ncdata()  #mean
        for var in d2.__dict__.keys():
            if var!='Areas':
                temppftsum=d2.__dict__[var]
                try:
                    temparea=temppftsum*d2.__dict__['Areas']
                except AttributeError:
                    raise KeyError("Error! Areas not in the file and cannot carry out spatial sum!")
                tempspa_sum=np.ma.sum(np.ma.sum(temparea,axis=-1),axis=-1)
                d3.__dict__[var]=tempspa_sum
                tempspa_mean=tempspa_sum/np.ma.sum(np.ma.sum(d2.__dict__['Areas'],axis=-1), axis=-1)
                d4.__dict__[var]=tempspa_mean
            else:
                pass
        self.spasum=d3
        self.spamean=d4

    def combine_vars(self,varlist,pftsum=False):
        """
        combine variales into a ndarray by adding a new dimension.
        """
        out_array = []
        if pftsum == False:
            data = self.d1
        else:
            data = self.pftsum
        for varname in varlist:
            out_array.append(data.__dict__[varname])
        return np.ma.array(out_array)

    def list_var(self,keyword=None):
        """
        supply keyword for show var names matching the pattern.
        supply "excludedim" to keyword to exclude the dim varnames in the
        output list.
        """
        if keyword==None:
            return self.d1.__dict__.keys()
        else:
            if keyword == 'excludedim':
                return pb.StringListAnotB(self.d1.__dict__.keys(),
                                          self.dimvar_name_list+['Areas','VEGET_MAX'])
            else:
                return pb.FilterStringList(keyword,self.d1.__dict__.keys())

    def list_var_attr(self,varname):
        """
        list var attributes.
        """
        for attr_name in self.d0.__dict__[varname].__dict__.keys():
            print "{0} : {1}".format(attr_name,self.d0.__dict__[varname].__dict__[attr_name])
        print "Reduced dimension: {0}".format(self.d1.__dict__[varname].shape)

    def map(self,mapvarname=None,forcedata=None,mapdim=None,agremode=None,pyfunc=None,mask=None,unit=None,title=None,pftsum=False,mask_value=None,\
             projection='cyl',mapbound='all',gridstep=(30,30),shift=False,cmap=None,map_threshold=None,colorbarlabel=None,levels=None,data_transform=False,ax=None,\
             colorbardic={}):
        """
        This is an implementation of ncdatamap
        """
        if pftsum==False:
            d0,d1=self.d0,self.d1
        else:
            d0,d1=self.d0,self.pftsum

        mlist= ncdatamap(d0,d1,latvarname=self.latvar_name,lonvarname=self.lonvar_name,mapvarname=mapvarname,forcedata=forcedata, mapdim=mapdim, \
                         agremode=agremode,pyfunc=pyfunc,\
                         mask=mask,mask_value=mask_value,unit=unit,title=title,\
                         projection=projection,mapbound=mapbound,gridstep=gridstep,shift=shift,cmap=cmap,map_threshold=map_threshold,\
                         colorbarlabel=colorbarlabel,levels=levels,data_transform=data_transform,ax=ax,colorbardic=colorbardic)
        self.m = mlist[-3]
        self.cbar = mlist[-2]
        self.drawdata = mlist[-1]

    def imshowmap(self,varname,forcedata=None,pftsum=False,ax=None,projection='cyl',mapbound='all',gridstep=(30,30),shift=False,colorbar=True,
                  colorbarlabel=None,*args,**kwargs):
        """
        A temporary wrapper of bamp.imshowmap
        """
        if pftsum==False:
            d1=self.d1
        else:
            d1=self.pftsum

        if forcedata != None:
            data = forcedata
        else:
            data = d1.__dict__[varname]
        m,cs,cbar = bmap.imshowmap(self.lat,self.lon,data,ax=ax,projection=projection,mapbound=mapbound,gridstep=gridstep,shift=shift,colorbar=colorbar,
                                   colorbarlabel=colorbarlabel,*args,**kwargs)
        self.m=m
        self.cbar=cbar
        self.cs=cs

    def scatter(self,lon,lat,*args,**kwargs):
        """
        Draw scatter points on the map.

        Parameters:
        -----------
        lon/lat: single value for 1D ndarray.
        index keyword: lon/lat could be the indexes when index == True,
            when making scatter plots, index is removed from the keys
            before it's passed to scatter function.
        """
        index = kwargs.get('index',False)
        if index == True:
            lon = self.lon[lon]
            lat = self.lat[lat]
            del kwargs['index']
        elif index == False:
            pass
        else:
            raise TypeError("keyword index must boolean type!")

        x,y = self.m(lon,lat)
        self.m.scatter(x,y,*args,**kwargs)

    def add_Rectangle(self,(lat1,lon1),(lat2,lon2),index=False,
                      **kwargs):
        """
        Add a rectangle of by specifing (lat1,lon1) and (lat2,lon2).

        Parameters:
        -----------
        (lat1,lon1): lowerleft coordinate
        (lat2,lon2): upperright coordinate
        """
        if index:
            lat1,lat2 = self.lat[lat1],self.lat[lat2]
            lon1,lon2 = self.lon[lon1],self.lon[lon2]
        else:
            pass
        (x1,x2),(y1,y2) = self.m([lon1,lon2],[lat1,lat2])
        rec = mat.patches.Rectangle((x1,y1),x2-x1,y2-y1,**kwargs)
        self.m.ax.add_patch(rec)
        return rec

    def add_text(self,lat,lon,s,fontdict=None,**kwargs):
        """
        Add text s on the position (lat,lon)
        """
        x,y = self.m(lon,lat)
        self.m.ax.text(x,y,s,fontdict=fontdict,**kwargs)

    def add_Rectangle_list_coordinates(self,coordlist,textlist=None,fontdict=None,textkw={},**kwargs):
        """
        Add a list of [(lat1,lon1),(lat2,lon2)] to add a series of rectangles at one time.

        Parameters:
        -----------
        coordlist: a nested list with [(lat1,lon1),(lat2,lon2)] as its members.
        """
        for i,coord in enumerate(coordlist):
            self.add_Rectangle(coord[0],coord[1],**kwargs)
            if textlist != None:
                lat = (coord[0][0]+coord[1][0])/2.
                lon = (coord[0][1]+coord[1][1])/2.
                self.add_text(lat,lon,textlist[i],fontdict=fontdict,**textkw)

    def add_Rectangle_list_by_dataframe(self,dataframe,label=False,
                                        fontdict=None,textkw={},**kwargs):
        """
        Add rectangle list by using dataframe.

        Notes:
        ------
        1. the dataframe should have "lat1,lat2,lon1,lon2,region" as column
            names.
        """
        for name in ['lat1','lat2','lon1','lon2','region']:
            if name not in dataframe.columns:
                raise ValueError("{0} not a column name of dataframe").format(name)
        region_name_list = []
        coordlist = []
        for index,row in dataframe.iterrows():
            region_name_list.append(row['region'])
            coordlist.append([(row['lat1'],row['lon1']),(row['lat2'],row['lon2'])])
        if label == False:
            textlist = None
        else:
            textlist = region_name_list
        self.add_Rectangle_list_coordinates(coordlist,textlist=textlist,
                        fontdict=fontdict,textkw=textkw,**kwargs)


    def Get_PointValue(self,var,(vlat,vlon)):
        return pb.Get_PointValue(self.d0,var,(self.latvar_name,vlat),(self.lonvar_name,vlon))

    def find_index_by_vertex(self,(vlat1,vlat2),(vlon1,vlon2)):
        """
        Purpose: find the index specified by (vla1,vlat2),(vlon1,vlon2)
        Return: (lon_index_min, lon_index_max, lat_index_min, lat_index_max)
        Note: This is a direct application of gnc.find_index_by_vertex.
        """
        return find_index_by_vertex(self.lonvar, self.latvar, (vlon1,vlon2), (vlat1,vlat2))

    def find_index_latlon_each_point_in_vertex(self,rlat=None,rlon=None):
        """
        Return a list of nested tuples, the nested tuple is like:
            ((index_lat,index_lon),(vlat,vlon)), suppose var is of shape
            (time,lat,lon), var[time,index_lat,index_lon] will allow
            to retrieve var value corresponding to (vlat,vlon)

        parameters:
        -----------
        rlat: range of lat, (valt1,vlat2)
        rlon: range of lon, (vlon1,vlon2)

        Notes:
        ------
        1. Suppose the rlat include m gridcells and rlon include n gridcells,
            the result will be mXn length list and it will allow to retrieve
            the variable values for each of the point which are within
            the rectangle given by rlat,rlon.

        See also,
        ---------
        Add_Single_Var_Pdata
        """
        rlat,rlon = self._return_defautl_latlon_range(rlat,rlon)
        (lon_index_min, lon_index_max, lat_index_min, lat_index_max) = \
            self.find_index_by_vertex(rlat,rlon)
        index_latlon_tuple_list = []
        for inlat in range(lat_index_min, lat_index_max+1):
            for inlon in range(lon_index_min, lon_index_max+1):
                index_latlon_tuple_list.append((
                                    (inlat,inlon),
                                    (self.lat[inlat],self.lon[inlon])
                                    ))
        return index_latlon_tuple_list


    def Get_GridValue(self,var,(vlat1,vlat2),(vlon1,vlon2), pftsum=False):
        (lon_index_min, lon_index_max, lat_index_min, lat_index_max) = find_index_by_vertex(self.lonvar, self.latvar, (vlon1,vlon2), (vlat1,vlat2))
        if pftsum:
            try:
                return self.pftsum.__dict__[var][..., lat_index_min:lat_index_max+1, lon_index_min:lon_index_max+1]
            except AttributeError:
                raise ValueError("please do pftsum operation first!")
        else:
            return self.d1.__dict__[var][..., lat_index_min:lat_index_max+1, lon_index_min:lon_index_max+1]
        #return pb.Get_GridValue(self.d0,var,(self.lat_name,vlat1,vlat2),(self.lon_name,vlon1,vlon2))


    #def check_longname_unit()



    def Plot_PointValue(self,var,(vlat,vlon),ax=None,ylab=False,pyfunc=None,pftsum=False,**kwargs):
        if ax==None:
            ax=plt.gca()

        data=self.Get_PointValue(var,(vlat,vlon))
        if pftsum==True:
            veget=np.ma.masked_array(self.Get_PointValue('VEGET_MAX',(vlat,vlon)),0.)
            data=np.ma.sum(data*veget,axis=1)
        if pyfunc!=None:
            if isfunction(pyfunc):
                data=pyfunc(data)
            else:
                data=data*pyfunc

        ax.plot(data,label=var,**kwargs)
        if ylab==True:
            try:
                ax.set_ylabel(self.d0.__dict__[var].getncattr('units'))
            except AttributeError:
                pass
        elif isinstance(ylab,str):
            ax.set_ylabel(ylab)
        elif ylab==False:
            pass
        else:
            raise TypeError("Incorrect ylab type")

    def Plot_PointValue_PFTsum(self,var,(vlat,vlon),ax=None,ylab=False,pyfunc=None):
        pass

    def _get_final_ncdata_by_flag(self,pftsum=False,spa=None):
        if spa=='sum':
            final_ncdata=self.spasum
        elif spa=='mean':
            final_ncdata=self.spamean
        elif spa==None:
            if pftsum==True:
                final_ncdata=self.pftsum
            else:
                final_ncdata=self.d1
        else:
            raise ValueError('''spatial operation '{0}' not expected!'''
                             .format(spa))
        return final_ncdata

    def _return_defautl_latlon_range(self,rlat=None,rlon=None):
        if rlat != None:
            pass
        else:
            rlat = self.geo_limit['lat']
        if rlon != None:
            pass
        else:
            rlon = self.geo_limit['lon']
        return (rlat,rlon)

    def Add_Vars_to_Pdata(self,varlist,npindex=np.s_[:],unit=True,
                          pd=None,pftsum=False,spa=None):
        """
        This will add the varnames in varlist to a Pdata object
            specified by npindex, note the npindex must be numpy
            index trick object (np.s_).
        """
        final_ncdata = self._get_final_ncdata_by_flag(pftsum=pftsum,
                            spa=spa)
        if pd==None:
            pd=Pdata.Pdata()
        for varname in varlist:
            data=final_ncdata.__dict__[varname][npindex]
            x=np.arange(len(data))
            pd.add_entry_noerror(x,data,varname)
            if unit==True:
                try:
                    pd.add_attr_by_tag(unit=[(varname,
                            self.d0.__dict__[varname].getncattr('units'))])
                except AttributeError:
                    pass
        return pd

    def Add_Single_Var_Pdata(self,varname,
                             rlat=None,
                             rlon=None,
                             pftsum=False,
                             pd=None,npindex=np.s_[:],):
        '''
        This is mainly to plot for each spatial point within the rectangle
            defined by rlat/rlon.

        Parameters:
        -----------
        lat: (lat1,lat2), lat range.
        lon: (lon1,lon2), lon ragne.

        Notes:
        ------
        npindex: npindex will be applied before slicing by lat/lon index.
            The data after npindex slicing should have len(rlat),len(rlon)
            as the last two dimensions.

        See also,
        ---------
        find_index_latlon_each_point_in_vertex
        '''
        if pd == None:
            pd = Pdata.Pdata()
        ydic = {}

        final_ncdata = self._get_final_ncdata_by_flag(pftsum=pftsum)
        data = final_ncdata.__dict__[varname][npindex]
        index_latlon_tuple_list = \
            self.find_index_latlon_each_point_in_vertex(rlat,rlon)
        for (inlat,inlon),lat_lon in index_latlon_tuple_list:
            ydic[str(lat_lon)] = data[...,inlat,inlon]
        pd.add_entry_sharex_noerror_by_dic(ydic)
        return pd

    def Add_Varlist_NestedPdata(self,varlist,rlat=None,rlon=None,
                           pftsum=False,npindex=np.s_[:]):
        '''
        Plot a list of variables for each point within the rectangle defined by
            rlat/rlon. Return a dictionry of Pdata.Pdata instances, with
            viriable names as keys.
        '''
        outdic = {}
        for varname in varlist:
            pd = self.Add_Single_Var_Pdata(varname,rlat=rlat,rlon=rlon,
                                           pftsum=pftsum,npindex=npindex)
            outdic[varname] = pd
        return Pdata.NestPdata(outdic)

    def Plot_Single_Var_Each_Point(self,varname,npindex=np.s_[:],
                                   pftsum=False,**legkw):
        '''
        Parameters:
        -----------
        legkw: as in plt.legend()
        '''
        pd = self.Add_Single_Var_Pdata(varname,npindex=npindex,pftsum=pftsum)
        pd.plot()
        pd.set_legend_all(taglab=True,**legkw)
        return pd

    def Plot_Varlist_Each_Point(self,varlist,npindex=np.s_[:],
                                pftsum=False,
                                legtag=None,
                                legtagseq=None,legkw={},
                                plotkw={},
                                **kwargs):
        '''
        Parameters:
        -----------
        kwargs: for Pdata.NestPdata.plot_split_parent_tag
        legkw: for plt.legend
        plotkw: for plt.plot
        '''
        npd = self.Add_Varlist_NestedPdata(varlist)
        npd.plot_split_parent_tag(plotkw=plotkw,
                                  legtag=legtag,
                                  legtagseq=legtagseq,
                                  legkw=legkw,
                                  **kwargs)
        pd0 = npd.child_pdata[varlist[0]]
        #pd0.set_legend_all(taglab=True,**legkw)

    def Add_Vars_to_Mdata(self,varlist,npindex=np.s_[:],
                          md=None,pftsum=False,spa=None):
        '''
        Add several vars to Mdata for easy mapping.
        '''
        final_ncdata = self._get_final_ncdata_by_flag(pftsum=pftsum,
                            spa=spa)
        if md==None:
            md=Pdata.Mdata()

        for varname in varlist:
            data=final_ncdata.__dict__[varname][npindex]
            md.add_tag(varname)
            md.add_array(data,varname)
            md.add_lat_lon(tag=varname,lat=self.lat,lon=self.lon)
        return md


    def Add_Single_Var_Mdata(self,varname=None,forcedata=None,
                             taglist=None,md=None,
                             pftsum=False,spa=None,npindex_list=None):
        '''
        Add single var to Mdata to explore the difference among
            different dimensions.

        TODO: allow default taglist and add a new key axis to allow select
        the dimension rather than only the first dimension.
        Parameters:
        -----------
        varname: the final data used for 'varname' should be of dim
            (tag_length,lat,lon),tag_length could be equal to dimension
            length.
        taglist: could be dimension name or other name.
        npindex_list: used to select only few dimensions for mapping.
            when npindex_list==None, the final array for varname must be
            of 3dimensions, the len(taglist) must be equal to the length
            of first dimension of the array.
        '''
        final_ncdata = self._get_final_ncdata_by_flag(pftsum=pftsum,spa=spa)
        if forcedata == None:
            maparray = final_ncdata.__dict__[varname]
        else:
            maparray = forcedata

        if md==None:
            md=Pdata.Mdata()
        if npindex_list == None:
            if len(taglist) != maparray.shape[0]:
                raise ValueError('''the length of tgalist and first
                    dimension of maparray not equal!''')
            else:
                for tag,array in zip(taglist,maparray):
                    md.add_tag(tag)
                    md.add_array_lat_lon(tag=tag,data=array,
                                         lat=self.lat,lon=self.lon)
        else:
            if len(taglist) != len(npindex_list):
                raise ValueError('''the length of tgalist not equal
                    to npindex_list''')
            else:
                for tag,nindex in zip(taglist,npindex_list):
                    md.add_tag(tag)
                    md.add_array_lat_lon(tag=tag,data=maparray[nindex],
                                         lat=self.lat,lon=self.lon)
        return md

    def break_by_region(self,varname,separation_array,
                        forcedata=None,
                        pyfunc=None,dimred_func=None,
                        pftsum=False):
        """
        Break the concerned variables into regional sum or avg or extracted array by specifying the separation_array.

        Parameters:
        -----------
        separation_array: The array that's used to separate different regions. the np.unique(separation_array) will be used as the keys for the dictionary which will be
            returned by the function.
        forcedata: used to force input data.
        pyfunc: function used to change the regional array data.
        dimred_func: functions that used to reduce the dimensions of regional array data, as long as it could be applied on numpy ndarray. eg., np.sum will get the sum
            of the regional array. If none, regional array will be returned. 
        pftsum: if True, the varname data from the pftsum will be used.

        Notes:
        ------
        1. varname could be any dimension as long as the last two dimensions are the same as separation_array.
        2. pyfunc and dimred_func work for np.ma functions

        Returns:
        --------
        region_dic: A dictionary with the region ID as keys and the region arrays or mean or sum as key values.
        """
        if pftsum == True:
            vardata = self.pftsum.__dict__[varname]
        else:
            vardata = self.d1.__dict__[varname]
        if forcedata != None:
            vardata = forcedata

        regdic={}
        reg_id_list = np.unique(separation_array)
        reg_id_list = reg_id_list.astype(int)
        for reg_id in reg_id_list:
            reg_valid = np.ma.masked_not_equal(separation_array,reg_id)
            annual_reg = mathex.ndarray_apply_mask(vardata,mask=reg_valid.mask)
            if np.any(np.isnan(annual_reg)) or np.any(np.isinf(annual_reg)):
                print "Warning! nan or inf values have been masked for variable {0} reg_id {1}".format(varname,reg_id)
                annual_reg = np.ma.masked_invalid(annual_reg)
            if pyfunc!=None:
                if isfunction(pyfunc):
                    data=pyfunc(annual_reg)
                else:
                    data=annual_reg*pyfunc
            else:
                data = annual_reg

            if dimred_func == None:
                regdic[reg_id] = data
            elif callable(dimred_func):
                if data.ndim >= 2:
                    regdic[reg_id] = dimred_func(dimred_func(data,axis=-1),axis=-1)
                else:
                    raise ValueError("strange the dimension of data is less than 2!")
            else:
                raise TypeError("dimred_func must be callable")
        return regdic

    def Plot_Vars(self,ax=None,varlist=None,npindex=np.s_[:],
                  unit=True,pftsum=False,spa=None,
                  split=False,splitkw={},
                  legkw={},
                  **plotkw):
        '''
        Plot for varlist as pd.
        '''

        pd=self.Add_Vars_to_Pdata(varlist=varlist,npindex=npindex,unit=unit,
                                  pftsum=pftsum,spa=spa)
        if not split:
            pd.plot(axes=ax,**plotkw)
            pd.set_legend_all(taglab=True,**legkw)
        else:
            pd.plot_split_axes(self,plotkw=plotkw,**splitkw)
        return pd

    def Add_Vars_to_Dict_Grid(self,varlist,grid=None,pftsum=False):
        """
        Add vars to dictionary.

        Parameters:
        -----------
        1.grid should be a tuple of (lat1,lon1,lat2,lon2)
        """
        final_ncdata = self._get_final_ncdata_by_flag(pftsum=pftsum)
        final_dict = {}
        for var in varlist:
            if grid == None:
                data = final_ncdata.__dict__[var]
            else:
                data = self.Get_GridValue(var,(grid[0],grid[2]),(grid[1],grid[3]),
                                         pftsum=pftsum)
            final_dict[var] = data
        return final_dict

def nc_get_var_value_grid(ncfile,varname,(vlat1,vlat2),(vlon1,vlon2),
                          pftsum=False):
    '''
    simple way to quickly retrieve grid value
    '''
    d = Ncdata(ncfile)
    return d.Get_GridValue(varname,(vlat1,vlat2),(vlon1,vlon2),
                    pftsum=pftsum)

def nc_get_var_value_point(ncfile,varname,(vlat,vlon)):
    '''
    simple way to quickly retrieve grid value
    '''
    d = Ncdata(ncfile)
    return d.Get_PointValue(self,varname,(vlat,vlon))



def nc_pftsum_trans(infile, outfile, length_one_time_dim=False):
    """
    Apply the pftsum transformation for infile and write to outfile. set length_one_time_dim to True to keep the time dimension whose lenght is 1 in the new nc file.
    """
    d = Ncdata(infile)
    d.get_pftsum()
    outfile_ob = NcWrite(outfile)
    if d.unlimited_dimlen == 1:
        if length_one_time_dim:
            outfile_ob.add_3dim_time_lat_lon(1, d.lat, d.lon)
            dim_flag = 2
        else:
            outfile_ob.add_2dim_lat_lon(d.lat, d.lon)
            dim_flag = 1
    else:
        outfile_ob.add_3dim_time_lat_lon(d.unlimited_dimlen, d.lat, d.lon)
        dim_flag = 3

    for varname in d.pftsum.__dict__.keys():
        if varname != 'Areas':
            if dim_flag == 1:
                outfile_ob.add_var_2dim_lat_lon('PFTSUM_'+varname, d.pftsum.__dict__[varname], attr_copy_from=d.d0.__dict__[varname])
            else:
                outfile_ob.add_var_3dim_time_lat_lon('PFTSUM_'+varname, d.pftsum.__dict__[varname], attr_copy_from=d.d0.__dict__[varname])

    if 'Areas' in d.pftsum.__dict__.keys():
        outfile_ob.add_var_2dim_lat_lon('Areas', d.pftsum.__dict__['Areas'], attr_copy_from=d.d0.__dict__['Areas'])

    outfile_ob.close()


def nc_subgrid(infile, outfile, subgrid=[(None,None),(None,None)]):
    """
    Purpose : subgrid infile by vertex. Currently can only handle as much as 4 dimensions (time,PFT,lat,lon)
    Definition : nc_subgrid(infile, outfile, subgrid=[(vlat1,vlat2),(vlon1,vlon2)]):
    Arguments:
        subgrid : [(vlat1,vlat2),(vlon1,vlon2)]; latitude before latitude
    Test : Test has been done for 2dim/3dim/4dim variables. see test_nc_subgrid
    """
    print "subgrid for file --{0}-- begins".format(infile)
    d = Ncdata(infile)
    (lon_index_min, lon_index_max, lat_index_min, lat_index_max) = find_index_by_vertex(d.lonvar, d.latvar, subgrid[1], subgrid[0])
    outfile_ob = NcWrite(outfile)
    outfile_ob.add_dim([d.latdim_name, d.latvar_name, d.latvar_name, 'f4', d.latvar[lat_index_min:lat_index_max+1], 'None', False])
    outfile_ob.add_dim([d.londim_name, d.lonvar_name, d.lonvar_name, 'f4', d.lonvar[lon_index_min:lon_index_max+1], 'None', False])
    if hasattr(d,'unlimited_dimname'):
        outfile_ob.add_dim([d.unlimited_dimname, d.timevar_name, d.timevar_name, 'i4', d.timevar, 'None', True])
    if 'PFT' in d.dimensions:
        outfile_ob.add_dim_pft()

    for varname in pb.StringListAnotB(d.d1.__dict__.keys(), d.dimvar_name_list):
        ordata = d.d0.__dict__[varname]
        assert ordata.ndim>=2,"variable '{0}' has ndim <2".format(varname)
        outfile_ob.add_var(varinfo_value=[varname, ordata.dimensions, ordata.dtype, ordata[:][...,lat_index_min:lat_index_max+1,lon_index_min:lon_index_max+1]], attr_copy_from=ordata)
        print "subgrid for --{0}-- done".format(varname)

    glob_attr_dic = d.global_attributes
    #add time stamp for this operation
    if 'history' in glob_attr_dic:
        glob_attr_dic['history'] = 'file created at ' + str(datetime.datetime.today()) + 'by subgriding file --' + infile + '--\n' + glob_attr_dic['history']
    else:
        glob_attr_dic['history'] = 'file created at ' + str(datetime.datetime.today()) + 'by subgriding file --' + infile

    #write the global attributes
    outfile_ob.add_global_attributes(glob_attr_dic)
    outfile_ob.close()
    print "subgrid for file --{0}-- done".format(outfile)

def nc_subgrid_csv(infile,csv_file=None,bound_name = ['south_bound','north_bound','west_bound','east_bound']):
    """
    Separate a (global) nc file by regions as specified in the csv_file.

    Parameters:
    -----------
    csv_file: the title should be region names, the first column should be the four bound names. see example in base_data/region.csv
    bound_name: the names that are used to distinguish the four bounds: (vlat1,vlat2),(vlon1,vlon2), must be in sequence of south,north,west,east

    Outputs:
    --------
    Several regional nc files, with file names indicating each region.

    Test
    ----
    nc_subgrid_csv and nc_merge_files are tested against each other in the gnc_test.py.
    """
    region=pa.read_csv(csv_file,index_col=0)
    region_dic=region.to_dict()
    for name in region_dic.keys():
        data=region_dic[name]
        (vlat1,vlat2) = (data[bound_name[0]], data[bound_name[1]])
        (vlon1,vlon2) = (data[bound_name[2]], data[bound_name[3]])
        nc_subgrid(infile,infile[:-3]+'_'+name+'.nc',[(vlat1,vlat2),(vlon1,vlon2)])

def compare_nc_file(file1,file2,varlist,npindex=np.s_[:]):
    """
    compare if the input netcdf files have the same value for the same variable.

    Parameters:
    -----------
    npindex: it could be a 2-lenght tuple (which will be broadcast to number of varialbes) or a n-length list with each element as a 2-length tuple
    """
    d1 = Ncdata(file1)
    d2 = Ncdata(file2)
    reslist = []
    if isinstance(npindex,list):
        if len(npindex) != len(varlist):
            raise ValueError("npindex is provided as list but length not equal to varlist!")
        else:
            npindex_final = npindex
    elif isinstance(npindex,tuple):
        if len(npindex) != 2:
            raise ValueError("npindex is provided as tuple and lenght is not 2!")
        else:
            npindex_final = [npindex] * len(varlist)
    elif npindex == np.s_[:]:
        npindex_final = [(npindex,npindex)] * len(varlist)
    else:
        raise ValueError("unknown npindex value!")

    for varname,comindex in zip(varlist,npindex_final):
        if np.ma.allclose(d1.d1.__dict__[varname][comindex[0]],d2.d1.__dict__[varname][comindex[1]]):
            print "variable --{0}-- equal".format(varname)
            reslist.append(True)
        else:
            print "variable --{0}-- NOT equal".format(varname)
            reslist.append(False)

def arithmetic_ncfiles_var(filelist,varlist,func,npindex=np.s_[:]):
    """
    Get arithetic operation form varialbes in nc files by a quick way. This is to avoid the tedious way to try to get some simple arithmetic operation on
        seleted varialbes from some nc files.

    Parameters:
    -----------
    filelist: nc file list.
    varlist: variable name list, when the variable names for all the nc files are the same, it will be broadcast.
    func: function that used the corresponding ndarrays as arguments. np.ma functions are tested to work well.
    npindex: it could be given by a list to specify the indexing of ndarrays before applying the function.
    """

    if isinstance(varlist,str):
        varlist = [varlist]*len(filelist)

    if isinstance(npindex,list):
        if len(npindex) != len(varlist):
            raise ValueError("npindex is provided as list but length not equal to varlist!")
        else:
            npindex_final = npindex
    elif npindex == np.s_[:]:
        npindex_final = [npindex] * len(varlist)
    else:
        raise ValueError("unknown npindex value!")

    data_list = [Ncdata(filename) for filename in filelist]
    ndarray_list = [d.d1.__dict__[varname][comindex] for d,varname,comindex in zip(data_list,varlist,npindex_final)]
    return func(*ndarray_list)

def nc_add_Mdata_by_DictFilenameVarlist(dict_filename_varlist,npindex=np.s_[:]):
    """
    Simple wrapper of Pdata.Mdata

    Notes:
    ------
    1. npindex applies to all filename and varname, thus the usage is
        very limited.

    2. dict_filename_varlist: a dictionary of (filename,varlist) pairs.
    """
    if pb.List_Duplicate_Check(pb.iteflat(dict_filename_varlist.values())):
        raise ValueError("Duplicate varnames found")
        sys.exit()
    md = Pdata.Mdata()
    for filename in dict_filename_varlist.keys():
        d = Ncdata(filename)
        for varname in dict_filename_varlist[filename]:
            tag = varname
            md.add_tag(tag)
            md.add_lat_lon(tag=tag,lat=d.latvar,lon=d.lonvar)
            md.add_array(d.d1.__dict__[varname][npindex],tag)
            if hasattr(d.d0.__dict__[varname],'units'):
                unit = d.d0.__dict__[varname].getncattr('units')
            elif hasattr(d.d0.__dict__[varname],'unit'):
                unit = d.d0.__dict__[varname].getncattr('unit')
            else:
                unit = None
            md.add_attr_by_tag(unit={tag:unit})
    return md

def nc_add_Mdata_mulitfile(filelist,taglist,varlist,npindex=np.s_[:]):
    '''
    Simple wrapper of Pdata.Mdata.add_entry_share_latlon_bydic.

    Parameters:
    -----------
    varlist: varlist length should be equal to filelist. Broadcast of list
        when single string value is provided.
    Notes:
    ------
    All the input nc files should have the same geo_limit
    '''
    if isinstance(varlist,list):
        pass
    elif isinstance(varlist,str):
        varlist = [varlist]*len(filelist)

    ydic = {}
    d0 = Ncdata(filelist[0])
    for filename,tag,varname in zip(filelist,taglist,varlist):
        dt = Ncdata(filename)
        ydic[tag] = dt.d1.__dict__[varname][npindex]
    md = Pdata.Mdata()
    md.add_entry_share_latlon_bydic(ydic, lat=d0.lat, lon=d0.lon)
    return md

def ncfile_get_varlist(filename):
    """
    This is a loose application of gnc.Ncdata.list_var()

    Arguments:
    ----------
    filename: when it's a single string, return varlist of the file;
        when it's a list of file names, return a dictionary of
        (filename,filelist).
    """
    if isinstance(filename,str):
        d = Ncdata(filename)
        return d.list_var(keyword='excludedim')
    elif isinstance(filename,list):
        return dict([(f,ncfile_get_varlist(f)) for f in filename])
        #tlist = [(f,ncfile_get_varlist(f)) for f in filename]
        #return dict(tlist)



def nc_add_Pdata_mulitfile(filelist,taglist,varlist,npindex=np.s_[:]):
    '''
    Simple wrapper of Pdata.Pdata.add_entry_sharex_noerror_by_dic

    Parameters:
    -----------
    varlist: varlist length should be equal to filelist. Broadcast of list
        when single string value is provided.
    Notes:
    ------
    1. all the entry use default indexed x
    '''
    if isinstance(varlist,list):
        pass
    elif isinstance(varlist,str):
        varlist = [varlist]*len(filelist)

    ydic = {}
    d0 = Ncdata(filelist[0])
    for filename,tag,varname in zip(filelist,taglist,varlist):
        dt = Ncdata(filename)
        ydic[tag] = dt.d1.__dict__[varname][npindex]
    pd = Pdata.Pdata()
    pd.add_entry_sharex_noerror_by_dic(ydic)
    return pd

@append_doc_of(_set_default_ncfile_for_write)
def nc_creat_ncfile_by_ncfiles(outfile,varname,input_file_list,input_varlist,
                               npindex=np.s_[:],
                               pyfunc=None,
                               Ncdata_latlon_dim_name=None,
                               attr_kwargs={},
                               **kwargs):
    '''
    Write ncfile by calculating from variables of other nc file.

    Parameters:
    -----------
    varname: the varname for writting into nc file.
    input_file_list: input file list.
    input_varlist: varname list from input files, broadcast when necessary.
    npindex: could be a list otherwise broadcast.
    Ncdata_latlon_dim_name: used in Ncdata function
    attr_kwargs: used in Ncwrite.add_var

    Notes:
    ------
    This is intended to write only ONE variable for the nc file.

    '''
    ncfile = NcWrite(outfile)
    _set_default_ncfile_for_write(ncfile,**kwargs)
    if isinstance(input_varlist,str):
        input_varlist = [input_varlist]*len(input_file_list)
    if not isinstance(npindex,list):
        npindex = [npindex]*len(input_file_list)

    data_list = [Ncdata(filename,latlon_dim_name=Ncdata_latlon_dim_name) for
               filename in input_file_list]
    ndarray_list = [d.d1.__dict__[var][comindex] for d,var,comindex
                    in zip(data_list,input_varlist,npindex)]
    ndim = len(ncfile.dimensions)
    data = pyfunc(ndarray_list)
    ncfile.add_var_smart_ndim(varname,ndim,data,**attr_kwargs)
    ncfile.close()

@append_doc_of(_set_default_ncfile_for_write)
def nc_merge_ncfiles(outfile,input_file_list,
                               Ncdata_latlon_dim_name=None,
                               attr_kwargs={},
                               **kwargs):
    '''
    Merge all the variables in the input_file_list as a single file. The
        duplicate vriables that come later will be discast.

    Parameters:
    -----------
    input_file_list: input file list.
    Ncdata_latlon_dim_name: used in Ncdata function
    attr_kwargs: used in Ncwrite.add_var

    Notes:
    ------
    This is intended to write only ONE variable for the nc file.

    '''
    ncfile = NcWrite(outfile)
    _set_default_ncfile_for_write(ncfile,**kwargs)
    ndim = len(ncfile.dimensions)
    varlist_inside = []

    for filename in input_file_list:
        d = Ncdata(filename,latlon_dim_name=Ncdata_latlon_dim_name)
        varlist = pb.StringListAnotB(d.d1.__dict__.keys(),
                                     d.dimvar_name_list)
        for varname in varlist:
            if varname in varlist_inside:
            #if varname in varlist_inside or varname == 'VEGET_MAX':
            #if varname in varlist_inside or varname == 'Areas' or varname == 'VEGET_MAX':
                print """var --{0}-- in file --{1}-- is discast due to
                    duplication""".format(varname,filename)
            else:
                if varname == 'VEGET_MAX':
                    pftdim = True
                else:
                    pftdim = False
                data = d.d1.__dict__[varname]
                ncfile.add_var_smart_ndim(varname,ndim,data,pftdim=pftdim,
                    attr_copy_from=d.d0.__dict__[varname],**attr_kwargs)
                varlist_inside.append(varname)
    ncfile.close()


def test_nc_subgrid():
    #test 3dim nc data
    nc_subgrid('testdata/grid_test_3dim.nc','testdata/grid_test_3dim_sub.nc',subgrid=[(-35,25),(-50,20)])
    d0 = Ncdata('testdata/grid_test_3dim.nc')
    d1 = Ncdata('testdata/grid_test_3dim_sub.nc')
    lon = np.arange(-49.75,19.8,0.5)
    lat = np.arange(-34.75,24.8,0.5)
    #test for 2dim variable
    assert np.array_equal(d1.d1.xx2dim,np.tile(lon,(len(lat),1)))
    assert np.array_equal(d1.d1.yy2dim,np.tile(lat[::-1][:,np.newaxis],(1,len(lon))))
    #test for 3dim variable
    assert np.array_equal(d1.d1.xx3dim,np.tile(lon,(d0.unlimited_dimlen, len(lat), 1)))
    assert np.array_equal(d1.d1.yy3dim,np.tile(lat[::-1][:,np.newaxis],(d0.unlimited_dimlen, 1, len(lon))))

    nc_subgrid('testdata/grid_test_4dim.nc','testdata/grid_test_4dim_sub.nc',subgrid=[(-35,25),(-50,20)])
    d0 = Ncdata('testdata/grid_test_4dim.nc')
    d1 = Ncdata('testdata/grid_test_4dim_sub.nc')
    lon = np.arange(-49.75,19.8,0.5)
    lat = np.arange(-34.75,24.8,0.5)
    #test for 4dim variable
    assert np.array_equal(d1.d1.xx4dim,np.tile(lon,(d0.unlimited_dimlen, 13, len(lat), 1)))
    assert np.array_equal(d1.d1.yy4dim,np.tile(lat[::-1][:,np.newaxis],(d0.unlimited_dimlen, 13, 1, len(lon))))

def test_Ncdata():
    def test_3dim():
        #0 prepare test data
        dt = Ncdata('testdata/stomate_history_AS_1240_TOTAL_M_LITTER_CONSUMP_3dim.nc')
        #0.1 test PFT sum operation
        nonzero_veget_max = np.ma.masked_equal(dt.d1.VEGET_MAX, 0)
        pft_total_m = pb.MaskArrayByNan(dt.d1.TOTAL_M) * nonzero_veget_max
        pftsum_total_m = np.ma.sum(pft_total_m, axis=0)
        #0.2 test spatial operation
        pftsum_total_m_area = pftsum_total_m * dt.d1.Areas
        spasum_total_m = np.ma.sum(np.ma.sum(pftsum_total_m_area,axis=-1), axis=-1)
        spamean_total_m = spasum/np.ma.sum(dt.d1.Areas)

        #1 do the test
        dnew = Ncdata('testdata/stomate_history_AS_1240_TOTAL_M_LITTER_CONSUMP_3dim.nc')
        dnew.get_spa()
        assert np.array_equal(dnew.pftsum.TOTAL_M, pftsum_total_m)
        assert np.array_equal(dnew.pftsum.LITTER_CONSUMP, dt.d1.LITTER_CONSUMP)
        pdb.set_trace()
        assert dnew.spasum.TOTAL_M == spasum_total_m
        assert dnew.spamean.TOTAL_M == spamean_total_m

    def test_4dim():
        #0 prepare test data
        dt = Ncdata('testdata/stomate_history_TOTAL_M_LITTER_CONSUMP_4dim.nc')
        #0.1 test PFT sum operation
        nonzero_veget_max = np.ma.masked_equal(dt.d1.VEGET_MAX, 0)
        pft_total_m = pb.MaskArrayByNan(dt.d1.TOTAL_M) * nonzero_veget_max
        pftsum_total_m = np.ma.sum(pft_total_m, axis=1)
        #0.2 test spatial operation
        pftsum_total_m_area = pftsum_total_m * dt.d1.Areas
        spasum_total_m = np.ma.sum(np.ma.sum(pftsum_total_m_area,axis=-1), axis=-1)
        spamean_total_m = spasum/np.ma.sum(dt.d1.Areas)

        #1 do the test
        dnew = Ncdata('testdata/stomate_history_TOTAL_M_LITTER_CONSUMP_4dim.nc')
        dnew.get_spa()
        assert np.array_equal(dnew.pftsum.TOTAL_M, pftsum_total_m)
        assert np.array_equal(dnew.pftsum.LITTER_CONSUMP, dt.d1.LITTER_CONSUMP)
        pdb.set_trace()
        assert dnew.spasum.TOTAL_M == spasum_total_m
        assert dnew.spamean.TOTAL_M == spamean_total_m





