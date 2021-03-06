import bsite as bsite
import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import pickle as pk
import mathex as mathex
import os as os
import re as re
import scipy as sp
import mpl_toolkits.basemap as bmp
from mpl_toolkits.basemap import cm
import pdb
import netCDF4 as nc
from matplotlib.backends.backend_pdf import PdfPages
import copy as pcopy
import g
import pb
import tools

def near5even(datain):
    if datain%5==0:
        dataout=datain
    else:
        if datain/5<0:
            dataout=np.ceil(datain/5)*5
        else:
            dataout=np.floor(datain/5)*5
    return dataout

class gmap(object):
    """
    Purpose: plot the map used for later contour or image plot.
    Note:
        return m,lonpro,latpro,latind,lonind
        1. return m --> map drawed; lonpro/latpro --> lat/lon transferred
            to projection coords; latind/lonind --> index for lat/lon
            falling with mapbound
        2. lat must be descending and lon must be ascending.

    Parameters:
    -----------
    kwargs: used for basemap.Basemap method.

    Example:
        >>> fig,ax=g.Create_1Axes()
        >>> m,lonpro,latpro,lonind,latind=bmap.gmap(ax,'cyl',mapbound='all',lat=np.arange(89.75,-89.8,-0.5),lon=np.arange(-179.75,179.8,0.5),gridstep=(30,30))
        >>> x,y=m(116,40)  #plot Beijing
        >>> m.scatter(x,y,s=30,marker='o',color='r')
    """
    def __init__(self,ax=None,projection='cyl',mapbound='all',lat=None,lon=None,
                 gridstep=(30,30),**kwargs):

        ax = tools._replace_none_axes(ax)
        lat = tools._replace_none_by_given(lat,np.arange(89.75,-89.8,-0.5))
        lon = tools._replace_none_by_given(lon,np.arange(-179.75,179.8,0.5))

        latstep = lat[0] - lat[1]
        if latstep <= 0:
            raise TypeError("lat input is increasing!")
        else:
            if latstep == 0.5:
                half_degree = True
            else:
                half_degree = False

        if projection=='cyl':
            if isinstance(mapbound,dict):
                raise ValueError('cannot use dict for cyl projection')
            elif mapbound=='all':
                lat1=lat[-1]
                lat2=lat[0]
                lon1=lon[0]
                lon2=lon[-1]
                #when the data is of half degree resolution, often the lat1 and
                #lat2 is in the center of the half degree cell, so we need to
                #adjust for the vertices.
                if half_degree == True:
                    if lat1%0.25 == 0:
                        lat1 = lat1-0.25
                    if lat2%0.25 == 0:
                        lat2 = lat2+0.25
                    if lon1%0.25 == 0:
                        lon1 = lon1-0.25
                    if lon2%0.25 == 0:
                        lon2 = lon2+0.25
                if lat1<-85:
                    lat1=-90.
                if lat2>85:
                    lat2=90.
                if lon1<-175:
                    lon1=-180.
                if 185>lon2>175:
                    lon2=180.
                if lon2>355:
                    lon2=360.
            else:
                lat1=mapbound[0]
                lat2=mapbound[1]
                lon1=mapbound[2]
                lon2=mapbound[3]
            #draw the map, parallels and meridians
            m=bmp.Basemap(projection=projection,llcrnrlat=lat1,urcrnrlat=lat2,
                          llcrnrlon=lon1,urcrnrlon=lon2,resolution='l',ax=ax,
                          **kwargs)
            m.drawcoastlines(linewidth=0.7)
            if gridstep!=None and gridstep!=False:
                para_range=np.arange(near5even(lat1),near5even(lat2)+0.1,gridstep[0])
                meri_range=np.arange(near5even(lon1),near5even(lon2)+0.1,gridstep[1])
                m.drawparallels(para_range,labels=[1,0,0,0])
                m.drawmeridians(meri_range,labels=[0,0,0,1])
            #make the grid
            latind=np.nonzero((lat>lat1)&(lat<lat2))[0]
            lonind=np.nonzero((lon>lon1)&(lon<lon2))[0]
            numlat=len(np.nonzero((lat>lat1)&(lat<lat2))[0])
            numlon=len(np.nonzero((lon>lon1)&(lon<lon2))[0])
            lonm,latm=m.makegrid(numlon,numlat)
            latm=np.flipud(latm)
            lonpro,latpro=m(lonm,latm)

        elif projection=='npstere':
            if not isinstance(mapbound,dict):
                raise ValueError('please use dict to specify')
            else:
                m=bmp.Basemap(projection='npstere',boundinglat=mapbound['blat'],
                              lon_0=mapbound['lon_0'],resolution='l',ax=ax,
                              **kwargs)
                m.drawcoastlines(linewidth=0.7)
                m.fillcontinents(color='0.8',zorder=0)
                if gridstep!=None and gridstep!=False:
                    m.drawparallels(np.arange(mapbound['para0'],91.,gridstep[0]),
                                    labels=[1,0,0,0],fontsize=10)
                    m.drawmeridians(np.arange(-180.,181.,gridstep[1]),
                                    labels=[0,0,0,0],fontsize=10)
                #make the grid
                lat1=mapbound['blat']
                latind=np.nonzero(lat>lat1)[0]
                lonind=np.arange(len(lon))
                latnew=np.linspace(90, lat1, num=len(latind), endpoint=True)
                if lon[-1]>180:
                    lonnew=np.linspace(0,360,num=len(lonind),endpoint=True)
                else:
                    lonnew=np.linspace(-180,180,num=len(lonind),endpoint=True)
                lonm,latm=np.meshgrid(lonnew,latnew)
                lonpro,latpro=m(lonm,latm)
        else:
            raise ValueError('''projection '{0}' not supported'''
                             .format(projection))
        self.m = m
        self.lonpro = lonpro
        self.latpro = latpro
        self.latind = latind
        self.lonind = lonind



def _transform_data(pdata,levels,data_transform):
    '''
    Return [pdata,plotlev,plotlab,extend,trans_base_list];
    if data_transform == False, trans_base_list = None.

    Notes:
    ------
    pdata: data used for contourf plotting.
    plotlev: the levels used in contourf plotting.
    extend: the value for parameter extand in contourf.
    trans_base_list: cf. mathex.plot_array_transg
    '''
    if levels==None:
        ftuple = (pdata,None,None,'neither')
        if data_transform==True:
            raise Warning("Strange levels is None but data_transform is True")
    else:
        if data_transform==True:
            #make the data transform before plotting.
            pdata_trans,plotlev,plotlab,trans_base_list = \
                mathex.plot_array_transg(pdata, levels, copy=True)
            if np.isneginf(plotlab[0]) and np.isposinf(plotlab[-1]):
                ftuple = (pdata_trans,plotlev[1:-1],plotlab,'both')
            elif np.isneginf(plotlab[0]) or np.isposinf(plotlab[-1]):
                raise ValueError('''only one extreme set as infinitive, please
                    set both as infinitive if arrow colorbar is wanted.''')
            else:
                ftuple = (pdata_trans,plotlev,plotlab,'neither')
        #data_transform==False
        else:
            plotlev=pb.iteflat(levels)
            plotlab=pb.iteflat(levels)
            if np.isneginf(plotlab[0]) and np.isposinf(plotlab[-1]):
                #here the levels would be like [np.NINF,1,2,3,np.PINF]
                #in following contourf, all values <1 and all values>3 will be
                #automatically plotted in the color of two arrows.
                #easy to see in this example:
                #a=np.tile(np.arange(10),10).reshape(10,10);
                #fig,ax=g.Create_1Axes();
                #cs=ax.contourf(a,levels=np.arange(2,7),extend='both');
                #plt.colorbar(cs)
                ftuple = (pdata,plotlev[1:-1],plotlab,'both')
            elif np.isneginf(plotlab[0]) or np.isposinf(plotlab[-1]):
                raise ValueError('''only one extreme set as infinitive, please
                    set both as infinitive if arrow colorbar is wanted.''')
            else:
                ftuple = (pdata,plotlev,plotlab,'neither')
    datalist = list(ftuple)

    if data_transform == True:
        datalist.append(trans_base_list)
    else:
        datalist.append(None)
    return datalist

def _generate_colorbar_ticks_label(data_transform=False,
                                   colorbarlabel=None,
                                   trans_base_list=None,
                                   forcelabel=None,
                                   plotlev=None,
                                   plotlab=None):
    '''
    Return (colorbar_ticks,colorbar_labels)
    '''
    #data_transform==True and levels!=None
    if data_transform==True:
        if colorbarlabel != None:
            colorbarlabel=pb.iteflat(colorbarlabel)
            transformed_colorbarlabel_ticks,x,y,trans_base_list = \
                mathex.plot_array_transg(colorbarlabel, trans_base_list,
                                         copy=True)

        #Note if/else blocks are organized in 1st tire by check if the two
        #ends are -inf/inf and 2nd tire by check if colorbarlabel is None
        if np.isneginf(plotlab[0]) and np.isposinf(plotlab[-1]):
            if colorbarlabel!=None:
                ftuple = (transformed_colorbarlabel_ticks,colorbarlabel)
            else:
                ftuple = (plotlev,plotlab[1:-1])
        elif np.isneginf(plotlab[0]) or np.isposinf(plotlab[-1]):
            raise ValueError("It's strange to set only side as infitive")
        else:
            if colorbarlabel!=None:
                ftuple = (transformed_colorbarlabel_ticks,colorbarlabel)
            else:
                ftuple = (plotlev,plotlab)

    #data_transform==False
    else:
        if np.isneginf(plotlab[0]) and np.isposinf(plotlab[-1]):
            #if colorbarlabel is forced, then ticks and ticklabels will be forced.
            if colorbarlabel!=None:
                ftuple = (colorbarlabel,colorbarlabel)
            #This by default will be done, it's maintained here only for clarity.
            else:
                ftuple = (plotlab[1:-1],plotlab[1:-1])
        elif np.isneginf(plotlab[0]) or np.isposinf(plotlab[-1]):
            raise ValueError("It's strange to set only side as infitive")
        else:
            if colorbarlabel!=None:
                ftuple = (colorbarlabel,colorbarlabel)
            else:
                ftuple = (plotlab,plotlab)

    ftuple = list(ftuple)
    if forcelabel != None:
        if len(forcelabel) != len(ftuple[1]):
            raise ValueError('''the length of the forcelabel and the
                length of labeled ticks is not equal!''')
        else:
            ftuple[1] = forcelabel

    return ftuple

def _generate_smartlevel(pdata):
    """
    generate smart levels by using the min, percentiles from 5th
        to 95th with every 5 as the step, and the max value.
    """

    def even_num(num):
        if num >= 10:
            return int(num)
        else:
            return round(num,4)

    def extract_percentile(array,per):
        return even_num(np.percentile(array,per))

    def generate_smartlevel_from_1Darray(array):
        vmax = even_num(np.max(array))
        vmin = even_num(np.min(array))
        per_level = map(lambda x:extract_percentile(array,x),
                        np.arange(5,96,5))
        return np.array([vmin]+per_level+[vmax])

    if np.isnan(np.sum(pdata)):
        pdata = np.ma.masked_invalid(pdata)

    if np.ma.isMA(pdata):
        array1D = pdata[np.nonzero(~pdata.mask)]
    else:
        array1D = pdata.flatten()

    return generate_smartlevel_from_1Darray(array1D)

def _generate_map_prepare_data(data=None,lat=None,lon=None,
                               projection='cyl',
                               mapbound='all',
                               gridstep=(30,30),
                               shift=False,
                               map_threshold=None,
                               levels=None,
                               cmap=None,
                               smartlevel=None,
                               data_transform=False,
                               gmapkw={},
                               ax=None):
    """
    This function makes the map, and transform data for ready
        use of m.contourf or m.imshow
    """
    if shift==True:
        data,lon=bmp.shiftgrid(180,data,lon,start=False)
    mgmap=gmap(ax,projection,mapbound,lat,lon,gridstep,**gmapkw)
    m,lonpro,latpro,latind,lonind = (mgmap.m, mgmap.lonpro, mgmap.latpro,
                                    mgmap.latind, mgmap.lonind)

    pdata = data[latind[0]:latind[-1]+1,lonind[0]:lonind[-1]+1]
    #mask by map_threshold
    pdata = mathex.ndarray_mask_by_threshold(pdata,map_threshold)
    #generate the smartlevel
    if smartlevel == True:
        if levels != None:
            raise ValueError("levels must be None when smartlevel is True!")
        else:
            levels = _generate_smartlevel(pdata)
            data_transform = True
    #prepare the data for contourf
    pdata,plotlev,plotlab,extend,trans_base_list = \
        _transform_data(pdata,levels,data_transform)
    return (mgmap,pdata,plotlev,plotlab,extend,
            trans_base_list,data_transform)

def _set_colorbar(m,cs,colorbardic={},
                  levels=None,
                  data_transform=False,
                  colorbarlabel=None,
                  trans_base_list=None,
                  forcelabel=None,
                  show_colorbar=True,
                  plotlev=None,
                  plotlab=None,
                  cbarkw={}):
    """
    Wrap the process for setting colorbar.
    """
    #handle the colorbar attributes by using dictionary which flexibility.
    if show_colorbar == False:
        cbar = None
    else:
        location = colorbardic.get('location','right')
        size = colorbardic.get('size','3%')
        pad = colorbardic.get('pad','2%')
        cbar=m.colorbar(cs,location=location, size=size, pad=pad,**cbarkw)
        #set colorbar ticks and colorbar label
        if levels==None:
            pass
        else:
            ticks,labels = \
                _generate_colorbar_ticks_label(data_transform=data_transform,
                                               colorbarlabel=colorbarlabel,
                                               trans_base_list=trans_base_list,
                                               forcelabel=forcelabel,
                                               plotlev=plotlev,
                                               plotlab=plotlab)
            cbar.set_ticks(ticks)
            cbar.set_ticklabels(labels)
    return cbar

class mapcontourf(object):
    """
    Purpose: plot a map on 'cyl' or 'npstere' projection.
    Arguments:
        ax --> An axes instance
        projection --> for now two projections have been added:
            1. 'cyl' -- for global and regional mapping
            2. 'npstere' -- for north polar centered mapping.
        lat,lon --> geographic coordinate variables; lat must be in
            desceding order and lon must be ascending.
        mapbound --> specify the bound for mapping;
            1. 'cyl'
                tuple containing (lat1,lat2,lon1,lon2); lat1 --> lower
                parallel; lat2 --> upper parallel; lon1 --> left meridian;
                lon2 --> right meridian; default 'all' means plot
                the extent of input lat, lon coordinate variables;
                for global mapping, set (-90,90,-180,180) or (-90,90,0,360).
            2. 'npstere'
                mapbound={'blat':45,'lon_0':0,'para0':40}
                blat --> boundinglat in the bmp.Basemap method.
                    The souther limit for mapping.
                lon_0 --> center of desired map domain.
                para0 --> souther boundary for parallel ticks, the default
                    norther limit is 90; default longitude 0-360 (or -180-180)
        gridstep --> the step for parallel and meridian grid for the map,
            tuple containing (parallel_step, meridian_step).
        levels --> default None; levels=[-5,-2,-1,0,1,2,5] ;
            or levels=[(-10,-4,-2,-1,-0.4),(-0.2,-0.1,0,0.1,0.2),
                        (0.4,1,2,4,10)].
            1.  Anything that can work as input for function pb.iteflat()
                will work.
            2.  If the first and last element of pb.iteflat(levels) is
                np.NINF and np.PINF, the colorbar of contourf plot will
                use the 'two-arrow' shape.
            3.  If data_transform==True, the input data will be transformed
                from pb.iteflat(levels) to
                np.linspace(1,len(pb.iteflat(interval_original)). this can
                help to create arbitrary contrasting in the plot.
                cf. mathex.plot_array_transg
        smartlevel:
            1. when True, a "smart" level will be generated by
               using the min,max value and the [5th, 10th, ..., 95th]
               percentile of the input array.
            2. it will be applied after applying the mask_threshold.
        data_transform:
            1. set as True if increased contrast in the plot is desired.
                In this case the function mathex.plot_array_transg will
                be called and pb.iteflat(levels) will be used as original
                interval for data transformation.
            2. In case of data_transform==False, pb.iteflat(levels)
                will be used directly in the plt.contour function for
                ploting and hence no data transformation is made. The
                treatment by this way allows very flexible
                (in a mixed way) to set levels.
            3. In any case, if np.NINF and np.PINF as used as two
                extremes of levels, arrowed colorbar will be returned.

        colorbarlabel:
            1. used to put customized colorbar label and this will override
                using levels as colorbar. IF colorbarlabel!=None,
                colorbar ticks and labels will be set using colorbarlabel.
                so this means colorbarlabel could only be array or
                list of numbers.
            2. If data_transform==True, colorbar will also be transformed
                accordingly. In this case, the colorbar ticks will use
                transformed colorbarlabel data, but colorbar ticklables
                will use non-transformed colorbarlabel data. This means
                the actual ticks numbers and labels are not the same.

        forcelabel --> to force the colorbar label as specified by forcelabel.
            This is used in case to set the labels not in numbers but in
            other forms (eg. strings).
            In case of data_transform = True, levels will be used to
            specifiy levels for the original colorbar, colorbarlabel will
            be used to create ticks on colrobar which will be labeled,
            if forcelabel=None, then colorbarlabel will agined be used
            to label the ticks, otherwise forcelabel will be used to
            label the ticks on the colorbar. So this means forcelabel will
            mainly be list of strings.

        data --> numpy array with dimension of len(lat)Xlen(lon)
        map_threshold --> dictionary like {'lb':2000,'ub':5000}, data
            less than 2000 and greater than 5000 will be masked.
            Note this will be applied before data.
               transform.
        shift --> boolean value. False for longtitude data ranging [-180,180];
            for longtitude data ranging [0,360] set shift to True if a
            180 east shift is desired. if shift as True, the mapbound
            range should be set using shifted longtitude
            (use -180,180 rather than 0,360).
        colorbardic --> dictionary to specify the attributes for colorbar,
            translate all the keys in function bmp.Basemap.colorbar()
            into keys in colorbardic to manipulation.

    Note:
        1.  lat must be descending, and lon must be ascending.
        2*. NOTE use both data_transform=True and impose unequal
            colorbarlabel could be very confusing! Because normaly in
            case of data_transform as True the labels are ALREADY
            UNEQUALLY distributed!

            an example to use colorbarlabel and forcelabel:
            data_transform=True,
            levels=[0,1,2,3,4,5,6,7,8]
            colorbarlabel=[0,2,4,6,8]
            forcelabel=['extreme low','low','middle','high','extreme high']

            So colorbarlabel will set both ticks and labels, but forcelabel
            will further overwrite the labels.

        3. This function has been test using data, the script and
            generated PNG files are availabe at ~/python/bmaptest
    See also:
        mathex.plot_array_transg
    """
    def __init__(self,data=None,lat=None,lon=None,ax=None,
                 projection='cyl',mapbound='all',
                 gridstep=(30,30),shift=False,map_threshold=None,
                 cmap=None,colorbarlabel=None,forcelabel=None,
                 show_colorbar=True,
                 smartlevel=False,
                 levels=None,data_transform=False,
                 colorbardic={},
                 cbarkw={},
                 gmapkw={}
                 ):


        (mgmap,pdata,plotlev,plotlab,extend,
         trans_base_list,data_transform) = \
            _generate_map_prepare_data(data=data,lat=lat,lon=lon,
                                       projection=projection,
                                       mapbound=mapbound,
                                       gridstep=gridstep,
                                       shift=shift,
                                       map_threshold=map_threshold,
                                       levels=levels,
                                       cmap=cmap,
                                       smartlevel=smartlevel,
                                       data_transform=data_transform,
                                       gmapkw=gmapkw,
                                       ax=ax)

        #make the contourf plot
        cs=mgmap.m.contourf(mgmap.lonpro,mgmap.latpro,pdata,
                            levels=plotlev,extend=extend,cmap=cmap)
        ##handle colorbar
        cbar = _set_colorbar(mgmap.m,cs,
                             colorbardic=colorbardic,
                             levels=plotlev,
                             data_transform=data_transform,
                             colorbarlabel=colorbarlabel,
                             trans_base_list=trans_base_list,
                             forcelabel=forcelabel,
                             plotlev=plotlev,
                             plotlab=plotlab,
                             cbarkw=cbarkw,
                             show_colorbar=show_colorbar)
        #return
        self.m = mgmap.m
        self.cs = cs
        self.cbar = cbar
        self.plotlev = plotlev
        self.plotlab = plotlab
        self.ax = mgmap.m.ax
        self.trans_base_list = trans_base_list
        self.gmap = mgmap

        if levels == None:
            pass
        else:
            cbar_ticks,cbar_labels = \
                _generate_colorbar_ticks_label(data_transform=data_transform,
                                               colorbarlabel=colorbarlabel,
                                               trans_base_list=trans_base_list,
                                               forcelabel=forcelabel,
                                               plotlev=plotlev,
                                               plotlab=plotlab)

            self.cbar_ticks = cbar_ticks
            self.cbar_labels = cbar_labels

    def colorbar(self,cax=None,**kwargs):
        """
        set colorbar on specified cax.

        kwargs applies for plt.colorbar
        """
        cbar = plt.colorbar(self.cs,cax=cax,**kwargs)
        cbar.set_ticks(self.cbar_ticks)
        cbar.set_ticklabels(self.cbar_labels)
        return cbar



class mapimshow(object):
    """
    Purpose: plot a map on cyl projection.
    Arguments:
        ax --> An axes instance
        lat,lon --> geographic coordinate variables;
        mapbound --> tuple containing (lat1,lat2,lon1,lon2);
                lat1 --> lower parallel; lat2 --> upper parallel;
                lon1 --> left meridian; lon2 --> right meridian;
                default 'all' means plot the extent of input lat, lon
                coordinate variables;
        gridstep --> the step for parallel and meridian grid for the map,
            tuple containing (parallel_step, meridian_step).
        vmin,vmax --> as in plt.imshow function
        data --> numpy array with dimension of len(lat)Xlen(lon)
        shift --> boolean value. False for longtitude data ranging [-180,180];
            for longtitude data ranging [0,360] set shift to True if
            a 180 east shift is desired.
    """
    def __init__(self,data=None,lat=None,lon=None,ax=None,
                 projection='cyl',mapbound='all',
                 gridstep=(30,30),shift=False,map_threshold=None,
                 cmap=None,colorbarlabel=None,forcelabel=None,
                 smartlevel=False,
                 levels=None,data_transform=False,
                 colorbardic={},
                 cbarkw={},
                 gmapkw={},
                 *args,
                 **kwargs):

        (mgmap,pdata,plotlev,plotlab,extend,
         trans_base_list,data_transform) = \
            _generate_map_prepare_data(data=data,lat=lat,lon=lon,
                                       projection=projection,
                                       mapbound=mapbound,
                                       gridstep=gridstep,
                                       shift=shift,
                                       map_threshold=map_threshold,
                                       levels=levels,
                                       cmap=cmap,
                                       smartlevel=smartlevel,
                                       data_transform=data_transform,
                                       gmapkw=gmapkw,
                                       ax=ax)

        cs=mgmap.m.imshow(pdata,origin='upper',*args,**kwargs)

        cbar = _set_colorbar(mgmap.m,cs,
                             colorbardic=colorbardic,
                             levels=plotlev,
                             data_transform=data_transform,
                             colorbarlabel=colorbarlabel,
                             trans_base_list=trans_base_list,
                             forcelabel=forcelabel,
                             plotlev=plotlev,
                             plotlab=plotlab,
                             cbarkw=cbarkw)

        self.m = mgmap.m
        self.cs = cs
        self.cbar = cbar
        self.plotlev = plotlev
        self.plotlab = plotlab
        self.ax = mgmap.m.ax
        self.trans_base_list = trans_base_list
        self.gmap = mgmap




































































































































#############################################################################
##################      Below deprecated functions      #####################

def makemap(ax,projection='cyl',mapbound='all',lat=None,lon=None,
            gridstep=(30,30),half_degree=True,**kwargs):
    """
    Purpose: plot the map used for later contour or image plot.
    Note:
        return m,lonpro,latpro,latind,lonind
        1. return m --> map drawed; lonpro/latpro --> lat/lon transferred
            to projection coords; latind/lonind --> index for lat/lon
            falling with mapbound
        2. lat must be descending and lon must be ascending.
    Example:
        >>> fig,ax=g.Create_1Axes()
        >>> m,lonpro,latpro,lonind,latind=bmap.makemap(ax,'cyl',mapbound='all',lat=np.arange(89.75,-89.8,-0.5),lon=np.arange(-179.75,179.8,0.5),gridstep=(30,30))
        >>> x,y=m(116,40)  #plot Beijing
        >>> m.scatter(x,y,s=30,marker='o',color='r')
    """
    print '''!!Deprecate Warning: please use gmap instead'''
    if projection=='cyl':
        if isinstance(mapbound,dict):
            raise ValueError('cannot use dict for cyl projection')
        elif mapbound=='all':
            lat1=lat[-1]
            lat2=lat[0]
            lon1=lon[0]
            lon2=lon[-1]
            #when the data is of half degree resolution, often the lat1 and
            #lat2 is in the center of the half degree cell, so we need to
            #adjust for the vertices.
            if half_degree == True:
                if lat1%0.25 == 0:
                    lat1 = lat1-0.25
                if lat2%0.25 == 0:
                    lat2 = lat2+0.25
            if lat1<-85:
                lat1=-90.
            if lat2>85:
                lat2=90.
            if lon1<-175:
                lon1=-180.
            if 185>lon2>175:
                lon2=180.
            if lon2>355:
                lon2=360.
        else:
            lat1=mapbound[0]
            lat2=mapbound[1]
            lon1=mapbound[2]
            lon2=mapbound[3]
        #draw the map, parallels and meridians
        m=bmp.Basemap(projection=projection,llcrnrlat=lat1,urcrnrlat=lat2,
                      llcrnrlon=lon1,urcrnrlon=lon2,resolution='l',ax=ax,
                      **kwargs)
        m.drawcoastlines(linewidth=0.7)
        if gridstep!=None:
            para_range=np.arange(near5even(lat1),near5even(lat2)+0.1,gridstep[0])
            meri_range=np.arange(near5even(lon1),near5even(lon2)+0.1,gridstep[1])
            m.drawparallels(para_range,labels=[1,0,0,0])
            m.drawmeridians(meri_range,labels=[0,0,0,1])
        #make the grid
        latind=np.nonzero((lat>lat1)&(lat<lat2))[0]
        lonind=np.nonzero((lon>lon1)&(lon<lon2))[0]
        numlat=len(np.nonzero((lat>lat1)&(lat<lat2))[0])
        numlon=len(np.nonzero((lon>lon1)&(lon<lon2))[0])
        lonm,latm=m.makegrid(numlon,numlat)
        latm=np.flipud(latm)
        lonpro,latpro=m(lonm,latm)
        return m,lonpro,latpro,latind,lonind

    elif projection=='npstere':
        if not isinstance(mapbound,dict):
            raise ValueError('please use dict to specify')
        else:
            m=bmp.Basemap(projection='npstere',boundinglat=mapbound['blat'],
                          lon_0=mapbound['lon_0'],resolution='l',ax=ax,
                          **kwargs)
            m.drawcoastlines(linewidth=0.7)
            m.fillcontinents(color='0.8',zorder=0)
            if gridstep not in [None,False]:
                m.drawparallels(np.arange(mapbound['para0'],91.,gridstep[0]),
                                labels=[1,0,0,0],fontsize=10)
                m.drawmeridians(np.arange(-180.,181.,gridstep[1]),
                                labels=[0,0,0,0],fontsize=10)
            #make the grid
            lat1=mapbound['blat']
            latind=np.nonzero(lat>lat1)[0]
            lonind=np.arange(len(lon))
            latnew=np.linspace(90, lat1, num=len(latind), endpoint=True)
            if lon[-1]>180:
                lonnew=np.linspace(0,360,num=len(lonind),endpoint=True)
            else:
                lonnew=np.linspace(-180,180,num=len(lonind),endpoint=True)
            lonm,latm=np.meshgrid(lonnew,latnew)
            lonpro,latpro=m(lonm,latm)
            return m,lonpro,latpro,latind,lonind


def imshowmap(lat,lon,indata,ax=None,projection='cyl',mapbound='all',
              gridstep=(30,30),shift=False,colorbar=True,
              colorbarlabel=None,*args,**kwargs):
    """
    Purpose: plot a map on cyl projection.
    Arguments:
        ax --> An axes instance
        lat,lon --> geographic coordinate variables;
        mapbound --> tuple containing (lat1,lat2,lon1,lon2);
                lat1 --> lower parallel; lat2 --> upper parallel;
                lon1 --> left meridian; lon2 --> right meridian;
                default 'all' means plot the extent of input lat, lon
                coordinate variables;
        gridstep --> the step for parallel and meridian grid for the map,
            tuple containing (parallel_step, meridian_step).
        vmin,vmax --> as in plt.imshow function
        indata --> numpy array with dimension of len(lat)Xlen(lon)
        shift --> boolean value. False for longtitude data ranging [-180,180];
            for longtitude data ranging [0,360] set shift to True if
            a 180 east shift is desired.
    """
    print "Deprecate Warning! imshowmap replaced by mapimshow"
    #handle the case ax==None:
    if ax==None:
        fig,axt=g.Create_1Axes()
    else:
        axt=ax

    if shift==True:
        indata,lon=bmp.shiftgrid(180,indata,lon,start=False)

    #make the map and use mapbound to cut the data
    m,lonpro,latpro,latind,lonind=makemap(axt, projection, mapbound,
                                          lat, lon, gridstep)
    pdata=indata[latind[0]:latind[-1]+1,lonind[0]:lonind[-1]+1]

    cs=m.imshow(pdata,origin='upper',*args,**kwargs)
    if colorbar==True:
        cbar=m.colorbar(cs)
        if colorbarlabel!=None:
            cbar.set_label(colorbarlabel)
    else:
        cbar=None
    return m,cs,cbar



def contourfmap2(lat,lon,indata,projection='cyl',mapbound='all',
                 gridstep=(30,30),shift=False,cmap=None,
                 map_threshold=None,colorbarlabel=None,
                 levels=None,data_transform=False,
                 ax=None,colorbardic={}):
    """
    contourfmap2 is a wrapper of contourfmap. NO need to set up
        a figure and axes before drawing map, return fig,ax,m,cbar.
    """
    print "Deprecate Warning 'contourfmap2', use mapcontourf instead."
    if ax==None:
        fig,axt=g.Create_1Axes()
    else:
        axt=ax
    m,bar=contourfmap(axt,lat,lon,indata,projection=projection,
                      mapbound=mapbound,gridstep=gridstep,
                      shift=shift,cmap=cmap,
                      map_threshold=map_threshold,
                      colorbarlabel=colorbarlabel,
                      levels=levels,
                      data_transform=data_transform,
                      colorbardic=colorbardic)
    if ax==None:
        return fig,axt,m,bar
    else:
        return m,bar


def contourfmap(ax=None,lat=None,lon=None,indata=None,projection='cyl',mapbound='all',
                gridstep=(30,30),shift=False,map_threshold=None,
                cmap=None,colorbarlabel=None,forcelabel=None,
                levels=None,data_transform=False,
                return_lev_lab=False,colorbardic={},
                cbarkw={}):
    """
    Purpose: plot a map on 'cyl' or 'npstere' projection.
    Arguments:
        ax --> An axes instance
        projection --> for now two projections have been added:
            1. 'cyl' -- for global and regional mapping
            2. 'npstere' -- for north polar centered mapping.
        lat,lon --> geographic coordinate variables; lat must be in
            desceding order and lon must be ascending.
        mapbound --> specify the bound for mapping;
            1. 'cyl'
                tuple containing (lat1,lat2,lon1,lon2); lat1 --> lower
                parallel; lat2 --> upper parallel; lon1 --> left meridian;
                lon2 --> right meridian; default 'all' means plot
                the extent of input lat, lon coordinate variables;
                for global mapping, set (-90,90,-180,180) or (-90,90,0,360).
            2. 'npstere'
                mapbound={'blat':45,'lon_0':0,'para0':40}
                blat --> boundinglat in the bmp.Basemap method.
                    The souther limit for mapping.
                lon_0 --> center of desired map domain.
                para0 --> souther boundary for parallel ticks, the default
                    norther limit is 90; default longitude 0-360 (or -180-180)
        gridstep --> the step for parallel and meridian grid for the map,
            tuple containing (parallel_step, meridian_step).
        levels --> default None; levels=[-5,-2,-1,0,1,2,5] ;
            or levels=[(-10,-4,-2,-1,-0.4),(-0.2,-0.1,0,0.1,0.2),
                        (0.4,1,2,4,10)].
            1.  Anything that can work as input for function pb.iteflat()
                will work.
            2.  If the first and last element of pb.iteflat(levels) is
                np.NINF and np.PINF, the colorbar of contourf plot will
                use the 'two-arrow' shape.
            3.  If data_transform==True, the input data will be transformed
                from pb.iteflat(levels) to
                np.linspace(1,len(pb.iteflat(interval_original)). this can
                help to create arbitrary contrasting in the plot.
                cf. mathex.plot_array_transg
        data_transform:
            1. set as True if increased contrast in the plot is desired.
                In this case the function mathex.plot_array_transg will
                be called and pb.iteflat(levels) will be used as original
                interval for data transformation.
            2. In case of data_transform==False, pb.iteflat(levels)
                will be used directly in the plt.contour function for
                ploting and hence no data transformation is made. The
                treatment by this way allows very flexible
                (in a mixed way) to set levels.
            3. In any case, if np.NINF and np.PINF as used as two
                extremes of levels, arrowed colorbar will be returned.

        colorbarlabel:
            1. used to put customized colorbar label and this will override
                using levels as colorbar. IF colorbarlabel!=None,
                colorbar ticks and labels will be set using colorbarlabel.
                so this means colorbarlabel could only be array or
                list of numbers.
            2. If data_transform==True, colorbar will also be transformed
                accordingly. In this case, the colorbar ticks will use
                transformed colorbarlabel data, but colorbar ticklables
                will use non-transformed colorbarlabel data. This means
                the actual ticks numbers and labels are not the same.

        forcelabel --> to force the colorbar label as specified by forcelabel.
            This is used in case to set the labels not in numbers but in
            other forms (eg. strings).
            In case of data_transform = True, levels will be used to
            specifiy levels for the original colorbar, colorbarlabel will
            be used to create ticks on colrobar which will be labeled,
            if forcelabel=None, then colorbarlabel will agined be used
            to label the ticks, otherwise forcelabel will be used to
            label the ticks on the colorbar. So this means forcelabel will
            mainly be list of strings.

        indata --> numpy array with dimension of len(lat)Xlen(lon)
        map_threshold --> dictionary like {'lb':2000,'ub':5000}, data
            less than 2000 and greater than 5000 will be masked.
            Note this will be applied before data.
               transform.
        shift --> boolean value. False for longtitude data ranging [-180,180];
            for longtitude data ranging [0,360] set shift to True if a
            180 east shift is desired. if shift as True, the mapbound
            range should be set using shifted longtitude
            (use -180,180 rather than 0,360).
        colorbardic --> dictionary to specify the attributes for colorbar,
            translate all the keys in function bmp.Basemap.colorbar()
            into keys in colorbardic to manipulation.

    Note:
        1.  lat must be descending, and lon must be ascending.
        2*. NOTE use both data_transform=True and impose unequal
            colorbarlabel could be very confusing! Because normaly in
            case of data_transform as True the labels are ALREADY
            UNEQUALLY distributed!

            an example to use colorbarlabel and forcelabel:
            data_transform=True,
            levels=[0,1,2,3,4,5,6,7,8]
            colorbarlabel=[0,2,4,6,8]
            forcelabel=['extreme low','low','middle','high','extreme high']

            So colorbarlabel will set both ticks and labels, but forcelabel
            will further overwrite the labels.

        3. This function has been test using data, the script and
            generated PNG files are availabe at ~/python/bmaptest
    See also:
        mathex.plot_array_transg
    """
    print "Deprecating Warning!: use mapcontourf instead!"
    if ax==None:
        fig,ax=g.Create_1Axes()
    if shift==True:
        indata,lon=bmp.shiftgrid(180,indata,lon,start=False)
    m,lonpro,latpro,latind,lonind=makemap(ax, projection, mapbound,
                                          lat, lon, gridstep)
    pdata=indata[latind[0]:latind[-1]+1,lonind[0]:lonind[-1]+1]
    #mask by map_threshold
    pdata = mathex.ndarray_mask_by_threshold(pdata,map_threshold)
    #prepare the data for contourf
    pdata,plotlev,plotlab,extend,trans_base_list = \
        _transform_data(pdata,levels,data_transform)
    #make the contourf plot
    if cmap==None:
        cmap=mat.cm.jet
    cs=m.contourf(lonpro,latpro,pdata,levels=plotlev,extend=extend,
                  cmap=cmap)

    ##handle colorbar
    #handle the colorbar attributes by using dictionary which flexibility.
    location = colorbardic.get('location','right')
    size = colorbardic.get('size','3%')
    pad = colorbardic.get('pad','2%')
    cbar=m.colorbar(cs,location=location, size=size, pad=pad,**cbarkw)
    #set colorbar ticks and colorbar label
    if levels==None:
        pass
    else:
        ticks,labels = \
            _generate_colorbar_ticks_label(data_transform=data_transform,
                                           colorbarlabel=colorbarlabel,
                                           trans_base_list=trans_base_list,
                                           forcelabel=forcelabel,
                                           plotlev=plotlev,
                                           plotlab=plotlab)
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(labels)
    #return
    if return_lev_lab==True:
        return m,cbar,plotlev,plotlab
    else:
        return m,cbar



