import bsite as bsite
import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import pupynere as pu
import pandas as pa
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

def near5even(datain):
    if datain%5==0:
        dataout=datain
    else:
        if datain/5<0:
            dataout=np.ceil(datain/5)*5
        else:
            dataout=np.floor(datain/5)*5
    return dataout

def makemap(ax,projection='cyl',mapbound='all',lat=None,lon=None,gridstep=(30,30),half_degree=True):
    """
    Purpose: plot the map used for later contour or image plot.
    Note:
        return m,lonpro,latpro,latind,lonind
        1. return m --> map drawed; lonpro/latpro --> lat/lon transferred to projection coords; latind/lonind --> index for lat/lon falling with mapbound
        2. lat must be descending and lon must be ascending.
    Example:
        >>> fig,ax=g.Create_1Axes()
        >>> m,lonpro,latpro,lonind,latind=bmap.makemap(ax,'cyl',mapbound='all',lat=np.arange(89.75,-89.8,-0.5),lon=np.arange(-179.75,179.8,0.5),gridstep=(30,30))
        >>> x,y=m(116,40)  #plot Beijing
        >>> m.scatter(x,y,s=30,marker='o',color='r')
    """
    if projection=='cyl':
        if isinstance(mapbound,dict):
            raise ValueError('cannot use dict for cyl projection')
        elif mapbound=='all':
            lat1=lat[-1]
            lat2=lat[0]
            lon1=lon[0]
            lon2=lon[-1]
            #when the data is of half degree resolution, often the lat1 and lat2 is in the center of the half degree cell, so we need to adjust for the vertices.
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
        m=bmp.Basemap(projection=projection,llcrnrlat=lat1,urcrnrlat=lat2,\
                    llcrnrlon=lon1,urcrnrlon=lon2,resolution='l',ax=ax)
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
            m=bmp.Basemap(projection='npstere',boundinglat=mapbound['blat'],lon_0=mapbound['lon_0'],resolution='l',ax=ax)
            m.drawcoastlines(linewidth=0.7)
            m.fillcontinents(color='0.8',zorder=0)
            if gridstep!=None:
                m.drawparallels(np.arange(mapbound['para0'],91.,gridstep[0]),labels=[1,0,0,0],fontsize=10)
                m.drawmeridians(np.arange(-180.,181.,gridstep[1]),labels=[0,0,0,0],fontsize=10)
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


def contourfmap(ax,lat,lon,indata,projection='cyl',mapbound='all',gridstep=(30,30),shift=False,map_threshold=None,cmap=None,colorbarlabel=None,forcelabel=None,levels=None,data_transform=False,return_lev_lab=False,colorbardic=None):
    """
    Purpose: plot a map on 'cyl' or 'npstere' projection.
    Arguments:
        ax --> An axes instance
        projection --> for now two projections have been added:
            1. 'cyl' -- for global and regional mapping
            2. 'npstere' -- for north polar centered mapping.
        lat,lon --> geographic coordinate variables; lat must be in desceding order and lon must be ascending.
        mapbound --> specify the bound for mapping;
            1. 'cyl'
                tuple containing (lat1,lat2,lon1,lon2); lat1 --> lower parallel; lat2 --> upper parallel; lon1 --> left meridian; lon2 --> right meridian; default
                'all' means plot the extent of input lat, lon coordinate variables;
                for global mapping, set (-90,90,-180,180) or (-90,90,0,360).
            2. 'npstere'
                mapbound={'blat':45,'lon_0':0,'para0':40}
                blat --> boundinglat in the bmp.Basemap method. The souther limit for mapping.
                lon_0 --> center of desired map domain.
                para0 --> souther boundary for parallel ticks, the default norther limit is 90; default longitude 0-360 (or -180-180)
        gridstep --> the step for parallel and meridian grid for the map, tuple containing (parallel_step, meridian_step).
        levels --> default None; levels=[-5,-2,-1,0,1,2,5] ; or levels=[(-10,-4,-2,-1,-0.4),(-0.2,-0.1,0,0.1,0.2),(0.4,1,2,4,10)]. 
            1. Anything that can work as input for function pb.iteflat() will work. 
            2. If the first and last element of pb.iteflat(levels) is np.NINF and np.PINF, the colorbar of contourf plot will use the 'two-arrow' shape.
            3. If data_transform==True, the input data will be transformed from pb.iteflat(levels) to np.linspace(1,len(pb.iteflat(interval_original)). this can 
               help to create arbitrary contrasting in the plot. cf. mathex.plot_array_transg
        data_transform: 
            1. set as True if increased contrast in the plot is desired. In this case the function mathex.plot_array_transg will be called and pb.iteflat(levels)
               will be used as original interval for data transformation. 
            2. In case of data_transform==False, pb.iteflat(levels) will be used directly in the plt.contour function for ploting and hence no data transformation
               is made. The treatment by this way allows very flexible (in a mixed way) to set levels. 
            3. In any case, if np.NINF and np.PINF as used as two extremes of levels, arrowed colorbar will be returned.
        colorbarlabel:
            1. used to put customized colorbar label and this will override using levels as colorbar. IF colorbarlabel!=None, colorbar ticks and labels will be set
               using colorbarlabel. 
            2. If data_transform==True, colorbar will also be transformed accordingly. In this case, the colorbar ticks will use transformed colorbarlabel data, but
               colorbar ticklables will use non-transformed colorbarlabel data.
        indata --> numpy array with dimension of len(lat)Xlen(lon)
        map_threshold --> dictionary like {'lb':2000,'ub':5000}, data less than 2000 and greater than 5000 will be masked. Note this will be applied before data 
               transform.
        shift --> boolean value. False for longtitude data ranging [-180,180]; for longtitude data ranging [0,360] set shift to True if a 180 east shift is desired. 
            if shift as True, the mapbound range should be set using shifted longtitude (use -180,180 rather than 0,360)
        colorbardic --> dictionary to specify the attributes for colorbar, translate all the keys in function bmp.Basemap.colorbar() into keys in colorbardic to 
            manipulation.
        forcelabel --> to force the colorbar label as specified by forcelabel. This is used in case to set the labels not in numbers but in other forms (eg. strings).
            In case of data_transform = True, levels will be used to specifiy levels for the original colorbar, colorbarlabel will be used to create ticks on colrobar
            which will be labeled, if forcelabel=None, then colorbarlabel will agined be used to label the ticks, otherwise forcelabel will be used to label the ticks
            on the colorbar.

    Note:
        1. lat must be descending, and lon must be ascending.
        2*. NOTE use both data_transform=True and impose unequal colorbarlabel could be very confusing! Because normaly in case of data_transform as True the labels
            are ALREADY UNEQUALLY distributed!
        3. This function has been test using data, the script and generated PNG files are availabe at ~/python/bmaptest
    See also:
        mathex.plot_array_transg
    """
    #parameter description: 
    #plotlev: real levels used in plot, incase of data_transform is True, return the pb.iteflat(levels)
    #plotlab: real lable used in colorbar
    if shift==True:
        indata,lon=bmp.shiftgrid(180,indata,lon,start=False)
    m,lonpro,latpro,latind,lonind=makemap(ax, projection, mapbound, lat, lon, gridstep)
    pdata=indata[latind[0]:latind[-1]+1,lonind[0]:lonind[-1]+1]
    #mask by map_threshold
    if map_threshold==None:
        pass
    else:
        if not isinstance(map_threshold,dict):
            raise ValueError('please provide a dictionary for map_threshold')
        else:
            for bound_key in map_threshold.keys():
                if bound_key not in ['lb','ub']:
                    raise ValueError ('Incorrect key is used.')
            if len(map_threshold.keys())==1:
                if map_threshold.keys()[0]=='lb':
                    lower_bound=map_threshold['lb']
                    pdata=np.ma.masked_less(pdata,lower_bound)
                elif map_threshold.keys()[0]=='ub':
                    upper_bound=map_threshold['ub']
                    pdata=np.ma.masked_greater(pdata,upper_bound)
            else:
                lower_bound=map_threshold['lb']
                upper_bound=map_threshold['ub']
                pdata=np.ma.masked_where(np.logical_not((pdata>lower_bound)&(pdata<upper_bound)),pdata)

    if cmap==None:
        cmap=mat.cm.jet
    ###################make the contourf plot
    if levels==None:
        cs=m.contourf(lonpro,latpro,pdata,cmap=cmap)
        if data_transform==True:
            raise Warning("the levels is None but data_transform is True, this is quite strange!")
    else:
        if data_transform==True:
            pdata_trans,plotlev,plotlab,trans_base_list=mathex.plot_array_transg(pdata, levels, copy=True)  #make the data transform before plot.
            if np.isneginf(plotlab[0]) and np.isposinf(plotlab[-1]):
                cs=m.contourf(lonpro,latpro,pdata_trans,levels=plotlev[1:-1],extend='both',cmap=cmap)
            elif np.isneginf(plotlab[0]) or np.isposinf(plotlab[-1]):
                raise ValueError("only one extreme set as infinitive, please set both as infinitive if arrow colorbar is wanted.")
            else:
                cs=m.contourf(lonpro,latpro,pdata_trans,levels=plotlev,cmap=cmap)
        #data_transform==False
        else:
            plotlev=pb.iteflat(levels)
            plotlab=pb.iteflat(levels)
            if np.isneginf(plotlab[0]) and np.isposinf(plotlab[-1]):
                #here the levels would be like [np.NINF,1,2,3,np.PINF]
                #in following contourf, all values <1 and all values>3 will be automatically plotted in the color of two arrows.
                #easy to see in this example:
                #a=np.tile(np.arange(10),10).reshape(10,10);fig,ax=g.Create_1Axes();cs=ax.contourf(a,levels=np.arange(2,7),extend='both');plt.colorbar(cs)
                cs=m.contourf(lonpro,latpro,pdata,levels=plotlev[1:-1],extend='both',cmap=cmap) 
            elif np.isneginf(plotlab[0]) or np.isposinf(plotlab[-1]):
                raise ValueError("only one extreme set as infinitive, please set both as infinitive if arrow colorbar is wanted.")
            else:
                cs=m.contourf(lonpro,latpro,pdata,levels=plotlev,cmap=cmap)
    #################handle colorbar
    #handle the colorbar attributes by using dictionary which allow flexibility.
    colorbardic_default = dict(location='right', size='3%', pad='2%')
    if colorbardic == None:
        colorbardic_used = colorbardic_default
    else:
        if not isinstance(colorbardic,dict):
            raise ValueError("colorbardic must be provided as a dict!")
        else:
            colorbardic_default.update(colorbardic)
            colorbardic_used = colorbardic_default.copy()
    colorbar_kwarg = pb.Dic_Remove_By_Subkeylist(colorbardic_used,['location','size','pad'])
    cbar=m.colorbar(cs,location=colorbardic_used['location'], size=colorbardic_used['size'], pad=colorbardic_used['pad'],**colorbar_kwarg)

    #########set colorbar ticks and colorbar label

    ##here we define a function allowing the easy comparison of forcelable with the number of labeled ticks (ticks that are needed to be labeled)
    def check_forcelabel(forcelabel,tick_array_to_compare):
        """
        This fuction allows to check if the forcelable is given in a correct way when we want to force the label, i.e., in case of using non-numbers,
        or we just want to use strings to denote the colorbar levels.
        """
        if forcelabel == None:
            return tick_array_to_compare
        else:
            if len(forcelabel) != len(tick_array_to_compare):
                raise ValueError("the length of the forcelabel and the length of labeled ticks on the colorbar is not equal!")
            else:
                return forcelabel

    ###here we begins to set colorbar ticks and colorbar labels.
    #this case allows to use default set by matplotlib.
    if levels==None:
        pass
    #data_constrast==True and levels!=None
    elif data_transform==True:
        ##**NOTE THE FOLLOWING IS VERY RARELY USED**##
        #because in this case, data has been transformed, so if we want to impose colorbarlabel, we need to do transformation of colorbarlabel in exactly
        #the same way as data itself. then we set colorbar ticks using transformed colorbarlabel, but set the colorbar ticklabels with original colorbarlabel.
        #*yet notice in this case, the string type of colorbarlabel is not allowed. this can only achieved by explicitly calling cbar.set_ticklabel method.
        #*Then the colorbarlable is actually used as a method to fix the position where we would like to put the string type labels.
        #(comment: 2012-11-01 with the forcelabel flag, the forced label could be ['a','b','c','ff','gg'] etc., but test has not been done yet)
        if colorbarlabel!=None:
            colorbarlabel=pb.iteflat(colorbarlabel)
            transformed_colorbarlabel_ticks,x,y,trans_base_list=mathex.plot_array_transg(colorbarlabel, trans_base_list, copy=True)

        #Noth if/else blocks are organized in 1st tire by check if the two ends are -inf/inf and 2nd tire by check if colorbarlabel is None
        if np.isneginf(plotlab[0]) and np.isposinf(plotlab[-1]):
            if colorbarlabel!=None:
                cbar.set_ticks(transformed_colorbarlabel_ticks)
                cbar.set_ticklabels(check_forcelabel(forcelabel,colorbarlabel))
            else:
                cbar.set_ticks(plotlev[1:-1])
                cbar.set_ticklabels(check_forcelabel(forcelabel,plotlab[1:-1]))
        elif np.isneginf(plotlab[0]) or np.isposinf(plotlab[-1]):
            raise ValueError("It's strange to set only side as infitive")
        else:
            if colorbarlabel!=None:
                cbar.set_ticks(transformed_colorbarlabel_ticks)
                cbar.set_ticklabels(check_forcelabel(forcelabel,colorbarlabel))
            else:
                cbar.set_ticks(plotlev)
                cbar.set_ticklabels(check_forcelabel(forcelabel,plotlab))

    #data_transform==False
    else:
        if np.isneginf(plotlab[0]) and np.isposinf(plotlab[-1]):
            #if colorbarlabel is forced, then ticks and ticklabels will be forced.
            if colorbarlabel!=None:
                cbar.set_ticks(colorbarlabel)
                cbar.set_ticklabels(check_forcelabel(forcelabel,colorbarlabel))
            #This by default will be done, it's maintained here only for clarity.
            else:
                cbar.set_ticks(plotlab[1:-1])
                cbar.set_ticklabels(check_forcelabel(forcelabel,plotlab[1:-1]))
        elif np.isneginf(plotlab[0]) or np.isposinf(plotlab[-1]):
            raise ValueError("It's strange to set only side as infitive")
        else:
            if colorbarlabel!=None:
                cbar.set_ticks(colorbarlabel)
                cbar.set_ticklabels(check_forcelabel(forcelabel,colorbarlabel))
            else:
                cbar.set_ticks(plotlab)
                cbar.set_ticklabels(check_forcelabel(forcelabel,plotlab))
    ##################return
    if return_lev_lab==True:
        return m,cbar,plotlev,plotlab
    else:
        return m,cbar

def contourfmap2(lat,lon,indata,projection='cyl',mapbound='all',gridstep=(30,30),shift=False,cmap=None,map_threshold=None,colorbarlabel=None,levels=None,data_transform=False,ax=None,colorbardic=None):
    """
    contourfmap2 is a wrapper of contourfmap. NO need to set up a figure and axes before drawing map, return fig,ax,m,cbar.
    """
    if ax==None:
        fig,axt=g.Create_1Axes()
    else:
        axt=ax
    m,bar=contourfmap(axt,lat,lon,indata,projection=projection,mapbound=mapbound,gridstep=gridstep,shift=shift,cmap=cmap,map_threshold=map_threshold,\
                      colorbarlabel=colorbarlabel,levels=levels,\
                      data_transform=data_transform,colorbardic=colorbardic)
    if ax==None:
        return fig,axt,m,bar
    else:
        return m,bar

def imshowmap(lat,lon,indata,ax=None,projection='cyl',mapbound='all',gridstep=(30,30),shift=False,colorbar=True,colorbarlabel=None,*args,**kwargs):
    """
    Purpose: plot a map on cyl projection.
    Arguments:
        ax --> An axes instance
        lat,lon --> geographic coordinate variables;
        mapbound --> tuple containing (lat1,lat2,lon1,lon2); lat1 --> lower parallel; lat2 --> upper parallel; lon1 --> left meridian; lon2 --> right meridian; default
            'all' means plot the extent of input lat, lon coordinate variables;
        gridstep --> the step for parallel and meridian grid for the map, tuple containing (parallel_step, meridian_step).
        vmin,vmax --> as in plt.imshow function
        indata --> numpy array with dimension of len(lat)Xlen(lon)
        shift --> boolean value. False for longtitude data ranging [-180,180]; for longtitude data ranging [0,360] set shift to True if a 180 east shift is desired. 
    """
    #handle the case ax==None:
    if ax==None:
        fig,axt=g.Create_1Axes()
    else:
        axt=ax

    if shift==True:
        indata,lon=bmp.shiftgrid(180,indata,lon,start=False)

    #make the map and use mapbound to cut the data
    m,lonpro,latpro,latind,lonind=makemap(axt, projection, mapbound, lat, lon, gridstep)
    pdata=indata[latind[0]:latind[-1]+1,lonind[0]:lonind[-1]+1]

    cs=m.imshow(pdata,origin='upper',*args,**kwargs)
    if colorbar==True:
        cbar=m.colorbar(cs)
        if colorbarlabel!=None:
            cbar.set_label(colorbarlabel)
    else:
        cbar=None
    return m,cs,cbar


