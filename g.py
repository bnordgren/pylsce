import bsite as bsite
import matplotlib
import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
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
import pb
import rpy2
import copy
#import rpy


'''
    1. Many auxiliary plot functions are defined.
    2. generally used labels are defined.

'''
pline=['b-','g-','r-','c-','m-','y-','k-']
pdot=['b.','g.','r.','c.','m.','y.','k.']
pcircledot=['bo','go','ro','co','mo','yo','ko']
pcolor=['b','g','r','c','m','y','k','#cd5c5c','#ffa500','#ff69b4','#4682b4','#191970','#7cfc00','#ffff00','#bdb76b']
spt8=[421,422,423,424,425,426,427,428]
spt4=[221,222,223,224]
c2={    'db':'#000080',
        'b':'#0000BF',
        'lb':'#8080FF',
        'dr':'#800000',
        'r':'#BF0000',
        'lr':'#FF8080',
        'dor':'#804000',
        'or':'#BF6000',
        'lor':'#FFBF80',
        'dg':'#008000',
        'g':'#00BF00',
        'lg':'#80FF80'
       }




chtml=['Red','White','Cyan','Silver','Blue','Grey','DarkBlue','Black','LightBlue','Orange','Purple','Brown','Yellow','Maroon','Lime','Green','Fuchsia','Olive']
def show_html_color():
    fig,ax=Create_1Axes()
    for i in range(len(chtml)):
        ax.plot(np.arange(10)+2*i,chtml[i],lw=3)
        ax.text(1,1+2*i,chtml[i])
    plt.show()

def show_html_color_all():
    fig,axt=plt.subplots(ncols=5,nrows=1,subplot_kw=dict(xticks=[],yticks=[]),figsize=(11,9))
    for i,(name,color) in enumerate(mat.colors.cnames.items()):
        ax=axt[i/29]
        ax.plot(np.ones(10)+3*i,color,lw=5)
        ax.text(1,3*i,name,alpha=0.5)
    plt.show()

def show_cmap(cmap=mat.cm.jet,levnum=10):
    """
    show the effect of a colormap
    """
    dt = np.arange(1,levnum+1)
    dt = np.tile(dt[:,np.newaxis][::-1],(1,levnum))
    plt.imshow(dt,cmap=cmap)
    plt.colorbar()

def Set_Axes_Cross_Line(ax,color='r',linewidth=None):
    l1 = mat.lines.Line2D([0,1],[1,0],transform=ax.transAxes,color=color,linewidth=linewidth,axes=ax)
    l2 = mat.lines.Line2D([0,1],[0,1],transform=ax.transAxes,color=color,linewidth=linewidth,axes=ax)
    ax.lines.extend([l1,l2])


def c2show():
    for i in range(len(c2)):
        cname=c2.keys()[i]
        plt.plot(np.arange(10)+i,c2[cname])
        plt.text(1,1+i,cname)
    plt.show()

def pcolor_show():
    for i in range(len(pcolor)):
        plt.plot(np.arange(10)+i,pcolor[i])
        plt.text(1,1+i,str(i))
    plt.show()
def rc():
    """
    get a random color from pcolor list
    """
    return pcolor[np.random.random_integers(1,len(pcolor))-1]

bluecmp=mat.colors.LinearSegmentedColormap.from_list('bluecmp',['#ADD6FF','#3366FF','#3333CC','#000066'])
def blue_level(levnum):
    a=bluecmp(np.linspace(0,1,num=levnum,endpoint=True))
    return a


def extract_color_list_cm(cmap,levnum=10):
    if not isinstance(cmap, mat.colors.Colormap):
        raise TypeError("provided colormap is not mat.colors.Colormap object!")
    else:
        rgba_array = cmap(np.linspace(0,1,num=levnum,endpoint=True))
        return rgba_array[...,0:3]


def cm_extract(cmap,(start,end),levnum=10,levels=None):
    """
    Create a new cmap by using the matplotlib provided colormaps, but using
        only the part specified by (start,end).
    -----------
    Parameters:
    cmap : a mat.cm instances.
    (start,end) : the start and end number of level that will be used to
        generate the new colormap. For example, with (start,end)= (3,10)
        and levnum=10, the original colormap will be devided into 10 levels,
        and the colors between level3 and level10 will be used to
        construct the new colormap.
        Note the 'start' level is always included and it's 1-based.
    levnum : the number of levels that original colormap will be divided into.
    levels : cf. function g.rgb2cm
    --------
    Doctest:
        >>> g.show_cmap(mat.cm.gist_ncar,30)
        >>> g.show_cmap(g.cm_extract(mat.cm.gist_ncar,(15,27),30))
        >>> g.show_cmap(g.cm_extract(mat.cm.gist_ncar,(6,18),30))

    """
    if not isinstance(cmap, mat.colors.Colormap):
        raise TypeError("provided colormap is not mat.colors.Colormap object!")
    else:
        rgba_array = cmap(np.linspace(0,1,num=levnum,endpoint=True))
        extract_rgba_array_255 = rgba_array[start-1:end,0:3]*255.
        return rgb2cm(extract_rgba_array_255, cmname='tempcm', levels=levels)

def cm_concat_multiple_cm(concat_cmap_list):
    """
    Create a new customized colormap by concatenating each part of multiple
        colormaps.

    Arguments:
    ----------
    concat_cmap_list: a list of tuples, with each tuple as
        (cmap,frac_start,frac_end,slice_number); in which,
        cmap: name of colormap that's to be extracted from.
        frac_start/frac_end: start and end positions of the camp that's to
            be extracted and concatenated with others. frac_start and frac_end
            must be fractions (0-1), and frac_end is always included.
        slice_number: the number of slices to be extracted from the
            [frac_start,frac_end] of the concerned colormap.

    Notes:
    ------
    For further example, see $PYLIB/function_doc_notebook/colorbar_color_change_around_zero.ipynb
    """
    final_rgba_array_list = []
    for (cmap,frac_start,frac_end,slice_number) in concat_cmap_list:
        if not isinstance(cmap, mat.colors.Colormap):
            raise TypeError("{0} is not mat.colors.Colormap object!".format(cmap))
        else:
            rgba_array = cmap(np.linspace(frac_start, frac_end,
                                          num=slice_number, endpoint=True))
            extract_rgba_array_255 = rgba_array[:,0:3]*255.
            final_rgba_array_list.append(extract_rgba_array_255)
    final_rgb_array = np.concatenate(final_rgba_array_list,axis=0)
    return rgb2cm(final_rgb_array, cmname='tempcm')


redcmp=mat.colors.LinearSegmentedColormap.from_list('redcmp',['#FFCCCC','#330000'])
def red_level(levnum):
    a=redcmp(np.linspace(0,1,num=levnum,endpoint=True))
    return a

greencmp=mat.colors.LinearSegmentedColormap.from_list('greencmp',['#ADEBAD','#0F3D0F'])
def green_level(levnum):
    a=greencmp(np.linspace(0,1,num=levnum,endpoint=True))
    return a

#define 7 blue clolors
#pcblue7=[(0.0039215686274509803, 0.082352941176470587, 0.32549019607843138),
# (0.0078431372549019607, 0.16470588235294117, 0.66666666666666663),
# (0.011764705882352941, 0.23137254901960785, 0.93725490196078431),
# (0.18431372549019609, 0.37647058823529411, 0.99215686274509807),
# (0.46666666666666667, 0.59215686274509804, 0.99215686274509807),
# (0.72156862745098038, 0.78823529411764703, 0.99215686274509807),
# (0.90588235294117647, 0.92156862745098034, 1.0)]
pcblue7=['#0066FF','#0099CC','#3333CC','#6600CC','#003399','#000066','#66CCFF']

#define 3 green colors
pcgreen3=['#009933','#33CC33','#00FF00']
#pcgreen3=[(0.054901960784313725, 0.22745098039215686, 0.10196078431372549),
# (0.16078431372549021, 0.63921568627450975, 0.28627450980392155),
# (0.60392156862745094, 0.90196078431372551, 0.68235294117647061)]
#define 3 red colors
pcred3=['#FF0000','#FF6666','#FF3399']
#pcred3=[(0.41960784313725491, 0.066666666666666666, 0.058823529411764705),
# (0.8784313725490196, 0.13725490196078433, 0.11764705882352941),
# (0.95294117647058818, 0.66274509803921566, 0.65490196078431373)]
#define 7 red colors
pcred7=[(0.41960784313725491, 0.066666666666666666, 0.058823529411764705),
 (0.8784313725490196, 0.13725490196078433, 0.11764705882352941),
 (0.95294117647058818, 0.66274509803921566, 0.65490196078431373),
 (0.95294117647058818, 0.66274509803921566, 0.65490196078431373),
 (0.95294117647058818, 0.66274509803921566, 0.65490196078431373),
 (0.95294117647058818, 0.66274509803921566, 0.65490196078431373),
 (0.95294117647058818, 0.66274509803921566, 0.65490196078431373)]
#temp 13 colors for 13 sites
pctemp13=[(0.0039215686274509803, 0.082352941176470587, 0.32549019607843138),
 (0.0078431372549019607, 0.16470588235294117, 0.66666666666666663),
 (0.011764705882352941, 0.23137254901960785, 0.93725490196078431),
 (0.18431372549019609, 0.37647058823529411, 0.99215686274509807),
 (0.46666666666666667, 0.59215686274509804, 0.99215686274509807),
 (0.72156862745098038, 0.78823529411764703, 0.99215686274509807),
 (0.90588235294117647, 0.92156862745098034, 1.0),
 (0.054901960784313725, 0.22745098039215686, 0.10196078431372549),
 (0.16078431372549021, 0.63921568627450975, 0.28627450980392155),
 (0.60392156862745094, 0.90196078431372551, 0.68235294117647061),
 (0.41960784313725491, 0.066666666666666666, 0.058823529411764705),
 (0.8784313725490196, 0.13725490196078433, 0.11764705882352941),
 (0.95294117647058818, 0.66274509803921566, 0.65490196078431373)]

#define a color dictionary for the 13 sites:CA-NS-->blue; CA-SF1-->green; US-Bn-->red
pc13=dict()
for i,site in enumerate(bsite.NSlist):
    pc13[site]=pcblue7[i]
for i,site in enumerate(bsite.SFlist):
    pc13[site]=pcgreen3[i]
for i,site in enumerate(bsite.Bnlist):
    pc13[site]=pcred3[i]


def set_default_ax(ax=None):
    if ax == None:
        fig,ax=Create_1Axes()
    else:
        pass
    return ax


def showrgb(rgblistin,plot=True,raw=True,getdecimal=False,copy=True,ax=None):
    """
    Purpose: show color of RGB(A) list or array. input as RGB list [(0.4,0.1,0.7)], or raw RGB list [(110,123,139)] (raw as True)
             set getdecimal as True to return the 0-1 RGB list; set plot as False to suppress the plot.
    Arguments:
        rbglistin: 1. a list of RGB tuples, eg. [(0.4,0.1,0.7),(0.4,0.5,0.7)], when raw==True,[(110,123,139),(110,123,255)] could be used.
                   2. a numpy array of RGB values (nX3) array, with R/G/B as 1st/2nd/3rd column.
                   3. a numpy array of RGBA values (nX3) array, with R/G/B/A as 1st/2nd/3rd/4th column.
                   3. when raw=False, rgblistin could be a nX4 ndarray with rows as RGBA values.
    """

    if isinstance(rgblistin, np.ndarray):
        #RGB array might be provided as 2dim ndarray, either in form of RGB array (nX3) or in form of RGBA array(nX4)
        if rgblistin.ndim ==2:
            #to make it able to recieve a nX4 RGBA array, we need to only extract the nX3 RGB values from the RGBA array.
            if rgblistin.shape[1] == 4:
                rgbarray = rgblistin[:,0:3]
            elif rgblistin.shape[1] == 3:
                rgbarray = rgblistin[:]
            else:
                raise ValueError("The provided color array is 2dim, but the second dim length is not 3 or 4, but {0}".format(rgblistin.shape[1]))
            # to recieve the rgb value as a nX3 ndarray with each row as RGB value, we need to change the nx3 array into rgblist
            newlist = [tuple(row) for row in rgbarray]
            rgblistin = newlist[:]
        else:
            raise ValueError("the provided color array is not 2dim, but rather {0} dim".format(rgblistin.ndim))
    elif isinstance(rgblistin, list):
        pass
    else:
        raise TypeError("rgblistin must be provided as list or ndarray, the provided type is {0}".format(type(rgblistin)))

    #set default axes
    if plot:
        ax = set_default_ax(ax)

    if copy==True:
        rgblist=rgblistin[:]
    else:
        rgblist=rgblistin

    n=len(rgblist)
    if raw==False:
        if plot:
            for i in range(1,n+1):
                ax.plot(i+0.5,i+0.5,color=rgblist[i-1],marker='s',ms=20)
                ax.text(i+0.5,i+0.7,str(i))
    else:
        for i in range(len(rgblist)):
            rgblist[i]=tuple([round(num/255.,2) for num in rgblist[i]])
        if plot:
            for i in range(1,n+1):
                ax.plot(i+0.5,i+0.5,color=rgblist[i-1],marker='s',ms=20)
                ax.text(i+0.5,i+0.7,str(i))
    if raw==True and getdecimal==True:
        return rgblist

def rgb2cmdic(rgblist,levels=None):
    """
    Purpose: transfer RGB colorlist to mat.cm dictionary as indicated by levels. input should be RGB values 0-255.
    Arguments:
        levels: the levels that used to corresond the rgblist. if None is provided, np.linspace(0,1,len(rgblist),endpoint=True) is used.
        rgblist: RGB tuple list or nX3 ndarray, cf. showrgb for more information.
    """
    rgblist=showrgb(rgblist,plot=False,raw=True,getdecimal=True)
    rgbarray=np.array(rgblist)
    if levels == None:
        levels=np.linspace(0,1,len(rgblist),endpoint=True)
    elif levels[0] !=0 or levels[-1] !=1:
        levels = mathex.lintrans(np.array(levels),(levels[0],levels[-1]),(0,1))
    else:
        pass
    def getcolordic(levels,b,colnum):
        if len(b[:,colnum])<3:
            raise ValueError('Need at least 3 points to do the transformation')
        else:
            red0=b.copy()
            red0[:,0]=np.array(levels)
            red0[0][1]=0.
            red0[-1][2]=0.
            red0[0][2]=b[:,colnum][0]
            red0[-1][1]=b[:,colnum][-1]
            for i in range(1,len(b[:,colnum])-1):
                red0[i,:][1:3]=b[:,colnum][i]
            redlist=[]
            for i in red0:
                redlist.append(tuple([round(j,2) for j in i]))
            return redlist
    colorlist=['red','green','blue']
    cmdict=dict()
    for colnum,col in enumerate(colorlist):
        cmdict[col]=getcolordic(levels,rgbarray,colnum)
    return cmdict

def rgb2cm(rgblist,cmname='tempcm',levels=None):
    """
    Purpose: use a rgblist to construct a customised colormap. len(levels) must be equal to len(rgblist).
    Arguments:
        rgblist: RGB tuple list or nX3 ndarray (range 0-255), cf. showrgb for more information.
    """
    cmdic = rgb2cmdic(rgblist,levels=levels)
    cm = mat.colors.LinearSegmentedColormap(cmname, cmdic, 256)
    return cm

#most frequently used labels
glab={  'nep':'NEP (gC m$^{-2}$yr$^{-1}$)',
        'gpp':'GPP (gC m$^{-2}$yr$^{-1}$)',
        're':'Ecosystem respiration (gC m$^{-2}$yr$^{-1}$)',
        'npp':'NPP (gC m$^{-2}$yr$^{-1}$)',
        'anpp':'Aboveground NPP (gC m$^{-2}$yr$^{-1}$)',
        'totalb':'Total biomass carbon (gC m$^{-2}$)',
        'belowb':'Belowground Biomass Carbon (gC m$^{-2}$)',
        'aboveb':'Aboveground Biomass Carbon (gC m$^{-2}$)',
        'forestf':'Forest floor carbon (gC m$^{-2}$)',
        'abovel':'aboveground litter carbon (gcm$^{-2}$)',
        'cwd':'Coarse Woody Debris carbon (gcm$^{-2}$)',
        'minerals':'Mineral soil carbon (gC m$^{-2}$)',
        'NEP_YEAR':'Net Ecosystem Production (gCm$^{-2}$yr$^{-1}$)',
        'NEP':'Net Ecosystem Production (gCm$^{-2}$yr$^{-1}$)',
        'GPP_YEAR':'GPP (gC m$^{-2}$yr$^{-1}$)',
        'RECO_YEAR':'Ecosystem respiration (gC m$^{-2}$yr$^{-1}$)',
        'HET_RESP_YEAR':'Heterotrophic respiration (gC m$^{-2}$yr$^{-1}$)',
        'AUTO_RESP_YEAR':'Autotrophic respiration (gC m$^{-2}$yr$^{-1}$)',
        'MAINT_RESP_YEAR':'Maintenance respiration (gC m$^{-2}$yr$^{-1}$)',
        'NPP_YEAR':'NPP (gC m$^{-2}$yr$^{-1}$)',
        'anpp':'Aboveground NPP (gC m$^{-2}$yr$^{-1}$)',
        'TOTAL_M':'Total biomass carbon (gC m$^{-2}$)',
        'ROOT_M':'Root biomass carbon (gC m$^{-2}$)',
        'TOTAL_M_BE':'Belowground Biomass Carbon (gC m$^{-2}$)',
        'TOTAL_M_AB':'Aboveground Biomass Carbon (gC m$^{-2}$)',
        'FOREST_FLOOR_CARB':'Forest Floor Carbon (gC m$^{-2}$)',
        'LITTER_AB':'Aboveground Litter Carbon (gCm$^{-2}$)',
        'LITTER_WOD_AB':'Coarse Woody Debris Carbon (gCm$^{-2}$)',
        'MINERAL_SOIL_CARB':'Mineral soil carbon (gC m$^{-2}$)',
        'TOTAL_SOIL_CARB':'Total soil carbon (gC m$^{-2}$)',
        'SNAG_AB':'Snag aboveground (gC m$^{-2}$)',
        'SNAG_BE':'Snag belowground (gC m$^{-2}$)',
        'IND':'Individual density (1 ha$^{-1}$)',
        'LAI':'Leaf Area Index (m2m$^{-2}$)',
        'LEAF_M':'Leaf Biomass Carbon (gC m$^{-2}$)',
        'BA':'Basal area(m2ha$^{-1}$)',
        'DBH':'Diameter at Breath Height (cm)'
     }

def FilterStringList(keyword,input_list):
    return [x for x in input_list if re.search(keyword,x)]

def FilterArtistLabel(ObjectInstanceList,*labelvalues):
    """
    Purpose: Filter a list of object instances by their attribute value; supply with (attribute_name,attribute_value) tuples
    """
    select_bool=[True]*len(ObjectInstanceList)
    for i,artist in enumerate(ObjectInstanceList):
        for labelvalue in labelvalues:
            try:
                if labelvalue in artist._label:
                    pass
                else:
                    select_bool[i]=False
                    break
            except:
                pass
    return [obinstance for boolflag,obinstance in zip(select_bool,ObjectInstanceList) if boolflag==True]

def ListSingleValue(d):
    flag=True
    for i in d:
        if i!=d[0]:
            flag=False
            break
    return flag

def StringListAnotB(listA,listB):
    return [i for i in listA if i not in listB]

def StringListAandB(listA,listB):
    return [i for i in listA if i in listB]

def PdfFile(filename):
    """
    Creat and return a pdf backend for ploting.
    """
    pp = PdfPages(filename)
    return pp

def Remove_dupTL(input_list_or_tuple,remove_element):
    if isinstance(input_list_or_tuple,tuple):
        templist=list(input_list_or_tuple)
        for i in range(templist.count(remove_element)):
            templist.remove(remove_element)
        return tuple(templist)
    elif isinstance(input_list_or_tuple,list):
        templist=input_list_or_tuple
        for i in range(templist.count(remove_element)):
            templist.remove(remove_element)
        return templist

def Remove_dupdim(input_numpy_array):
    or_shape=input_numpy_array.shape
    new_reshape=Remove_dupTL(or_shape,1)
    return np.reshape(input_numpy_array,new_reshape)
class ncdata(object):
    pass

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

def ncread(filename):
    """
    Purpose: read a .nc file using netCDF4 package and read the data into netCDF4.Variable object
    Definition: ncread(file)
    Arguments:
        file--> file name
    Return: return a list of ncdata object; the 1st one contains original nctCDF4 objects and the 2nd one contains data with duplicate dimensions removed.
    Note:
        1. NEP and AUTO_RESP, BM_ALLOC_TOTAL, if not existing and can be calculated from input files, are automatically calculated.

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


    #calculate NEP and AUTO_RESP
    if 'NEP' not in datanew.__dict__.keys():
        if 'NPP' in f1.variables.keys() and 'HET_RESP' in f1.variables.keys():
            datanew.__dict__['NEP']=datanew2.__dict__['NPP']-datanew2.__dict__['HET_RESP']
            datanew2.__dict__['NEP']=datanew2.__dict__['NPP']-datanew2.__dict__['HET_RESP']
    if 'AUTO_RESP' not in datanew.__dict__.keys():
        if 'MAINT_RESP' in f1.variables.keys() and 'GROWTH_RESP' in f1.variables.keys():
            datanew.__dict__['AUTO_RESP']=datanew2.__dict__['MAINT_RESP']+datanew2.__dict__['GROWTH_RESP']
            datanew2.__dict__['AUTO_RESP']=datanew2.__dict__['MAINT_RESP']+datanew2.__dict__['GROWTH_RESP']

    #make the total allocated biomass
    alloc_flag=True
    alloc_list=['BM_ALLOC_FRUIT','BM_ALLOC_HEART_AB','BM_ALLOC_ROOT','BM_ALLOC_SAP_AB','BM_ALLOC_HEART_BE','BM_ALLOC_LEAF','BM_ALLOC_SAP_BE','BM_ALLOC_RES']
    alloc_shape=[]
    for var in alloc_list:
        if var not in f1.variables.keys():
            print var,'not in input file!'
            alloc_flag=False
            break
        else:
            alloc_shape.append(datanew2.__dict__[var].shape)
    if alloc_flag and ListSingleValue(alloc_shape):
        tempvar=np.zeros_like(datanew2.__dict__['BM_ALLOC_FRUIT'])
        for var in alloc_list:
            tempvar=tempvar+datanew2.__dict__[var]
        datanew.__dict__['BM_ALLOC_TOTAL']=tempvar
        datanew2.__dict__['BM_ALLOC_TOTAL']=tempvar
    return [datanew,datanew2]
    f1.close()

def duplicates(lst, item):
    return [i for i,x in enumerate(lst) if x==item]
def Get_PointValue(data,var,(lat,vlat),(lon,vlon)):
    """
    Purpose: get the value of variable 'var' for the point(vlat,vlon) from the grid.
    Definition: Get_PointValue(data,var,(lat,vlat),(lon,vlon))
    Arguments:
        data --> a ncdata object; var --> variable name for which value needs to be retrieved, data.var should be a netCDF4.Variable object.
        lat --> latitude coordinate; lon --> longitude coordinate
    Return: a ndarray of var specified by location
    Note: 
        1. the lat and lon must be the two last dimensions of variable 'var'.
        2. the lat must be a decreasing sequence while the lon must a increasing sequence.
        3. only applies to extract value for a single point
        4. var, lat, and lon must be provided as str
    Example:
        >>> data
         <g.ncdata object at 0x32b6990>
        >>> data.air
         <netCDF4.Variable object at 0x18b3cd0>
        >>> data.air.dimensions
         (u'time', u'lat', u'lon')
        >>> data.air.shape
         (1464, 94, 192)
        >>> temp=g.Get_GridValue(data,'air',('lat',54.49),('lon',150.65))
        >>> temp.shape
         (1464,)
    """

    lat=data.__dict__[lat][:]
    lon=data.__dict__[lon][:]
    if lat.ndim == 2:
        lat=lat[:,0]
    else:
        pass
    if lon.ndim == 2:
        lon=lon[0,:]
    else:
        pass

    index_more=np.where(lat>vlat)[0][-1]
    index_less=np.where(lat<vlat)[0][0]
    if abs(lat[index_more]-vlat) >= abs(lat[index_less]-vlat):
        index_lat=index_less
    else:
        index_lat=index_more
    index_more=np.where(lon>vlon)[0][0]
    index_less=np.where(lon<vlon)[0][-1]
    if abs(lon[index_more]-vlon) >= abs(lon[index_less]-vlon):
        index_lon=index_less
    else:
        index_lon=index_more
    target=data.__dict__[var][:][...,index_lat,index_lon]
    return target

def Get_GridValue(data,var,(lat,vlat1,vlat2),(lon,vlon1,vlon2)):
    """
    Purpose: get the value of variable 'var' for a subgrid (vlat1:vlat2,vlon1:vlon2) from the grid.
    Definition: Get_GridValue(data,var,(lat,vlat1,vlat2),(lon,vlon1,vlon2))
    Arguments:
        data --> a ncdata object; var --> variable name for which value needs to be retrieved, data.var should be a netCDF4.Variable object.
        lat --> latitude coordinate; lon --> longitude coordinate
        vlat1-->vlat2 (South-->North); vlon1-->vlon2(West-->East)
    Return: a ndarray of var specified by location
    Note: 
        1. the lat and lon must be the two last dimensions of variable 'var'.
        2. the lat must be a decreasing sequence while the lon must a increasing sequence.
        3. only applies to extract value for a single point
        4. var, lat, and lon must be provided as str
    Example:
        >>> ba=g.ncread('/home/orchidee01/ychao/FIRE_DATA/GFED3/monthly/BA2009.nc')
        >>> ba.ba.shape
         (12, 360, 720)
        >>> subba=g.Get_GridValue(ba,'ba',('lat',55,56),('lon',-99,-98))
        >>> subba.shape
         (12, 2, 2)

    """
    lat=data.__dict__[lat][:]
    lon=data.__dict__[lon][:]
    if lat.ndim == 2:
        lat=lat[:,0]
    else:
        pass
    if lon.ndim == 2:
        lon=lon[0,:]
    else:
        pass
    index_lat=np.nonzero((lat[:]>=vlat1)&(lat[:]<=vlat2))[0]
    index_lon=np.nonzero((lon[:]>=vlon1)&(lon[:]<=vlon2))[0]
    target=data.__dict__[var][:][...,index_lat[0]:index_lat[-1]+1,index_lon[0]:index_lon[-1]+1]
    return target

def Set_Equal_XYlim(ax,diagline=True,c='k',ls='--',aspect_dic={'anchor':None},**kwargs):
    ax.set_aspect('equal',**aspect_dic)
    xlimit=ax.get_xlim()
    ylimit=ax.get_ylim()
    ax.set_xlim((min(xlimit[0],ylimit[0]),max(xlimit[1],ylimit[1])))
    ax.set_ylim((min(xlimit[0],ylimit[0]),max(xlimit[1],ylimit[1])))
    if diagline==True:
        x=np.linspace(ax.get_xlim()[0],ax.get_xlim()[1],20,endpoint=True)
        ax.set_autoscale_on(False)
        t=ax.plot(x,x,color=c,linestyle=ls,linewidth=0.8,**kwargs)
        ax.set_autoscale_on(True)
        return t

def Set_xaxis_month2year(ax,BeginYear):
    """
    Purpose: set xaxis from month to year, ie. the monthly data are used to draw plots but the xaxis is displayed as year.
             So in this case the xaxis must be integers. 
    Note:
        1. BeginYear, the beginning year used to draw the xtickslabels.
    """
    x0=ax.get_xlim()[0]
    x1=ax.get_xlim()[1]
    numyear=np.floor((x1-x0)/12)
    x1new=numyear*12-1
    ax.set_xlim((x0,x1new))
    ax.set_xticks(x0+5.5+np.arange(numyear)*12)
    xlabelarray=BeginYear+np.arange(numyear)
    ax.set_xticklabels(xlabelarray.astype(int))

def Set_Axes_Position(ax,pos,frac):
    """
    Purpose: Adjust the width or height of an axes instance
    Definition: Set_Axes_Position(ax,pos,frac)
    Arguments:
        ax --> An axes instance
        pos --> 'w' for width adjustment; 'h' for height adjustment
        frac --> fraction of original width or height

    """
    box=ax.get_position()
    if pos=='h':
        ax.set_position([box.x0, box.y0, box.width, box.height*frac])
    elif pos=='w':
        ax.set_position([box.x0, box.y0, box.width*frac, box.height])
    else:
        print "use 'w' for width adjustment and 'h' for height adjustment"

def Axes_Set_Axis_Locater(axes,major=None,minor=None,axis='x'):
    """
        Set x/y axis major and minor locators.
    """
    if axis == 'x':
        axes.xaxis.set_major_locator(plt.MultipleLocator(major))
        axes.xaxis.set_minor_locator(plt.MultipleLocator(minor))
    if axis == 'y':
        axes.yaxis.set_major_locator(plt.MultipleLocator(major))
        axes.yaxis.set_minor_locator(plt.MultipleLocator(minor))

def Set_FigText(fig,figtext):
    fig.text(0.5,0.95,figtext,horizontalalignment='center',va='center')

def Set_Figxlabel(fig,figtext,pos=(0.5,0.04),ha='center'):
    x=pos[0]
    y=pos[1]
    return fig.text(x,y,figtext,ha=ha,rotation='horizontal')
    
def Set_Figylabel(fig,figtext,pos=(0.08,0.5),va='center'):
    x=pos[0]
    y=pos[1]
    return fig.text(x,y,figtext,va=va,ha='center',rotation='vertical')

def Fig_Save_Close(fig,figname):
    fig.savefig(figname)
    plt.close(figname)


def Set_AxText(ax,text_string,pos='uc',ftdic={'size':12},**kwargs):
    """
    Purpose: add a text to axes. pos can be['lc','uc','ll','ul','lr','ur'].
    Example:
        >>> g.Set_AxText_eg()
        >>> plt.show()
    """
    if pos==None or pos==False:
        pass
    elif pos=='lc':
        ax.text(0.5, 0.05,text_string,horizontalalignment='center',verticalalignment='center',transform = ax.transAxes,fontdict=ftdic,**kwargs)
    elif pos=='ouc':
        ax.text(0.5, 1.02,text_string,horizontalalignment='center',verticalalignment='bottom',transform = ax.transAxes,fontdict=ftdic,**kwargs)
    elif pos=='uc':
        ax.text(0.5, 0.93,text_string,horizontalalignment='center',verticalalignment='center',transform = ax.transAxes,fontdict=ftdic,**kwargs)
    elif pos=='ll':
        ax.text(0.02, 0.05,text_string,horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontdict=ftdic,**kwargs)
    elif pos=='ul':
        ax.text(0.02, 0.95,text_string,horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontdict=ftdic,**kwargs)
    elif pos=='lr':
        ax.text(0.98, 0.05,text_string,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes,fontdict=ftdic,**kwargs)
    elif pos=='ur':
        ax.text(0.98, 0.93,text_string,horizontalalignment='right',verticalalignment='center',transform = ax.transAxes,fontdict=ftdic,**kwargs)
    else:
        ax.text(pos[0], pos[1],text_string,transform = ax.transAxes,fontdict=ftdic,**kwargs)

def Set_Insetaxes(fig,ax,pos='ul',xpad=None,ypad=None,width=None,height=None,inset='new',**kwargs):
    """
    Purpose: Create an inset axes (or set its position) within the axes 'ax' inside figure 'fig', return the created axes.
    Arguments:
        pos --> position of inset axes, available now 'ul','ur','ll','lr'
        xpad,ypad --> distance of new axes xaxis to its parent axes xaxis. and ypad new yaxis to parent yaxis. In case of upper location, xpad means upper xaxis of 
            parent axes; In case of right location, ypad means right yaxis of parent axes.
        width,height --> width and height of axes to be created.
        inset --> 'new' for create a new inset axes; or when isinstance (inset,matplotlib.axes.Axes), reset the position of this inset axes position.
        kwargs --> kwagrs in fig.add_axes function.
    Note:
        xpad,ypad,width,height are in fractions of axes that's going to host the inset axes.
    """
    w=width
    h=height
    #This is from http://old.nabble.com/Adding-custom-axes-within-a-subplot-td22159536.html
    #create inset axes at (l1,b1) with width w1 and height h1 inside axes "ax"; l1,b1,w1,h1 in fraction of ax axes.
    def create_inset_axes(fig,ax,l1,b1,w1,h1,**kwargs):
        Bbox = matplotlib.transforms.Bbox.from_bounds(l1,b1,w1,h1)
        trans = ax.transAxes + fig.transFigure.inverted()
        l, b, w, h = matplotlib.transforms.TransformedBbox(Bbox, trans).bounds
        print "the figure fraction for the new inset axes is, l:{0} b:{1} w:{2} h:{3}".format(l,b,w,h)
        axins = fig.add_axes([l, b, w, h],**kwargs)
        return axins
    def cal_coordinate(pos,xpad,ypad,w,h):
        """
        return the lowerleft coordinates (x,y) for the axes to be created.
        """
        x0=y0=0.
        x1=y1=1.
        if pos=='ul':
            return [x0+ypad,y1-xpad-h]
        elif pos=='ur':
            return [x1-ypad-w,y1-xpad-h]
        elif pos=='ll':
            return [x0+ypad,y0+xpad]
        elif pos=='lr':
            return [x1-ypad-w,y0+xpad]
        else:
            raise ValueError('position {0} is not included in current function').format(pos)
    inset_llcor=cal_coordinate(pos,xpad,ypad,w,h)
    if inset=='new':
        l1,b1,w1,h1=inset_llcor+[w,h]
        return create_inset_axes(fig,ax,l1,b1,w1,h1,**kwargs)
    elif isinstance(inset,matplotlib.axes.Axes):
        inset.set_position(inset_llcor+[w,h])
        return inset_llcor

def Calc_Newaxes_Fraction(posor, split_fraction, direction='vertical'):
    """
        To calculate the the x0,y0,width,height for the new axes that will take the place of the old axes, which position is given by argument posor.

    Parameters
    ----------
    posor : list
        a list, denoting the position of original axes that will be replaced by new axes. posor = [x0, y0, width, height]
    split_fraction : list
        indicate how the space of original axes is split among the new axes and the space between these new axes. its length should be number_of_new_axes
        + number_of_new_axes-1, as there are "number_of_new_axes-1" blank spaces between these new axes.

    Returns
    -------
    out : a nx4 numpy array, with "n" as the number of new axes that are inserted. Each row indicates the new axes position: [x0, y0, width, height]

    Example
    -------
    >>> Calc_Newaxes_Fraction([0,0,1,1],[0.2, 0.1, 0.3, 0.1, 0.3])
        array([[ 0. ,  0. ,  1. ,  0.2],
               [ 0. ,  0.3,  1. ,  0.3],
               [ 0. ,  0.7,  1. ,  0.3]])

    """
    x0, y0, width, height = posor
    newpos_list = []
    if np.sum(np.array(split_fraction)) != 1:
        raise ValueError("The split_fraction list sum is not 1")
    else:
        newfrac = np.cumsum(np.array(split_fraction))
        #newaxes_frac_array is a nX2 array, with each row indicating the fraction (to original axes that will be replaced) of the bottom and top axis.
        newaxes_frac_array = np.hstack((np.array([0]),newfrac)).reshape(-1,2)
        if direction == 'vertical':
            newaxes_figfrac = newaxes_frac_array * height
            for axes_frac_bottom_top in newaxes_figfrac:
                x0new = x0
                widthnew = width
                y0new = y0 + axes_frac_bottom_top[0]
                heightnew = axes_frac_bottom_top[1] - axes_frac_bottom_top[0]
                newpos_list.append([x0new, y0new, widthnew, heightnew])
        return np.array(newpos_list)


def Axes_Replace_Split_Axes(fig, axes_remove, split_fraction, direction='vertical'):
    """
        Add sevearl new axes in the place of the axes_remove(the axes which will be replaced by new ones)

    Parameters
    ----------
    axes_remove
        the axes who is to be removed and whose place will be used to build new axes
    split_fraction
        a list : the fraction of the new axes and the space between axes. the lenght should be "number_of_new_axes" + "number_of_new_axes-1". Note the fraction is
        in terms of the original axes but not the figure
    direction
        'vertical' for vertical split, 'horizontal' for horizontal split.

    Returns
    -------
    out :
        a list of axes, in order of bottom to top.
    """
    newaxes_list = []
    posor = axes_remove.get_position()
    height_or = posor.y1 - posor.y0
    width_or = posor.x1 - posor.x0
    print posor
    print width_or, height_or
    newaxes_pos = Calc_Newaxes_Fraction([posor.x0, posor.y0, width_or, height_or], split_fraction, direction=direction)
    print newaxes_pos
    fig.delaxes(axes_remove)
    for x0new,y0new,widthnew,heightnew in newaxes_pos:
        newaxes_list.append(fig.add_axes([x0new, y0new, widthnew, heightnew]))
    return newaxes_list

def Axes_Set_Breakaxis(bottom_ax,top_ax,h,v):
    # hide the spines between ax and bottom_ax
    top_ax.spines['bottom'].set_visible(False)
    bottom_ax.spines['top'].set_visible(False)
    top_ax.xaxis.tick_top()
    top_ax.tick_params(labeltop='off') # don't put tick labels at the top
    bottom_ax.xaxis.tick_bottom()

    fig = top_ax.figure

    def get_axes_height(axes):
        pos = axes.get_position()
        return pos.bounds[3]

    bottom_axheight = get_axes_height(bottom_ax)
    top_axheight = get_axes_height(top_ax)

    # arguments to pass plot, just so we don't keep repeating them
    h1=h; v1=v*bottom_axheight/top_axheight
    kwargs = dict(transform=top_ax.transAxes, color='k', clip_on=False)
    top_ax.plot((-h1,+h1),(-v1,+v1), **kwargs) # top-left diagonal
    top_ax.plot((1-h1,1+h1),(-v1,+v1), **kwargs) # top-right diagonal

    h2=h; v2=v
    kwargs.update(transform=bottom_ax.transAxes) # switch to the bottom axes
    bottom_ax.plot((-h2,+h2),(1-v2,1+v2), **kwargs) # bottom-left diagonal
    bottom_ax.plot((1-h2,1+h2),(1-v2,1+v2), **kwargs) # bottom-right diagonal

def set_breakaxis2(top_ax,bottom_ax,h,v):
    # hide the spines between ax and bottom_ax
    top_ax.spines['bottom'].set_visible(False)
    bottom_ax.spines['top'].set_visible(False)
    top_ax.xaxis.tick_top()
    top_ax.tick_params(labeltop='off') # don't put tick labels at the top
    bottom_ax.xaxis.tick_bottom()

    fig = top_ax.figure
    topaxpos = top_ax.get_position()
    print "topaxpos",topaxpos
    x0 = topaxpos.x0
    y0 = topaxpos.y0
    x1 = topaxpos.x1
    y1 = topaxpos.y1
    l1 = mat.lines.Line2D([x0-h, y0-v], [x0+h, y0+v],transform=fig.transFigure, figure=fig, color='k')
    l2 = mat.lines.Line2D([x1-h, y0-v], [x1+h, y0+v],transform=fig.transFigure, figure=fig, color='k')

    bottomaxpos = bottom_ax.get_position()
    print "bottomaxpos",bottomaxpos
    x0 = bottomaxpos.x0
    y0 = bottomaxpos.y0
    x1 = bottomaxpos.x1
    y1 = bottomaxpos.y1
    l3 = mat.lines.Line2D([x0-h, y1-v], [x0+h, y1+v],transform=fig.transFigure, figure=fig, color='k')
    l4 = mat.lines.Line2D([x1-h, y1-v], [x1+h, y1+v],transform=fig.transFigure, figure=fig, color='k')

    fig.lines.extend([l1,l2,l3,l4])

def Get_Axes_Ratio(ax):
    """
    Return the height/width ratio of an axes
    """
    box=ax.get_position()
    return (box.y1-box.y0)/(box.x1-box.x0)

def Create_1Axes():
    """
    Purpose: Creat a figure and an axes instance for ploting
    Definition: Create_Axes()
    Example: fig, ax=Create_Axes()
    """
    fig=plt.figure()
    ax=fig.add_subplot(111)
    return fig, ax

def Add_Rectangle(ax,data=None,**kwargs):
    """
    Add Rectangle to axes by specifying data limits.
    data=((x0,y0),(x1,y1))
    """
    height=data[1][1]-data[0][1]
    width=data[1][0]-data[0][0]
    rec=plt.Rectangle((data[0][0],data[0][1]),width,height,**kwargs)
    ax.add_patch(rec)


def Set_AxText_eg():
    fig,ax=Create_1Axes()
    for p in ['lc','uc','ll','ul','lr','ur']:
        Set_AxText(ax,p,pos=p)

def Create_4Axes():
    """
    Purpose: Creat a figure and an axes instance for ploting
    Example: fig, ax=Create_4Axes()
        the sequence:
            ax1  ax2
            ax3  ax4
    """
    fig=plt.figure()
    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3=fig.add_subplot(223)
    ax4=fig.add_subplot(224)
    ax=tuple([ax1,ax2,ax3,ax4])
    return fig, ax


def Create_4ShrinkAxes(pos,frac):
    """
    Purpose: Creat a figure and 4 axes with automatic shrinking
    Example: fig, ax=Create_4ShrinkAxes('w',0.8)
        the sequence:
            ax1  ax2
            ax3  ax4
    """
    fig=plt.figure()
    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3=fig.add_subplot(223)
    ax4=fig.add_subplot(224)

    if pos=='h':
        Set_Axes_Position(ax3,pos,frac)
        Set_Axes_Position(ax4,pos,frac)
        for ax in [ax1,ax2]:
            box=ax.get_position()
            ax.set_position([box.x0, box.y0-box.height*(1-frac), box.width, box.height*frac])
    elif pos=='w':
        Set_Axes_Position(ax1,pos,frac)
        Set_Axes_Position(ax3,pos,frac)
        for ax in [ax2,ax4]:
            box=ax.get_position()
            ax.set_position([box.x0-box.width*(1-frac), box.y0, box.width*frac, box.height])
    else:
        print "use 'w' for width adjustment and 'h' for height adjustment"

    ax=tuple([ax1,ax2,ax3,ax4])
    return fig, ax

def Create_2HAxes():
    """
    Purpose: Creat a figure and 2 horizontal axes instance for ploting
    """
    fig=plt.figure()
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)
    ax=tuple([ax1,ax2])
    return fig, ax

def Create_2VAxes():
    """
    Purpose: Creat a figure and 2 vertical axes instance for ploting
    """
    fig=plt.figure()
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)
    ax=tuple([ax1,ax2])
    return fig, ax

def Create_AxesMatrix(row_num,column_num):
    fig=plt.figure()
    ax=[]
    for i in range(1,row_num*column_num+1):
        axt=fig.add_subplot(row_num,column_num,i)
        ax.append(axt)
    return fig,tuple(ax)


def Create_Twin_Axes():
    from mpl_toolkits.axes_grid1 import host_subplot
    fig=plt.figure()
    host = host_subplot(111)
    par = host.twinx()
    return fig,host,par

def Tight_Legend(ax,loc='upper right',bbox_to_anchor=None,ncolumn=1,fontsize=12,numpoints=1):
    """
    Definition:
        ax.legend(loc=loc,borderpad=0.3,labelspacing=0.2,handletextpad=0.2,handlelength=2,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})
    Note: This is used for only one fig with only axes
    """
    if bbox_to_anchor==None: 
        ax.legend(loc=loc,borderpad=0.1,labelspacing=0.1,handletextpad=0.1,handlelength=2,borderaxespad=0.1,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})
    else:
        ax.legend(loc=loc,bbox_to_anchor=bbox_to_anchor,borderpad=0.1,labelspacing=0.1,handletextpad=0.1,handlelength=2,borderaxespad=0.1,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})

def Tight_Legend_ex(ax,handles,labels,loc='upper right',**kwargs):
    """
    Definition:
        ax.legend(loc=loc,borderpad=0.3,labelspacing=0.2,handletextpad=0.2,handlelength=2,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})
    Note: This is used for only one fig with only axes
    """
    ax.legend(handles,labels,loc=loc,borderpad=0.1,labelspacing=0.1,handletextpad=0.1,handlelength=2,borderaxespad=0.1,columnspacing=0.4,**kwargs)

def Tight_Legend_1Figtwinx(fig,loc='upper left',ncolumn=1,fontsize=12,numpoints=1):
    """
    Definition:
        plt.legend(tuple(pline),tuple(plab),loc=loc,bbox_to_anchor=legend_position,bbox_transform=fig.transFigure,borderpad=0.3,labelspacing=0.2,handletextpad=0.2,handlelength=2,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})

    Note: 1. It's especially used for set legend for 1 figure but with twin axes.
          2. when plot line/points, use plot(...,label='') or plot(...) freely (explicitely seting label or not is both OK), the command can extract label automatically
          3. the loc argument is designed exactly as normal loc in legend. by default legend is in 'upper left' of the overlapped axes.

    """
    ax=fig.get_axes() 
    pline=[]
    plab=[]
    for axt in ax:
        box=axt.get_position()
        dul=tuple([box.x0,box.y1])
        dur=tuple([box.x1,box.y1])
        dll=tuple([box.x0,box.y0])
        dlr=tuple([box.x1,box.y0])
        lines=axt.get_lines()    #from here get the lines with label
        for line in lines:
            if line.get_label()=='_line1':
                pass
            else:
                pline.append(line)
                plab.append(line.get_label())

    if loc=='upper left':
        legend_position=dul
    elif loc=='upper right':
        legend_position=dur
    elif loc=='lower left':
        legend_position=dll
    elif loc=='lower right':
        legend_position=dlr

    plt.legend(tuple(pline),tuple(plab),loc=loc,bbox_to_anchor=legend_position,bbox_transform=fig.transFigure,borderpad=0.3,labelspacing=0.2,handletextpad=0.2,handlelength=2,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})



def Tight_Legend_4ShrinkAxesW(fig,ncolumn=1,fontsize=12,numpoints=1):
    """
    Purpose: To put a legend on the right part of the figure. To be used when one figure has 4 subplots and subplots are shrinked in width.
    Note: The advantage of this method is that you can specify the line as labeled or not as you want (it extracts lables automatically). see the example below.
    Example:
        >>> fig,ax=g.Create_4ShrinkAxes('w',0.7)
        >>> for i,axt in enumerate(ax):
                axt.plot(np.arange(10),color=g.pcolor[i],label='lab'+str(i))
                axt.plot(np.arange(10)+i+1,color=g.pcolor[i+4])
        >>> g.Tight_Legend_4ShrinkAxesW(fig)
    """
    ax=fig.get_axes()
    d=[]
    pline=[]
    plab=[]
    for axt in ax:
        box=axt.get_position()
        d.append([box.x1+0.02,box.y0])   #plus 0.02 to move the legend a litter right
        lines=axt.get_lines()    #from here get the lines with label
        for line in lines:
            if '_line' in line.get_label():
                pass
            else:
                pline.append(line)

    line_col_dic=dict((line.get_color(),line) for line in pline)  #extract lines with unique color
    plab=[line.get_label() for line in line_col_dic.values()]
    if len(plab)==0:
        raise ValueError('Cannot find any line with label, please check if any label has been set.')
    else:
        fig.legend(tuple(line_col_dic.values()),tuple(plab),loc=(max(d)[0],min(d)[1]),borderpad=0.3,labelspacing=0.2,handletextpad=0.2,handlelength=2,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})
    
    #plt.legend(loc=loc,bbox_to_anchor = tuple(max(d)),bbox_transform=fig.transFigure,borderpad=0.3,labelspacing=0.2,handletextpad=0.2,handlelength=2,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})


def Tight_Legend_4ShrinkAxesH(fig,ncolumn=4,fontsize=12,numpoints=1):
    """
    Purpose: To put a legend on the right part of the figure. To be used when one figure has 4 subplots and subplots are shrinked in height.
    Note: The advantage of this method is that you can specify the line as labeled or not as you want (it extracts lables automatically). see the example below.
    Example:
        >>> fig,ax=g.Create_4ShrinkAxes('h',0.7)
        >>> for i,axt in enumerate(ax):
                axt.plot(np.arange(10),color=g.pcolor[i],label='lab'+str(i))
                axt.plot(np.arange(10)+i+1,color=g.pcolor[i+4])
        >>> g.Tight_Legend_4ShrinkAxesH(fig)
    """
    ax=fig.get_axes() 
    d=[]
    pline=[]
    plab=[]
    for axt in ax:
        box=axt.get_position()
        d.append([box.x0,box.y1])
        lines=axt.get_lines()    #from here get the lines with label
        for line in lines:
            if '_line' in line.get_label():
                pass
            else:
                pline.append(line)

    line_col_dic=dict((line.get_color(),line) for line in pline)
    plab=[line.get_label() for line in line_col_dic.values()]
    if len(plab)==0:
        raise ValueError('Cannot find any line with label, please check if any label has been set.')
    else:
        fig.legend(tuple(line_col_dic.values()),tuple(plab),loc=(min(d)[0],max(d)[1]+0.02),borderpad=0.3,labelspacing=0.2,handletextpad=0.2,handlelength=2,columnspacing=0.4,ncol=ncolumn,prop={'size':fontsize})

def Set_Fontsize_Figure(fig,fontsize):
    def match(artist):
        return artist.__module__ == "matplotlib.text"
    for textobj in fig.findobj(match=match):
        textobj.set_fontsize(fontsize)

def Plot_1Horizontal_Zero_Line():
    axc=plt.gca()
    axc.plot(np.arange(axc.get_xlim()[1]),np.zeros(axc.get_xlim()[1]),'r--')

def Plot_1Horizontal_Line(pos=0,lp='r--'):
    """
    pos is the position where needs a vertical line; "lp" receives the linestype string as same in plt.plot() function.
    """
    axc=plt.gca()
    axc.plot(list(axc.get_xlim()),[pos,pos],lp)

def Plot_1Vertical_Line(pos=0,lp='k--',**kwargs):
    """
    pos is the position where needs a vertical line
    """
    axc=plt.gca()
    axc.set_autoscale_on(False)
    axc.plot([pos,pos],list(axc.get_ylim()),lp,**kwargs)
    axc.set_autoscale_on(True)

def Plot_Vertical_Lines(ax,pos=[0],lp='k--',**kwargs):
    """
    pos is the position where needs a vertical line
    """
    lines=[]
    ax.set_autoscale_on(False)
    for pos1 in pos:
        lines.append(ax.plot([pos1,pos1],list(ax.get_ylim()),lp,**kwargs))
    ax.set_autoscale_on(True)
    return lines

def plot_mean_std(data,move_ave=None):
    """
    Purpose: plot for a 2D array the average with std
    Definition: plt.fill_between(np.arange(len(dmean)),dmean-dstd,dmean+dstd,color=c2['lb'])
    Note:
        1. the rows are variation of data and column number as length of data
        2. the move_ave only applies for non masked array. For masked array, usually movign average is not used so move_ave can be set as None to avoid moving average.
        3. moving average is done with mathex.move_ave

    """
    if data.ndim !=2:
        print 'please provide 2D array!'
    if move_ave!=None:
        data=mathex.move_ave2d(data,move_ave)
    dmean=np.ma.mean(data,axis=0)
    dstd=np.ma.std(data,axis=0)
    plt.plot(dmean,color=c2['db'])
    plt.fill_between(np.arange(len(dmean)),dmean-dstd,dmean+dstd,color=c2['lb'])

def plot_mean_std_ax(ax,data,move_ave=None,lab='label',colmean=c2['db'],colstd=c2['lb'],alph=0.4):
    """
    Purpose: plot for a 2D array the average with std
    Definition:
    Note:
        1. the rows are variation of data and column number as length of data
        2. the move_ave only applies for non masked array. For masked array, usually movign average is not used so move_ave can be set as None to avoid moving average.
        3. moving average is done with mathex.move_ave

    """
    if data.ndim !=2:
        print 'please provide 2D array!'
    if move_ave!=None:
        data=mathex.move_ave2d(data,move_ave)
    dmean=np.ma.average(data,axis=0)
    dstd=np.ma.std(data,axis=0)
    line=ax.plot(dmean,color=colmean,label=lab)
    poly=ax.fill_between(np.arange(len(dmean)),dmean-dstd,dmean+dstd,color=colstd,alpha=0.4)
    return line,poly

def pfread(filename):
    fob=file(filename,'r')
    data=pk.load(fob)
    fob.close()
    return data

def pfdump(data,filename):
    """
    Purpose: use python pickle to dump data object to a file
    Definition: pfdump(data,filename)
    """
    fob=file(filename,'w')
    data=pk.dump(data,fob)
    fob.close()

def ipsearch(string):
    """
    [str(l) for l in  _ih if l.startswith(string)]
    """
    [str(l) for l in  _ih if l.startswith(string)]

class cm(object):
    #red2greendict={'red': ((0.0, 0.0, 1.0),
    #                 (0.5, 0.3, 0.0),
    #                 (1.0, 1.0, 1.0)),
    #         'green': ((0.0, 0.0, 0.0),
    #                   (0.5, 0.0, 0.3),
    #                   (1.0, 1.0, 1.0)),
    #         'blue': ((0.0, 0.0, 0.0),
    #                  (0.5, 0.0, 0.0),
    #                  (1.0, 1.0, 0.0))}
    _red2green_rgblist=[(242,6,6),(246,248,159),(207,249,144),(66,128,25)]
    _red2greendict=rgb2cmdic(_red2green_rgblist,[0,0.5,0.5,1])
    red2green = mat.colors.LinearSegmentedColormap('red2green',_red2greendict,256)
    #
    red2bluedict= {'blue':[(0.0, 0.0, 0.0), 
                           (0.5, 0.0, 0.0),
                           (0.5, 0.5, 0.5),
                           (1.0, 1.0, 1.0)],
                 'green': [(0.0, 0.0, 0.0),
                           (0.5, 1.0, 1.0),
                           (0.5, 0.0, 0.0),
                           (1.0, 1.0, 1.0)],
                   'red': [(0.0, 1.0, 1.0),
                           (0.5, 1.0, 1.0),
                           (0.5, 0.0, 0.0),
                           (1.0, 0.0, 0.0)]}
    red2blue = mat.colors.LinearSegmentedColormap('red2blue',red2bluedict,256)
    
    red2blue_rgblist=[(255,51,0),(255,255,153),(175,238,238),(0,0,250)]
    red2bluedict2=rgb2cmdic([(255,51,0),(255,255,153),(175,238,238),(0,0,250)],[0,0.5,0.5,1])
    red2blue2 = mat.colors.LinearSegmentedColormap('red2blue2',red2bluedict2,256)
    red2bluedict2_r=rgb2cmdic(red2blue_rgblist[::-1],[0,0.5,0.5,1])
    red2blue2_r = mat.colors.LinearSegmentedColormap('red2blue2_r',red2bluedict2_r,256)


    ###
    cdict = {'red':   [(0.0,  0.0, 0.0),
                       (0.5,  1.0, 1.0),
                       (1.0,  1.0, 1.0)],
             'green': [(0.0,  0.0, 0.0),
                       (0.25, 0.0, 0.0),
                       (0.75, 1.0, 1.0),
                       (1.0,  1.0, 1.0)],
             'blue':  [(0.0,  0.0, 0.0),
                       (0.5,  0.0, 0.0),
                       (1.0,  1.0, 1.0)]}
    cm1 = mat.colors.LinearSegmentedColormap('cm1',cdict,256)
    ###
    clist=[(112,0,0),(148,0,148),(255,0,102),(245,245,0),(107,214,0),(0,43,87),(0,51,204)]
    PurpleGreenBlueDic=rgb2cmdic([(112,0,0),(148,0,148),(255,0,102),(245,245,0),(107,214,0),(0,43,87),(0,51,204)],[0,0.1,0.2,0.49,0.5,0.9,1])
    PurpleGreenBlue = mat.colors.LinearSegmentedColormap('PurpleGreenBlue',PurpleGreenBlueDic,256)
    PurpleGreenBlueDic_rev=rgb2cmdic(clist[::-1],[0,0.1,0.2,0.49,0.5,0.9,1])
    PurpleGreenBlue_rev = mat.colors.LinearSegmentedColormap('PurpleGreenBlue_rev',PurpleGreenBlueDic_rev,256)
    ### this bar is from Dai et al. 2011
    clist=[(0.062745098039215685, 0.3843137254901961, 0.68235294117647061),
     (0.0039215686274509803, 0.58431372549019611, 0.83529411764705885),
     (0.0, 0.70980392156862748, 0.88627450980392153),
     (0.0, 0.70588235294117652, 0.80000000000000004),
     (0.0, 0.70588235294117652, 0.50980392156862742),
     (0.38039215686274508, 0.74509803921568629, 0.29411764705882354),
     (0.94509803921568625, 0.92941176470588238, 0.40392156862745099),
     (1.0, 0.77647058823529413, 0.26666666666666666),
     (0.97647058823529409, 0.59999999999999998, 0.30196078431372547),
     (0.94901960784313721, 0.20784313725490197, 0.43529411764705883),
     (0.87058823529411766, 0.50980392156862742, 0.70980392156862748),
     (0.55686274509803924, 0.38823529411764707, 0.6705882352941176)]
    precip=mat.colors.ListedColormap(clist[::-1],'precip')
    precip_rev=mat.colors.ListedColormap(clist,'precip_rev')

    ##temp; the same sequence of color with precip_rev by the level can be maniplated with templevel.
    ##the twelve colors; yellow: 7; green:6. 
    templevel=[0.0, 0.09, 0.18, 0.27, 0.36, 0.53, 0.58, 0.64, 0.73, 0.82, 0.91, 1.0]
    tempcolor=clist[:]
    tempcolor=np.array(tempcolor)*255.
    tempc=[tuple(i) for i in tempcolor]
    tempcmdic=rgb2cmdic(tempc,templevel)
    tempcm=mat.colors.LinearSegmentedColormap('tempcm',tempcmdic,256)


class mapb(object):
    caal=(40,75,-170,-50)

def plot_OLS_reg(ax,x,y,c='k',ls='--'):
    """
    Purpose: plot OLS regression line for y~x on axes ax.
    Note:
        1. being able to automatically mask np.nan
    Return:
        return line[0],[slope, intercept, r_value, p_value, std_err] 

    """
    xnew=pb.MaskArrayByNan(x)
    ynew=pb.MaskArrayByNan(y)
    [slope, intercept, r_value, p_value, std_err] = sp.stats.mstats.linregress(xnew,ynew)
    xnew_plot=pb.linspace_array(xnew)
    line=ax.plot(xnew_plot,xnew_plot*slope+intercept,color=c,linestyle=ls)
    return line[0],[slope, intercept, r_value, p_value, std_err] 
        
def plot_RTO_OLSreg2var(ax,x,y,c='k',ls='--'):
    """
    Purpose: plot OLS regression line for y~x through origin on axes ax. RTO is short for "Regression Through Origin"
    Note:
        1. being able to automatically mask np.nan  
        2. only 2 variables are tested.
    Return:
        return line[0],[slope, r_value, p_value, std_err] 
    """
    #extract data that are not NaN or masked
    xnew=pb.MaskArrayByNan(x)
    ynew=pb.MaskArrayByNan(y)
    x_reg,y_reg=pb.shared_unmask_data(xnew,ynew)
    #use R to do regression
    rpy.r.assign('x_reg',x_reg)
    rpy.r.assign('y_reg',y_reg)
    lmRTO=rpy.r('lmRTO=lm(y_reg ~ x_reg-1)')  #Regression through origin
    summary_lmRTO=rpy.r('summary_lmRTO=summary(lmRTO)')
      #note the 1st, 2nd, 3rd, 4th elements for summary_lmRTO['coefficients'] is Estimate; Std. Error; t value; Pr(>|t|)
    slope=summary_lmRTO['coefficients'][0][0]
    p_value=summary_lmRTO['coefficients'][0][3]
    std_err=summary_lmRTO['coefficients'][0][1]
    r_value=np.sqrt(summary_lmRTO['r.squared'])
    #creat a new array for ploting regression line
    xnew_plot=pb.linspace_array(xnew)
    line=ax.plot(xnew_plot,xnew_plot*slope,color=c,linestyle=ls)
    return line[0],[slope, r_value, p_value, std_err] 


def legpro_points(colorlist,labellist):
    """
    Purpose: Creat a group of points/lines line2D used as proxy for legend, return tuple of (linelist,labellist).
    Use:
        proleg=g.legpro_points(['r','g'],['red','green'])
        ax.legend(proleg[0],proleg[1],**kwarg)
    """
    point_list=[mat.lines.Line2D([],[],marker='o',ms=5,mfc=c,mew=0,ls='none',color=c) for c in colorlist]
    return (point_list,labellist)

def legpro_lines(colorlist,labellist,ls='-',**kwargs):
    """
    Purpose: Creat a group of points/lines line2D used as proxy for legend, return tuple of (linelist,labellist).
    Use:
        proleg=g.legpro_points(['r','g'],['red','green'])
        ax.legend(proleg[0],proleg[1],**kwargs)
    """
    point_list=[mat.lines.Line2D([],[],color=c,ls=ls,**kwargs) for c in colorlist]
    return (point_list,labellist)

def Set_Cbar_Label_Parallel(cbar,label_list,leftpos=1.2,**kwargs):
    """
    This is to set colorbar label besie the colorbar. see the example at:/homel/ychao/python/script/set_label_parallel_colorbar.py
    """
    cbar.set_ticklabels([])
    cbar.ax.tick_params(right='off',left='off')
    yloc=(cbar.values-cbar.boundaries[0])/(cbar.boundaries[-1]-cbar.boundaries[0])
    if len(label_list) != len(yloc):
        raise ValueError("the lenght of cbar segments and label list are not equal!")
    else:
        for label,ypos in zip(label_list,yloc):
            cbar.ax.text(leftpos,ypos,label,ha='left',va='center',**kwargs)

def setp(*artist,**kwargs):
    """
    Purpose: set artist properties by kwagrs pairs in an easy and flexible way.
    Note:
        1.  artist will be flat using pb.iteflat so you can use mixed types of matplotlib artists as long as they have the same keyword properties.
        2.  when artist is a tuple or list,kwargs[key] can also be set as tuple or list, but when kwargs[key] is only one value, it will be broadcast 
            to the same length with artist automatically.
    """
    if len(artist)==1 and isinstance(artist,(tuple,list)):
        artist_list=pb.iteflat(artist[0])
    else:
        artist_list=pb.iteflat(artist)

    for key in kwargs:
        value=kwargs[key]
        if not isinstance(value,Iterable) or isinstance(value,str):
            value_list=[value]*len(artist_list)
        else:
            if len(value)==1:
                value_list=value*len(artist_list)
            else:
                value_list=pb.iteflat(value)
        if len(artist_list)!=len(value_list):
            raise ValueError('artist list lenght {0} is not equal to value list length {1}'.format(len(artist_list),len(value_list)))
        else:
            for art,val in zip(artist_list,value_list):
                plt.setp(art,key,val)
        print key,value_list,'has been set'
    return artist_list,[key]*len(artist_list),value_list

class ProxyLegend(object):
    """
    Tags are used as labels for proxy legend.
    """
    def __init__(self,newdata=None):
        '''
        newdata should be a dict with (tag,handler) pairs.
        '''
        if newdata == None:
            self.data = {}
        elif isinstance(newdata,list):
            self.data = copy.copy(dict.fromkeys(newdata))
        elif isinstance(newdata,dict):
            self.data = copy.copy(newdata)
        else:
            raise TypeError("only list and dict accepted!")

        if self.data == {}:
            self.tags = []
        else:
            self.tags = self.data.keys()

    @classmethod
    def merge_pleg(cls,*pleglist):
        '''
        Notes:
        ------
        merge_pleg will keep the sequence of tags for each pleg and
            and the sequence of pleglist.
        '''
        data = {}
        tags = []
        for pleg in pleglist:
            data.update(pleg.data)
            tags.extend(pleg.tags)
        plegnew = ProxyLegend(data)
        plegnew.set_tag_order(tags)
        return plegnew

    def set_tag_order(self,tagseq=None):
        """
        Set tag order and this order will be kept throughout all the class
        method when default taglist is used.
        """
        if sorted(self.tags) == sorted(tagseq):
            self.tags = tagseq[:]
        else:
            raise ValueError('ordered tag list not equal to present taglist')

    def add_line_by_tag(self,tag,**kwargs):
        line=mat.lines.Line2D([],[],**kwargs)
        self.data[tag]=line
        if tag in self.tags:
            pass
        else:
            self.tags.append(tag)

    def add_rec_by_tag(self,tag,**kwargs):
        """
        Add a mat.patches.Rectangle instance; but notice setting the
            edgecolor as 'none' will not lead to a none edgecolor
            in the legend handle.
        """
        rect=mat.patches.Rectangle((0, 0), 1, 1,**kwargs)
        self.data[tag]=rect
        if tag in self.tags:
            pass
        else:
            self.tags.append(tag)

    def add_scatter_by_tag(self, tag, s=20, c='k', marker='o', cmap=None,
                           norm=None,vmin=None, vmax=None, alpha=None,
                           linewidths=None,verts=None, **kwargs):
        scatter = mat.scatter([],[],s=s, c='k', marker='o', cmap=cmap,
                              norm=norm,vmin=vmin, vmax=vmax, alpha=alpha,
                              linewidths=linewidths, verts=verts, **kwargs)
        self.data[tag] = scatter
        if tag in self.tags:
            pass
        else:
            self.tags.append(tag)

    def add_proxy_by_dic(self,proxy_dic):
        if pb.Lists_Intersection_Check([self.tags,proxy_dic.keys()]):
            raise ValueError("Try to add new proxy to the same tag!")
        else:
            self.data.update(proxy_dic)
            self.tags = self.tags + proxy_dic.keys()

    def copy(self):
        proleg_new = ProxyLegend(self.data)
        return proleg_new

    @staticmethod
    def _check_void_handle_label(label):
        if label == '':
            return True
        elif label[0] == '_':
            return True
        elif label == 'None':
            return True

    def _get_handle_label(self,tagseq=None):
        handle_list=[]
        label_list=[]
        if tagseq==None:
            tagseq=self.tags
        for tag in tagseq:
            handle=self.data[tag]
            handle_list.append(handle)
            label=handle.get_label()
            if ProxyLegend._check_void_handle_label(label):
                label_list.append(tag)
            else:
                label_list.append(label)
        return handle_list,label_list


    def create_legend(self,ax,tagseq=None,**kwargs):
        handle_list,label_list=self._get_handle_label(tagseq)
        return ax.legend(handle_list,label_list,**kwargs)

    def add_lines_by_color(self,colorlist,labellist,**kwargs):
        for lab,c in zip(labellist,colorlist):
            self.add_line_by_tag(lab,color=c,ls='-',**kwargs)

    def add_recs_by_color(self,colorlist,labellist,**kwargs):
        for lab,c in zip(labellist,colorlist):
            self.add_rec_by_tag(lab,fc=c,**kwargs)

    def add_points_by_color(self,colorlist,labellist,**kwargs):
        for lab,c in zip(labellist,colorlist):
            self.add_line_by_tag(lab,color=c,marker='o',ls='none',**kwargs)

    def create_legend_top(self,ax,tagseq=None,expand=True,**kwargs):
        if expand==True:
            mode='expand'
        else:
            mode=None
        self.create_legend(ax,tagseq=tagseq,loc=(0,1.02),mode=mode,borderaxespad=0,**kwargs)

    def add_mixed_handle(self,lines=None,points=None):
        """
        lines=(['r','g'],['red','green'])
        points=(['r','g'],['red','green'])
        lines, points should be 2-length tuple with 1st as colorlist and 2nd as labellist.
        """
        try:
            self.add_lines_by_color(*lines)
        except TypeError:
            pass

        try:
            self.add_points_by_color(*points)
        except TypeError:
            pass


    def add_handle_by_dic(self,**nested_attr_tag_value_dic):
        """
        Purpose: add handles by specifying a attr_name/tag:value dictionary.
        """
        #initiate the dictionary (with tags as keys, and another subdic with attr/attr_value paris as values)
        plot_attr_tag_dic={}
        if self.tags==None:
            raise ValueError("Please provide tags when creating ProxyLeg instance before using this method")
        else:
            for tag in self.tags:
                plot_attr_tag_dic[tag]={}
        #do the dictionary permutation by keys
        for attr_name,tag_attr_value_dic in nested_attr_tag_value_dic.items():
            #This if/else build the final dic to be used.
            if not isinstance(tag_attr_value_dic,dict):
                #assume a list
                if isinstance(tag_attr_value_dic,list):
                    #tag_attr_value_dic is a list of (tag,attr_value) tuples.
                    if isinstance(tag_attr_value_dic[0],tuple):
                        final_dic=pb.Dic_by_List_of_Tuple(tag_attr_value_dic)
                    #a list of only values
                    else:
                        if len(tag_attr_value_dic)!=len(self.tags):
                            raise ValueError('existing tags have len {0} but input list len is {1}'.format(len(self.tags),len(tag_attr_value_dic)))
                        else:
                            final_dic=dict(zip(self.tags,tag_attr_value_dic))
                #assume a single value (number or string)
                else:
                    final_dic=dict(zip(self.tags,len(self.tags)*[tag_attr_value_dic]))
            #in case of dict
            else:
                final_dic=tag_attr_value_dic
            #apply value to tag
            for tag,attr_value in final_dic.items():
                if tag not in self.tags:
                    raise ValueError("tag '{0}' dose not exist in current tags, you may have a wrong tag name".format(tag))
                else:
                    plot_attr_tag_dic[tag][attr_name]=attr_value
        #add handles
        for tag,tag_attr_subdic in plot_attr_tag_dic.items():
            self.add_line_by_tag(tag,**tag_attr_subdic)

