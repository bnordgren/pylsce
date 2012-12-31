import bsite as bsite
import matplotlib
import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import pupynere as pu
import pandas as pa
import pickle as pk
import os as os
import re as re
import scipy as sp
import mpl_toolkits.basemap as bmp
from mpl_toolkits.basemap import cm
import pdb
import netCDF4 as nc
from matplotlib.backends.backend_pdf import PdfPages
import pb
import numpy as np
import weakref
from collections import Iterable
import g
import exio
import copy as copy

def gsetp(*artist,**kwargs):
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

def _replace_none_colorlist(colors=None,num=None):
    if colors == None:
        if num <= len(g.pcolor):
            return g.pcolor[0:num]
        else:
            raise ValueError("g.pcolor is not long enough when using default colorlist")
    else:
        return colors

def _replace_none_axes(ax):
    if ax==None:
        fig, axnew = g.Create_1Axes()
        return axnew
    else:
        return ax

def _replace_none_by_given(orinput,default):
    if orinput == None:
        return default
    else:
        return orinput

def Dic_Subset_End(indic,end_num):
    """
    subset a dictionary by retaining only the last "end_num" elements.
    Note:
        1. Test has been done for only the case that indic[key] is 1D ndarray.
    """
    outdic={}
    for key,value in indic.items():
        try:
            newvalue=value[-end_num:]
        except:
            newvalue=value
        outdic[key]=newvalue
    return outdic

def Dic_Subset_Begin(indic,end_num):
    """
    subset a dictionary by retaining only the beginning "end_num" elements.
    Note:
        1. Test has been done for only the case that indic[key] is 1D ndarray.
    """
    outdic={}
    for key,value in indic.items():
        try:
            newvalue=value[0:end_num]
        except:
            newvalue=value
        outdic[key]=newvalue
    return outdic

def Dic_Extract_By_Subkeylist(indic,keylist):
    """
    Return a new dic by extracting the key/value paris present in keylist
    """
    outdic={}
    for key in keylist:
        try:
            outdic[key]=indic[key]
        except KeyError:
            pass
    return outdic

def Dic_Remove_By_Subkeylist(indic,keylist):
    """
    Return a new dic, with key/value pairs present in keylist removed.
    """
    outdic=indic.copy()
    for key in outdic.keys():
        if key in keylist:
            del outdic[key]
    return outdic

def StringListAnotB(listA,listB):
    return [i for i in listA if i not in listB]

def FilterStringList(keyword,input_list):
    return [x for x in input_list if re.search(keyword,x)]

def Is_Nested_Dic(indic):
    for value in indic.values():
        if isinstance(value,dict):
            return True
        else:
            pass
    return False

def iteflat(inlist):
    """
    Purpose: flat all things that are member elements in a iteralbe object (eg. tuple,list) to a list.
            inlist can also be itself a iteralbe with non-iterable elements (eg. np array). 
    Note:
        1. Check the source code when confused by this function and the definition is very clear.
        2. 2nd level nested lists will not be flatted. eg. 
           ['a','b','c'] --> no nested list
           [['a','b','c'],['e','f']]  --> 1st nested list.
           [['a','b','c'],['e','f',['h','l']]]  --> 2nd nested list.
    Example:
        >>> inlist1=[(np.NINF),(3),(np.PINF)]
        >>> inlist2=[(-10,-3),np.arange(-2,2.1,0.5),(3,10)]
        >>> inlist3=[(-10,-3),(0),(3,10)]

        In [28]: pb.iteflat(inlist1)
        Out[28]: [-inf, 3, inf]

        In [29]: pb.iteflat(inlist2)
        Out[29]: [-10, -3, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3, 10]

        In [30]: pb.iteflat(inlist3)
        Out[30]: [-10, -3, 0, 3, 10]

        In [55]: pb.iteflat([(-10,-3),0,(3,10)])
        Out[55]: [-10, -3, 0, 3, 10]

        In [56]: pb.iteflat(np.arange(10))
        Out[56]: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        In [58]: pb.iteflat(((0,2,3),(1,2,3)))
        Out[58]: [0, 2, 3, 1, 2, 3]

        In [122]: pb.iteflat([['a','b','c'],['e','f',['h','l']]])
        Out[122]: ['a', 'b', 'c', 'e', 'f', ['h', 'l']]

        In [123]: pb.iteflat([['a','b','xyz'],['e','f',['h','l']]])
        Out[123]: ['a', 'b', 'xyz', 'e', 'f', ['h', 'l']]

    """
    outlist=[]
    for member in inlist:
        if not isinstance(member,Iterable) or isinstance(member,str):
            outlist.append(member)
        else:
            for submember in member:
                #return iteflat(member)
                outlist.append(submember)
    return outlist

def Dic_Test_Empty_Dic(indic):
    """
    check if indic is a (nested) empty dictionary.
    """
    if indic.keys()!=[]:
        for key,keydata in indic.items():
            if isinstance(keydata,dict):
                if keydata=={}:
                    pass
                else:
                    return Dic_Test_Empty_Dic(keydata)
            else:
                return False
        return True
    else:
        return False

def _get_attr_value_from_objins_dic(objins_dic,*attr_list):
    """
    get attribute value list from a group of object instances by specifying attribute names.
    Return a dict with attr name as key and attr_value list as as key values.
    objins: object instance
    """
    attr_dic={}
    for attr_name in attr_list:
        attr_dic[attr_name]=dict()
        for tag,objins in objins_dic.items():
            attr_dic[attr_name][tag]=objins.__getattribute__(attr_name)
    return attr_dic

class Pdata(object):
    """
    There are two ways to set customized keyword, it can be done using add_attr_by_tag function or passed when calling ploting function. 
    Yet the seting by add_attr_by_tag has high priority.
    Pdata is an object constructed for easy ploting in matplotlib.
    """
    #declare base keylist
    _data_base_keylist=['x','y','xerrl','xerrh','yerrl','yerrh']
    _extra_base_keylist=['label','bwidth','bbottom','bleftshift']
    _extra_nonplot_keylist=['unit','ylab']
    _all_nonplot_keylist = _extra_nonplot_keylist + _data_base_keylist + _extra_base_keylist
    _extra_base_keylist_default=dict(zip(_extra_base_keylist,[None,0.5,0,0]))

    _new_entry=dict(x=None,y=None,xerrl=None,xerrh=None,yerrl=None,yerrh=None)
    _new_entry.update(_extra_base_keylist_default)

    _scatter_attr_base_keylist=['ssize','scolor','smarker','scmap','snorm','svmin','svmax',]
    _error_attr_base_keylist=['efmt','ecolor','elinewidth','capsize','barsabove','lolims','uplims','xlolims','xuplims']
    _bar_attr_base_keylist=['bcolor']
    _plot_attr_dic=dict(scatter=_scatter_attr_base_keylist,errorbar=_error_attr_base_keylist,bar=_bar_attr_base_keylist)
    _plot_attr_keylist_all = _scatter_attr_base_keylist + _error_attr_base_keylist + _bar_attr_base_keylist

    #initiate default plot attritubtes
    _error_attr_default={             \
                 'efmt':None,\
                 'ecolor':'r',   \
                 'elinewidth':None, \
                 'capsize':3, \
                 'barsabove':False, \
                 'lolims':False,\
                 'uplims':False,\
                 'xlolims':False,\
                 'xuplims':False \
                 }
    _scatter_attr_default={             \
                 'ssize':20,   \
                 'scolor':'b', \
                 'smarker':'o',\
                 'scmap':None, \
                 'snorm':None, \
                 'svmin':None, \
                 'svmax':None, \
                 }
    _bar_attr_default={        \
                'bcolor':'b'
                }

    #set default extra plot attribute keyword argument values
    _plot_attr_default={}
    _plot_attr_default.update(_error_attr_default)
    _plot_attr_default.update(_scatter_attr_default)
    _plot_attr_default.update(_bar_attr_default)


    #not used
    def _get_plot_attr_list_except(self,plottype):
        return StringListAnotB(Pdata._plot_attr_keylist_all, Pdata._plot_attr_dic[plottype])

    #tested
    def _plot_attr_keylist_check(self):
        inlist = Pdata._plot_attr_keylist_all
        if any(inlist.count(x) > 1 for x in inlist):
            raise ValueError('plot base attribute lists have duplicates!')



    def __init__(self,newdata={}):
        self._plot_attr_keylist_check()
        self.data=newdata.copy()
        if self.data == {}:
            self._taglist = []
        else:
            self._taglist = self.data.keys()
            self._data_complete_check_all()

    def add_tag(self,tag=None):
        self.data[tag]=Pdata._new_entry.copy()
        self._taglist.append(tag)

    def _add_indata_by_tag_and_column(self,indata,tag,column):
        self.data[tag][column]=indata
    def addx(self,indata,tag):
        self._add_indata_by_tag_and_column(indata,tag,'x')
    def addy(self,indata,tag):
        self._add_indata_by_tag_and_column(indata,tag,'y')
    def _broadcast_scalar_byx(self,tag,scalar):
        return np.repeat(scalar,len(self.data[tag]['x']))
    def addxerrl(self,indata,tag):
        if not isinstance(indata,Iterable):
            indata=self._broadcast_scalar_byx(tag,indata)
        self._add_indata_by_tag_and_column(indata,tag,'xerrl')
    def addxerrh(self,indata,tag):
        if not isinstance(indata,Iterable):
            indata=self._broadcast_scalar_byx(tag,indata)
        self._add_indata_by_tag_and_column(indata,tag,'xerrh')
    def addyerrl(self,indata,tag):
        if not isinstance(indata,Iterable):
            indata=self._broadcast_scalar_byx(tag,indata)
        self._add_indata_by_tag_and_column(indata,tag,'yerrl')
    def addyerrh(self,indata,tag):
        if not isinstance(indata,Iterable):
            indata=self._broadcast_scalar_byx(tag,indata)
        self._add_indata_by_tag_and_column(indata,tag,'yerrh')

    #TESTED
    def add_entry_by_dic(self,**kwargs):
        """
        Add a complete entry by dictionary. note that kwargs[key] must
            have keys()=['x','y',...]
        Parameters:
        ----------
        kwargs: kwargs are tag/tag_value pairs, tag_value is again a dict
            with 'x'/'y'/'yerrh' ... etc as keys.

        """
        for tag,tag_value in kwargs.items():
            self.add_tag(tag)
            for attr_name,attr_value in tag_value.items():
                self.data[tag][attr_name] = attr_value

    def add_entry_noerror(self,x=None,y=None,tag=None):
        self.add_tag(tag)
        if x == None:
            self.addx(np.arange(len(y))+1,tag)
        else:
            if len(y)!=len(x):
                raise ValueError('''lenght of ydata for 'tag' {0} is {1},
                    not equal to length of xdata with length of {2}'''
                    .format(tag,len(y),len(x)))
            else:
                self.addx(x,tag)
        self.addy(y,tag)

    def add_entry_singleYerror(self,x,y,yerr,tag):
        self.add_tag(tag)
        self.addx(x,tag)
        self.addy(y,tag)
        self.addyerrl(yerr,tag)

    def add_entry_singleYerror3(self,data_array,tag):
        """
        add an entry by giving one single 3Xn np.ndarray, with 1st row as X,
            2nd row as Y, 3rd row as Yerr.
        """
        x=data_array[0];y=data_array[1];yerr=data_array[2]
        self.add_tag(tag)
        self.addx(x,tag)
        self.addy(y,tag)
        self.addyerrl(yerr,tag)

    def add_entry_doubleYerror(self,x,y,yerrl,yerrh,tag):
        self.add_tag(tag)
        self.addx(x,tag)
        self.addy(y,tag)
        self.addyerrl(yerrl,tag)
        self.addyerrh(yerrh,tag)

    def add_entry_sharex_noerror_by_dic(self,ydic,x=None):
        """
        Add several tags which share the same x data; ydic will be a dictionary
            (or pandas DataFrame) with tag:ydata pairs. the length of x must
            be equal to all that of ydic[tag]
        """
        if isinstance(ydic,pa.DataFrame):
            newydic={}
            for key in ydic.columns:
                newydic[key]=np.array(ydic[key])
            ydic=newydic

        for tag,ydata in ydic.items():
            self.add_entry_noerror(x=x,y=ydata,tag=tag)

    def add_entry_noerror_by_dic_default_xindex(self,ydic):
        print "Warning! this will be deprecated!"
        for tag,ydata in ydic.items():
            x=np.arange(len(ydata))+1
            self.add_entry_noerror(x=x,y=ydata,tag=tag)

    #TODO: no distinction has been made for keys 'x' and 'scolor'
    def subset_end(self,end_num):
        """
        Subset the pdata by retaining only the last 'end_num' number of
            last elements.
        """
        pdata=self.copy()
        pdatanew=Pdata()
        for tag in pdata.data.keys():
            pdatanew.data[tag]=Dic_Subset_End(pdata.data[tag],end_num)
        return pdatanew

    #TODO: same problem as subset_begin
    def subset_begin(self,end_num):
        """
        Subset the pdata by retaining only the beginning 'end_num' number of
            last elements.
        """
        pdata=self.copy()
        pdatanew=Pdata()
        for tag in pdata.data.keys():
            pdatanew.data[tag]=Dic_Subset_Begin(pdata.data[tag],end_num)
        return pdatanew

    def copy(self):
        data=copy.deepcopy(self.data)
        pdata=Pdata(data)
        return pdata

    #NON_TESTED
    def add_entry_df_groupby_column(self,indf,tag_col=None,**kwargs):
        """
        Purpose: Add tag:tag_dic pairs by grouping a pandas DataFrame by a
            column. The unique values in the column which is used as
            groupby will serve as tags. x,y,xerrl... etc. are specified
            by using kwargs. supported keywords: ['x','y','xerrl','xerrh',
            'yerrl','yerrh']
        Example:
            >>> comgpp_fluxdailyobe=pa.read_csv('/home/chaoyue/python/
                testdata/comgpp_fluxdailyobe_data.csv')
            >>> pdata=Pdata.Pdata()
            >>> pdata.add_entry_df_groupby_column(comgpp_fluxdailyobe,
                tag_col='site',x='mod',y='GEP_Amiro')
            >>> pdata.add_attr_by_tag(scolor=['g','b','r'])
            >>> fig,ax=g.Create_1Axes()
            >>> pdata.scatter(ax)
            >>> pdata.set_legend_scatter(ax,taglab=True)
        """
        indf_groupby=indf.groupby(tag_col)
        for tag,tag_df in indf_groupby:
            entry_dic={}
            for basic_key,df_colname in kwargs.items():
                entry_dic[basic_key]=np.array(tag_df[df_colname])
            self.add_entry_by_dic(**{tag:entry_dic})


    def list_tags(self,tagkw=None):
        """
        Method to formally retrieve Pdata._taglist
        """
        if tagkw==None:
            return self._taglist
        else:
            return FilterStringList(tagkw, self._taglist)

    def list_keys_for_tag(self,tag):
        """
        list the keys for the specified tag.
        """
        return self.data[tag].keys()


    def set_tag_order(self,tagseq=None):
        """
        Set tag order and this order will be kept throughout all the class
        method when default taglist is used.
        """
        if sorted(self._taglist) == sorted(tagseq):
            self._taglist = tagseq
        else:
            raise ValueError('ordered tag list not equal to present taglist')

    def _set_default_tag(self,taglist='all'):
        if taglist=='all':
            return self._taglist
        else:
            return taglist

    def list_attr(self,*attr_name):
        outdic={}
        for tag,tag_data in self.data.items():
            outdic[tag]=Dic_Extract_By_Subkeylist(tag_data,iteflat([attr_name]))
        return outdic

    def list_attr_extra_base(self):
        attr_extra_dic=self.list_attr(*self._extra_base_keylist)
        if Dic_Test_Empty_Dic(attr_extra_dic):
            return None
        else:
            return attr_extra_dic

    def list_attr_plot(self):
        attr_extra_dic=self.list_attr(*Pdata._plot_attr_keylist_all)
        if Dic_Test_Empty_Dic(attr_extra_dic):
            return None
        else:
            return attr_extra_dic

    def get_data_as_dic(self,attr_name,taglist='all'):
        """
        Get the spedcified x/y/yerr... or other attribute data as a
            dictionary with tags as keys
        """
        taglist = self._set_default_tag(taglist)
        data_dic={}
        for tag in taglist:
            data_dic[tag] = self.data[tag][attr_name]
        return data_dic

    def shift_ydata(self,shift=None):
        """
        Shift the y data by given shift value in a progressive way, this
            is mainly for comparing the data with the smae y value in a
            more sensible way. i.e., to shift the data a little bit for
            aoviding the overlapping of the lines.
        """
        for i,tag in enumerate(self._taglist):
            self.data[tag]['y']=self.data[tag]['y']-shift*i

    def apply_function(self, func=None, axis=None, taglist='all', copy=False):
        """
        Apply a function either 'x' or 'y' axis or 'both' or 'diff', if
            axis=='diff', func should be supplied with a dictionary by
            using ('x'/'y',x_func/y_func) pairs.

        Parameters:
        -----------
        copy: return if copy if copy==True.
        """
        if copy == True:
            pdtemp = self.copy()
        else:
            pdtemp = self

        taglist=pdtemp._set_default_tag(taglist)
        for tag in taglist:
            if axis == 'x':
                pdtemp.data[tag]['x']=func(pdtemp.data[tag]['x'])
            elif axis == 'y':
                pdtemp.data[tag]['y']=func(pdtemp.data[tag]['y'])
            elif axis == 'both':
                pdtemp.data[tag]['x']=func(pdtemp.data[tag]['x'])
                pdtemp.data[tag]['y']=func(pdtemp.data[tag]['y'])
            elif axis == 'diff':
                if isinstance(func,dict):
                    try:
                        pdtemp.data[tag]['x']=func['x'](pdtemp.data[tag]['x'])
                        pdtemp.data[tag]['y']=func['y'](pdtemp.data[tag]['y'])
                    except KeyError,error:
                        print error
                        print """func should be dictionary of ('x'/'y',
                            x_func/y_func) pairs."""
                else:
                    raise ValueError("""func should be a dictionary
                        by using ('x'/'y', x_func/y_func) pairs.""")
            else:
                raise ValueError("Unknown axis value")

        return pdtemp


    def _data_complete_check_by_tag(self,tag):
        """
        Check if all the 6 base keys with value not as None has the same
            length and return key list with valid and invalid(None) value.
        """
        tagdic=Dic_Extract_By_Subkeylist(self.data[tag],Pdata._data_base_keylist)
        #get valid key and value list
        valid_key_list=[]
        valid_value_list=[]
        for key,value in tagdic.items():
            if value==None:
                pass
            else:
                valid_key_list.append(key)
                valid_value_list.append(value)
        #check if valid_key_list have unequal length of value
        for key1 in valid_key_list:
            for key2 in valid_key_list:
                if key1!=key2:
                    len1=len(tagdic[key1])
                    len2=len(tagdic[key2])
                    if len1!=len2:
                        raise ValueError("length of '{0}' data is {1} but lenghth of '{2}' data is {3} in tag '{4}'".format(key1,len1,key2,len2,tag))
        invalid_key_list=StringListAnotB(tagdic.keys(),valid_key_list)
        return valid_key_list,invalid_key_list

    def _data_complete_check_all(self):
        for tag in self._taglist:
            _,_=self._data_complete_check_by_tag(tag)

    def _get_err_by_tag(self,tag):
        """
        Get xerr,yerr for tag; if both errl and errh is not None, the nx2
            array returned; if errl!=None & errh==None,nx1 array returned;
            otherwise None returned.

        Notes:
        ------
        We assume errl=None while errh!=None will not occur.
        """
        tag_data=self.data[tag]
        valid_key_list,invalid_key_list=self._data_complete_check_by_tag(tag)
        #here we assume when equal low and high end of error is used, it's
        #been assigned to x/yerrl with x/yerrh as None. if you want to set 
        #error bar for only side, then set the other one explicitly to 0.
        if 'xerrl' in invalid_key_list and 'xerrh' in invalid_key_list:
            xerr=None
        elif 'xerrl' in invalid_key_list and 'xerrh' not in invalid_key_list:
            raise ValueError("strange that 'xerrl' is None but 'xerrh' not for tag '{0}', set one side as zero if you want single-side errorbar".format(tag))
        elif 'xerrl' not in invalid_key_list and 'xerrh' in invalid_key_list:
            xerr=tag_data['xerrl']
        else:
            xerr=np.ma.vstack((tag_data['xerrl'],tag_data['xerrh']))

        if 'yerrl' in invalid_key_list and 'yerrh' in invalid_key_list:
            yerr=None
        elif 'yerrl' in invalid_key_list and 'yerrh' not in invalid_key_list:
            raise ValueError("strange that 'yerrl' is None but 'yerrh' not for tag '{0}', set one side as zero if you want single-side errorbar".format(tag))
        elif 'yerrl' not in invalid_key_list and 'yerrh' in invalid_key_list:
            yerr=tag_data['yerrl']
        else:
            yerr=np.ma.vstack((tag_data['yerrl'],tag_data['yerrh']))

        return xerr,yerr

    @staticmethod
    def _expand_by_keyword(taglist,list_tagkw_tagkwvalue_tuple):
        """
        Purpose: convert a list of [(tag_kw1,v1),(tag_kw2,v2)] tuples to a
            list like [(tag1_kw1:v1),(tag2_kw1,v1),(tag3_kw2,v2),
            (tag4_kw2,v2),...], where tag1_kw1,tag2_kw1 are tags containg
            keyword tag_kw1, tag3_kw2,tag4_kw2 are tags containg keyword
            tag_kw2,..
        list_tag_tagvalue_tuple is like [('dry','r'),('wet','b')] where
            'dry' and 'wet' are keywords used for classifying tags.

        Examples:
        ---------
        >>> tags = ['wet1', 'wet2', 'wet3', 'dry1', 'dry3', 'dry2']
        >>> Pdata._expand_by_keyword(tags,[('dry','r'),('wet','b')])
        [('dry1', 'r'), ('dry3', 'r'), ('dry2', 'r'), ('wet1', 'b'),
         ('wet2', 'b'), ('wet3', 'b')]
        """
        full_tagkw_list=[]
        for tagkw,tagkwvalue in list_tagkw_tagkwvalue_tuple:
            for tag in FilterStringList(tagkw,taglist):
                full_tagkw_list.append((tag,tagkwvalue))
        return full_tagkw_list


    @staticmethod
    def _expand_tag_value_to_dic(taglist,tag_attr_value,tagkw=False):
        '''
        Check notes for Pdata.add_attr_by_tag for more details.
        '''
        #This if/else build the final dic to be used.
        if not isinstance(tag_attr_value,dict):
            if isinstance(tag_attr_value,list):
                #tag_attr_value is a list of (tag,attr_value) tuples.
                if isinstance(tag_attr_value[0],tuple):
                    #tag is keyword
                    if tagkw==True:
                        expand_tag_attr_value_list = \
                            Pdata._expand_by_keyword(taglist,
                                                     tag_attr_value)
                        final_dic=dict(expand_tag_attr_value_list)
                    else:
                        final_dic=dict(tag_attr_value)
                #a list of values
                else:
                    if len(tag_attr_value)!=len(taglist):
                        raise ValueError('''taglist has len '{0}' but input
                            list len is {1}'''.format(len(taglist),
                            len(tag_attr_value)))
                    else:
                        final_dic=dict(zip(taglist,tag_attr_value))
            #assume a single value (number or string)
            else:
                final_dic=dict(zip(taglist, len(taglist)*
                               [tag_attr_value]))
        #tag_attr_value is a dict
        else:
            final_dic=tag_attr_value
        return final_dic

    def add_attr_by_tag(self,tagkw=False,**nested_attr_tag_value_dic):
        """
        Add extra base attribute or ploting attribute by using key/keyvalue
            pairs, keys can principally be any keys of self._taglist; but this
            method is suggested to be used to add attr_name in
            _extra_base_keylist or _plot_attr_keylist_all

        Notes:
        ------
        Note the keyvalue is very flexible:
        1. In case of a single value, it will be broadcast to all tags for the
            attr concerned.
        2. In case of a list with lenght equal to number of tags, it will be
            add to tags by sequence of Pdata._taglist
        3. In case of a dictionary of tag/attr_value pairs, add attr_value
            accordingly to the tag corresponded.
        4. In case of a list of (tag,value) tuples:
            4.1 if tagkw==True (tagkw is only for this purpose):
                treat the tag in tag/attr_value as tag_keyword to set the same
                attr_value for all tags that contains this tag_keyword.
            4.2 if tagkw==False:
                treat the tag in tag/attr_value as a full tag and will not do
                the keyword search, it will change the tuple directly to a
                dictionary and apply the dictionary in seting tag/attr_value
                directly.

        Available keys are:
        **extra base attribute:
            bleftshift --> bar plot left shift, when bars are very close to
                each other with identical x values, some could be shift
                leftward to seperate them.(default 0)
            bwidth --> barplot width,(default: 0.5)
            bbottom --> barplot bottom (default:0)
            label
        **scatter plot:
             ssize
             scolor
             smarker
             scmap
             snorm
             svmin
             svmax
             + all other kwargs in axes.scatter()
         **errorbar plot:
             'efmt':None
             'ecolor':'b'
             'elinewidth':None
             'capsize':3
             'barsabove':False
             'lolims':False
             'uplims':False
             'xlolims':False
             'xuplims':False
             + all other kwargs in axes.errorbar()
        **plot plot:
        **bar:
            barcolor
        """
        for attr_name,tag_attr_value in nested_attr_tag_value_dic.items():
            final_dic = Pdata._expand_tag_value_to_dic(self._taglist,
                                                     tag_attr_value, tagkw)
            #apply value to tag
            for tag,attr_value in final_dic.items():
                if tag not in self.data:
                    raise ValueError("tag '{0}' dose not exist".format(tag))
                else:
                    self.data[tag][attr_name]=attr_value

    def set_default_plot_attr(self):
        pass

    def _fill_errNone_with_Nan(self,tag):
        tag_data=self.data[tag].copy()
        data_len=len(tag_data['x'])
        if tag_data['xerrl']==None and tag_data['xerrh']==None:
            tag_data['xerrl']=tag_data['xerrh']=np.repeat(np.nan,data_len)
        elif tag_data['xerrl']!=None and tag_data['xerrh']==None:
            tag_data['xerrh']=tag_data['xerrl']
        elif tag_data['xerrl']==None and tag_data['xerrh']!=None:
            raise ValueError('''err low end is None but high end not for xerr
                in tag '{0}' '''.format(tag))
        else:
            pass
        if tag_data['yerrl']==None and tag_data['yerrh']==None:
            tag_data['yerrl']=tag_data['yerrh']=np.repeat(np.nan,data_len)
        elif tag_data['yerrl']!=None and tag_data['yerrh']==None:
            tag_data['yerrh']=tag_data['yerrl']
        elif tag_data['yerrl']==None and tag_data['yerrh']!=None:
            raise ValueError('''err low end is None but high end not for yerr
                in tag '{0}' '''.format(tag))
        else:
            pass
        return tag_data

    @staticmethod
    def Hstack_Dic_By_Key(dict_list,key_list):
        """
        Horizontal stack dic values from dict list for keys present in
            key_list. Return a dict.
        Dictionaries in dict_list must have exactly the same keys with
            value as np.ndarray type
        """
        outdic={}
        for key in key_list:
            value=[]
            for subdict in dict_list:
                value.append(subdict[key])
            outdic[key]=np.ma.hstack(tuple(value))
        return outdic

    def Vstack_By_Tag(self,tagseq=None,axis='y'):
        """
        Vertial stack the data into numpy array for the tags in tagseq, the
            1st tag' data is on bottom output array.

        Notes:
        ------
        1. implicite case to call this function is tags have sharex data.
        """
        data_list = [self.data[tag][axis] for tag in tagseq]
        #we need to reverse the data_list as
        #we want the first come-in to stay on the bottom
        return np.ma.vstack(data_list[::-1])

    def pool_data_by_tag(self,tagkw=False,**group_dic):
        """
        Pool the data together by specifying new_tag=[old_tag_list] pairs;
            when tagkw==True, it's allowed to use new_tag=old_tag_keyword to
            specify what old tags will be pooled together which will
            include old_tag_keyword. Other attributes outside the
            _data_base_keylist will be copied from the first old_tag of
            the old_tag_list.
        Example:
            pool_data_by_tag(alldry=['1stdry','2nddry','3rddry'],
                allwet=['wet1','wet2'])
            pool_data_by_tag(tagkw=True,alldry='dry')
        """
        if tagkw==True:
            group_dic_final={}
            tags=self._taglist
            for newtag,newtag_tagkw in group_dic.items():
                group_dic_final[newtag]=FilterStringList(newtag_tagkw,tags)
        else:
            group_dic_final=group_dic

        pdata=Pdata()
        for newtag in group_dic_final.keys():
            old_tag_list=group_dic_final[newtag]
            old_entry_list=[self._fill_errNone_with_Nan(tag) for
                                tag in old_tag_list]
            new_entry=Pdata.Hstack_Dic_By_Key(old_entry_list,
                                              Pdata._data_base_keylist)
            #for _extra_base_keylist attributes, their default value in
            #_new_entry will be supplied to the new tag (pooled) data.
            #all other ploting attributes in _plot_attr_dic.values(),
            #they're lost.
            first_old_tag=old_tag_list[0]
            new_entry.update(Dic_Remove_By_Subkeylist(self.data[first_old_tag],
                                                      Pdata._data_base_keylist))
            pdata.add_entry_by_dic(**{newtag:new_entry})
        return pdata

    def regroup_data_by_tag(self,taglist):
        """
        Subset data by "taglist" and return as a new Pdata instance.
        With all features for tag reserved.
        """
        targetdic=Dic_Extract_By_Subkeylist(self.data,taglist)
        pdata=Pdata(targetdic)
        return pdata

    def regroup_data_by_tag_keyword(self,tag_keyword):
        tags=self._taglist
        #return self.regroup_data_by_tag(FilterStringList(tag_keyword,tags))

    def leftshift(self,shift=0,taglist='all'):
        taglist=self._set_default_tag(taglist)
        for tag in taglist:
            self.data[tag]['x']=self.data[tag]['x']-shift

    #a wrapper of scatter plot function
    def _gscatter(self,axes,x,y,**kwargs):
        attr_dic={             \
                 'ssize':20,   \
                 'scolor':'k', \
                 'smarker':'o',\
                 'scmap':None, \
                 'snorm':None, \
                 'svmin':None, \
                 'svmax':None, \
                 }
        attr_dic.update(Dic_Extract_By_Subkeylist(kwargs,Pdata._scatter_attr_base_keylist))
        remain_kwargs=Dic_Remove_By_Subkeylist(kwargs,Pdata._plot_attr_keylist_all)
        print 'attr_dic kwargs passed to axes.scatter()',attr_dic
        print 'scatter kwargs passed to axes.scatter()',remain_kwargs
        return axes.scatter(x, y,s=attr_dic['ssize'], c=attr_dic['scolor'], marker=attr_dic['smarker'], cmap=attr_dic['scmap'], \
                   norm=attr_dic['snorm'], vmin=attr_dic['svmin'], vmax=attr_dic['svmax'], alpha=None, linewidths=None, faceted=True, verts=None, **remain_kwargs)

    def scatter(self,axes=None,erase=True,**kwargs):
        """
        Make a scatter plot.
        Arugments:
            erase : set erase as True is a scatter plot has already been made for this pd object and you want to make a new scatter plot and "erase" the existing one
        """
        axes=_replace_none_axes(axes)
        if erase==True:
            if hasattr(self,'Scatter_PathC') and self.Scatter_PathC != {}:
                self.remove_scatter_by_tag()
        self._data_complete_check_all()
        if Is_Nested_Dic(kwargs):
            raise ValueError('tag specific plot attributes should be set by using add_attr_by_tag')
        else:
            self.Scatter_PathC={}  #Scatter_PathC ---> Scatter PathCollection
            for tag,tag_data in self.data.items():
                tag_kwargs=kwargs.copy()
                tag_plot_attr_dic=Dic_Remove_By_Subkeylist(tag_data,Pdata._all_nonplot_keylist)
                tag_kwargs.update(tag_plot_attr_dic)
                print 'tag & tag_kwargs',tag,tag_kwargs
                self.Scatter_PathC[tag]=self._gscatter(axes,tag_data['x'],tag_data['y'],**tag_kwargs)
        return self.Scatter_PathC

    def scattertag(self,axes,taglist='all',**kwargs):
        taglist = self._set_default_tag(taglist)
        if hasattr(self,'Scatter_PathC') and self.Scatter_PathC!={}:
            self.remove_scatter_by_tag(taglist)
        self._data_complete_check_all()
        if Is_Nested_Dic(kwargs):
            raise ValueError('tag specific plot attributes should be set by using add_attr_by_tag')
        else:
            self.Scatter_PathC={}  #Scatter_PathC ---> Scatter PathCollection
            for tag,tag_data in Dic_Extract_By_Subkeylist(self.data,taglist).items():
                tag_kwargs=kwargs.copy()
                tag_plot_attr_dic=Dic_Remove_By_Subkeylist(tag_data,Pdata._all_nonplot_keylist)
                tag_kwargs.update(tag_plot_attr_dic)
                print 'tag & tag_kwargs',tag,tag_kwargs
                self.Scatter_PathC[tag]=self._gscatter(axes,tag_data['x'],tag_data['y'],label=tag_data['label'],**tag_kwargs)
        return self.Scatter_PathC

    def plot(self,axes=None,**kwargs):
        axes=_replace_none_axes(axes)
        if hasattr(self,'Line2D') and self.Line2D!={}:
            self.remove_line_by_tag()
        self.Line2D={}
        for tag,tag_data in self.data.items():
            line2D=axes.plot(tag_data['x'],tag_data['y'],label=tag_data['label'],**kwargs)
            self.Line2D[tag]=line2D[0]
        return self.Line2D

    def plot_stackline(self,axes=None,tagseq=None,colors=None,bottom_fill=True,legend=True,legdic={},fillkw={}):
        """
        Make stacked line plot for sharex pdata.

        Parameters:
        -----------
        tagseq: the tag list for which stacked line plot will be made. Notice the sequece for tagseq is from bottom to the top.
        colors: color list, the length should be equal to the number of filled area. In case of bottom_fill == True, len(colors) should be equal to len(tagseq),
            otherwise should be equal to len(tagseq)-1.
        bottom_fill: set True if the area between xaxis and the bottom line (the first tag) is to be filled.
        legdic: the legend is set by using g.ProxyLegend(), legdic could be passed as kwargs for the function of g.ProxyLegend.create_legend
        kwargs: used for plt.fill_between functions.

        Returns:
        --------
        pdata,proleg
        pdata: the new pdata that stores the cumsum 'y' data.
        proleg: the g.ProxyLegend object that's used to creat legend for the stacked plots.

        Notes:
        ------
        1. the **kwargs is also passed directly as kwargs for the function mat.patches.Rectangle and this may lead to conflicts in some special cases.

        See also:
        ---------
        g.ProxyLegend, mat.patches.Rectangle
        """
        #build the new pdata that stores the cumsum data
        tagseq = _replace_none_by_given(tagseq,self._taglist)
        arr = self.Vstack_By_Tag(tagseq)
        arr_cumsum = arr[::-1].cumsum(axis=0)[::-1]  #we want the cumsum with direction from bottom-->top.
        pdnew = self.copy()
        for index,tag in enumerate(tagseq):
            pdnew.data[tag]['y'] = arr_cumsum[::-1][index]

        axes=_replace_none_axes(axes)
        #treat case of bottom_fill
        if bottom_fill == True:
            colors = _replace_none_colorlist(colors,num=len(tagseq))
            leg_tagseq = tagseq[:]
        else:
            colors = _replace_none_colorlist(colors,num=len(tagseq)-1)
            leg_tagseq = tagseq[1:]

        proleg = g.ProxyLegend()
        if bottom_fill == True:
            bottom_tag = tagseq[0]
            bottom_data = pdnew.data[bottom_tag]
            axes.fill_between(bottom_data['x'],np.zeros(len(bottom_data['x'])),bottom_data['y'],color=colors[0],**fillkw)
            proleg.add_rec_by_tag(bottom_tag,color=colors[0],**fillkw)
            colors_remain = colors[1:]
        else:
            colors_remain = colors

        for index in range(len(tagseq)-1):
            below_tag = tagseq[index]
            focus_tag = tagseq[index+1]
            below_data = pdnew.data[below_tag]
            focus_data = pdnew.data[focus_tag]
            axes.fill_between(focus_data['x'],below_data['y'],focus_data['y'],color=colors_remain[index],**fillkw)
            proleg.add_rec_by_tag(focus_tag,color=colors_remain[index],**fillkw)
            if index == len(tagseq)-2:
                axes.set_xlim(focus_data['x'][0],focus_data['x'][-1])

        if legend:
            proleg.create_legend(axes,tagseq=leg_tagseq[::-1],**legdic)
        else:
            pass
        self.stackline_proleg = proleg
        self.stackline_actual_pd = pdnew
        self.axes = axes

    def _gerrorbar(self,axes,x,y,yerr=None,xerr=None,**kwargs):
        attr_dic={             \
                 'efmt':None,\
                 'ecolor':'r',   \
                 'elinewidth':None, \
                 'capsize':3, \
                 'barsabove':False, \
                 'lolims':False,\
                 'uplims':False,\
                 'xlolims':False,\
                 'xuplims':False \
                 }
        attr_dic.update(Dic_Extract_By_Subkeylist(kwargs,Pdata._error_attr_base_keylist))
        remain_kwargs=Dic_Remove_By_Subkeylist(kwargs,Pdata._plot_attr_keylist_all)
        print 'attr_dic in passed to errorbar()',attr_dic
        print 'kwargs passed to errorbar()',remain_kwargs
        return axes.errorbar(x, y, yerr=yerr, xerr=xerr, fmt=attr_dic['efmt'], ecolor=attr_dic['ecolor'], elinewidth=attr_dic['elinewidth'], \
                             capsize=attr_dic['capsize'], barsabove=attr_dic['barsabove'], lolims=attr_dic['lolims'], uplims=attr_dic['uplims'], \
                             xlolims=attr_dic['xlolims'], xuplims=attr_dic['xuplims'], **remain_kwargs)

    def _set_ecolor_by_scatter(self):
        taglist=[]
        colorlist=[]
        for tag in self.data.keys():
            if tag in self.Scatter_PathC:
                taglist.append(tag)
                try:
                    colorlist.append(self.Scatter_PathC[tag].get_facecolor()[0])
                except IndexError:
                    print "the scatter facecolor for tag '{0}' is none, scatter edgecolor is used as errorbar color.".format(tag)
                    colorlist.append(self.Scatter_PathC[tag].get_edgecolor()[0])
        self.add_attr_by_tag(ecolor=dict(zip(taglist,colorlist)))
        return dict(zip(taglist,colorlist))

    def _set_ecolor_by_bar(self):
        taglist=[]
        colorlist=[]
        for tag in self.data.keys():
            if tag in self.Bar_Container:
                taglist.append(tag)
                try:
                    colorlist.append(self.Bar_Container[tag][0].get_facecolor())
                except IndexError:
                    print "the bar facecolor for tag '{0}' is none, bar edgecolor is used as errorbar color.".format(tag)
                    colorlist.append(self.Bar_Container[tag][0].get_edgecolor())
        self.add_attr_by_tag(ecolor=dict(zip(taglist,colorlist)))

    @staticmethod
    def _get_attr_value_from_handle(handle,attr_name):
        """
        Return an attribute from the specified handle.
        Example:
            _get_attr_from_handle(matplotlib.collections.PathCollection,'_facecolors') is the same as matplotlib.collections.PathCollection.get_facecolor()
        """
        return handle.__getattribute__(attr_name)

    @staticmethod
    def _get_single_attr_by_check_similar_attrname_and_condition(handle, func_attr_condition_tuple_list):
        """
        Purpose: Return a single attribute by checking sequentially the attribute name list and if the condition has been made.
        FirstCase: We want to make the program as intelligent as possible, if we want to set the errorbar color according to the facecolor or edgecolor of the scatters.
            Sometimes the facecolor of a scatter could be None, then we want to use the edgecolor as the errorbar color. In this case, we need to set some priority: if 
            the facecolor is not None, then the facecolor will be returned; otherwise the edgecolor will be returned.
        Arguments:
            1.  attr_condition_tuple_list is a list of 3-length tuple (func, attr_name, discard_condition_value). For the tuple, the first element of is a function, the second element is the attribute name, and and the third element is used to check if the attribute value should be discard. Sometimes the attributed value could be checked directly as equal to equal to the given (third) value or not (in this case the first element will be give as None), sometimes it's not possible to check this directly, so we need to use a function (first tuple element) to check if the condtion to discard the attribute value is met. If the condition is met (either through a function) or by equal relation operation directly, the
                attribute value will be discard.
            2   If all given conditions for all the attribute names are met, then an error will be raised.
        Example:
            _get_single_attr_by_check_similar_attrname_and_condition(matplotlib.collections.PathCollection, [(__builtin__.len,'_facecolors',0),(__builtin__.len,'_edgecolors',0)]). This means if both facecolor and edgecolor is np.array([]), then an error will be raised. But if facecolor is not np.array([]), then the color of the facecolor will be used (This valid attribute value that stays earlier has higher priority in determining the final return value).
        Doctest:
            >>> fig, ax = g.Create_1Axes()
            >>> a = np.arange(10)
            >>> ax.scatter(a,a,facecolor=np.array([ 0.,  0.,  1.,  1.]),edgecolor=np.array([ 1.,  0.,  0.,  1.]))
            >>> _get_single_attr_by_check_similar_attrname_and_condition(col,[(__builtin__.len,'_facecolors',0),(__builtin__.len,'_edgecolors',0)])[0] == np.array([ 0.,  0.,  1.,  1.])
            >>> _get_single_attr_by_check_similar_attrname_and_condition(col,[(__builtin__.len,'_edgecolors',0),(__builtin__.len,'_facecolors',0)]) == np.array([ 1.,  0.,  0.,  1.])
        """
        error_state = True
        for (func, attr_name, discard_condition_value) in func_attr_condition_tuple_list:
            valid_attribute_value = False
            attribute_value = Pdata._get_attr_value_from_handle(handle,attr_name)
            #attribute value and discard_condition_value could be compared directly
            if func == None:
                if attribute_value == discard_condition_value:
                    pass
                else:
                    valid_attribute_value = True
            #we needs a function to check the attribute value is valid
            else:
                if func(attribute_value) == discard_condition_value:
                    pass
                else:
                    valid_attribute_value = True

            #we return the valid attribute value
            if valid_attribute_value:
                error_state = False
                return attribute_value
        if error_state:
            raise ValueError("not valid attribute value found for the list of attributes: {0}".format(zip(*attr_condition_tuple_list)[0]))

    def errorbar(self,axes,ef='scatter',**kwargs):
        """
        ef is short for 'ecolor follow', ef could be 'scatter','bar', if None, ecolor in self.data[tag]['ecolor'] is used, or it not present in self.data,
        _error_attr_default value will be used.
        """
        self._data_complete_check_all()
        if ef=='scatter':
            self._set_ecolor_by_scatter()
        elif ef == 'bar':
            self._set_ecolor_by_bar()
        elif ef=='none':
            pass
        else:
            raise ValueError("incorrect ef value '{0}'".format(ef))
        if Is_Nested_Dic(kwargs):
            raise ValueError('tag specific plot attributes should be set by using add_attr_by_tag')
        else:
            self.Errorbar_Lines={}
            for tag,tag_data in self.data.items():
                tag_xerr,tag_yerr=self._get_err_by_tag(tag)
                if (tag_xerr==None and tag_yerr==None):
                    continue
                else:
                    print tag_xerr,tag_yerr
                    tag_kwargs=kwargs.copy()
                    tag_plot_attr_dic=Dic_Remove_By_Subkeylist(tag_data,Pdata._all_nonplot_keylist)
                    tag_kwargs.update(tag_plot_attr_dic)  #here, tag_plot_attr_dic contains keys for all ploting types; tag_kwargs get contaminated.
                    self.Errorbar_Lines[tag]=self._gerrorbar(axes,tag_data['x'],tag_data['y'],xerr=tag_xerr,yerr=tag_yerr,**tag_kwargs)
        return self.Errorbar_Lines

    def _gbar(self,axes,left,height,width=0.5,bottom=0,label=None,**kwargs):
        attr_dic=Pdata._bar_attr_default.copy()
        attr_dic.update(Dic_Extract_By_Subkeylist(kwargs,Pdata._bar_attr_base_keylist))
        kwargs=Dic_Remove_By_Subkeylist(kwargs,Pdata._plot_attr_keylist_all)
        print '_gbar',kwargs
        print '_bar',attr_dic
        return axes.bar(left, height, width=width, bottom=bottom,color=attr_dic['bcolor'], label=label, **kwargs)

    def bar(self,axes=None,**kwargs):
        axes=_replace_none_axes(axes)
        self._data_complete_check_all()
        if Is_Nested_Dic(kwargs):
            raise ValueError('tag specific plot attributes should be set by using add_attr_by_tag')
        else:
            self.Bar_Container={}  #Bar_Container --> <Container object of 10 artists> with matplotlib.patches.Rectangle as list members.
            for tag,tag_data in self.data.items():
                tag_kwargs=kwargs.copy()
                tag_plot_attr_dic=Dic_Remove_By_Subkeylist(tag_data,Pdata._all_nonplot_keylist)
                tag_kwargs.update(tag_plot_attr_dic)
                print 'tag & tag_kwargs',tag,tag_kwargs
                self.Bar_Container[tag]=self._gbar(axes,tag_data['x']-tag_data['bleftshift']-tag_data['bwidth']*0.5,tag_data['y'],width=tag_data['bwidth'],\
                                                   bottom=tag_data['bbottom'],\
                                                   label=tag_data['label'],**tag_kwargs)
        return self.Bar_Container


    def _get_plot_attr_value_from_data(self,*attr_list):
        """
        Get plot attribute value from data. attr_list must be subset of Pdata._plot_attr_keylist_all
        """
        attr_dic={}
        for attr_name in attr_list:
            attr_dic[attr_name]=dict()
            for tag,tag_data in self.data.items():
                try:
                    attr_dic[attr_name][tag] = tag_data[attr_name]
                except KeyError:
                    attr_dic[attr_name][tag] = Pdata._plot_attr_default[attr_name]
        return attr_dic

    def _get_scatter_marker_size_color(self):
        self.outdic={}
        self.outdic['color']=_get_attr_value_from_objins_dic(self.Scatter_PathC,'_facecolors').values()[0]
        self.outdic['size']=self._get_plot_attr_value_from_data('ssize').values()[0]
        self.outdic['marker']=self._get_plot_attr_value_from_data('smarker').values()[0]
        return self.outdic

    def _get_Scatter_PathC_marker_size_color(self):
        sc=self._get_scatter_marker_size_color()
        outdic={}
        for tag,ScatterP in self.Scatter_PathC.items():
            outdic[(tag,ScatterP)]=dict()
            for attr_name,attr_dic in sc.items():
                outdic[(tag,ScatterP)][attr_name]=attr_dic[tag]
        return outdic


    def _get_handler_label_by_artistdic(self,artist_dic,taglab=False,tag_seq=None):
        """
        Retrieve handler and label list from self.ScatterP/Line2D/BarContainer... for legend seting.
        When tag_seq is a list of tags, will return handler/label list as sorted by tags in tag_seq, note using tag_seq will restrict the scope of tags included in 
        legend.
        """
        handler_list=[]
        label_list=[]
        if tag_seq==None:
            tag_seq=artist_dic.keys()
        for tag in tag_seq:
            artist=artist_dic[tag]
            handler_list.append(artist)
            if taglab==True:
                label_list.append(tag)
            else:
                label_list.append(self.data[tag]['label'])
        return handler_list,label_list

    def _get_axes_by_artistdic(self,artist_dic):
        """
        Get which axes the artisti_dic is in
        """
        artist=artist_dic.values()[0]
        try:
            return artist.get_axes()
        except AttributeError:
            return artist[0].get_axes()

    def set_legend_scatter(self,axes=None,taglab=False,tag_seq=None,**kwargs):
        tag_seq = _replace_none_by_given(tag_seq,self._taglist)
        handler_list,label_list=self._get_handler_label_by_artistdic(self.Scatter_PathC,taglab=taglab,tag_seq=tag_seq)
        if axes==None:
            axes=self._get_axes_by_artistdic(self.Scatter_PathC)
        self.LegendScatter=axes.legend(handler_list,label_list,**kwargs)
        return self.LegendScatter

    def set_legend_line(self,axes=None,taglab=False,tag_seq=None,**kwargs):
        handler_list,label_list=self._get_handler_label_by_artistdic(self.Line2D,taglab=taglab,tag_seq=tag_seq)
        if axes==None:
            axes=self._get_axes_by_artistdic(self.Line2D)
        self.LegendLine=axes.legend(handler_list,label_list,**kwargs)
        return self.LegendLine

    def set_legend_bar(self,axes=None,taglab=False,tag_seq=None,**kwargs):
        handler_list,label_list=self._get_handler_label_by_artistdic(self.Bar_Container,taglab=taglab,tag_seq=tag_seq)
        if axes==None:
            axes=self._get_axes_by_artistdic(self.Bar_Container)
        self.LegendBar=axes.legend(handler_list,label_list,**kwargs)

    def set_data_void(self):
        for tag, tag_value in self.data.items():
            tag_value['x'] = []
            tag_value['y'] = []


    def set_legend_all(self,axes=None,taglab=False,tag_seq=None,**kwargs):
        handler_list=[]
        label_list=[]
        attr_list=['Scatter_PathC','Line2D','Bar_Container']
        artdic_list=[self.__getattribute__(attr) for attr in attr_list if hasattr(self,attr)]
        for artdic in artdic_list:
            try:
                sub_hand_list,sub_lab_list=self._get_handler_label_by_artistdic(artdic,taglab=taglab,tag_seq=tag_seq)
                handler_list.append(sub_hand_list)
                label_list.append(sub_lab_list)
            except AttributeError:
                pass
        if axes==None:
            axes=handler_list[0][0].get_axes()
        self.LegendAll=axes.legend(iteflat(handler_list),iteflat(label_list),**kwargs)

    def set_legend_select(self,axes=None,taglab=False,tag_seq=None,**kwargs):
        pass
    def _remove_artist_by_tag(self,artist_dic,taglist='all'):
        #this method from http://stackoverflow.com/questions/4981815/how-to-remove-lines-in-a-matplotlib-plot
        if taglist=='all':
            taglist=artist_dic.keys()
        for tag in taglist:
            if tag in artist_dic:
                art=artist_dic[tag]
                wl=weakref.ref(art)
                art.remove()
                del art
                del artist_dic[tag]

    def remove_line_by_tag(self,taglist='all'):
        self._remove_artist_by_tag(self.Line2D,taglist=taglist)

    def remove_scatter_by_tag(self,taglist='all'):
        self._remove_artist_by_tag(self.Scatter_PathC,taglist=taglist)

    def setp_by_tag(self,artist_dic,**nested_attr_tag_value_dic):
        """
        Set artist property by attr_name=single value, or attr_name={tag:attr_value} pairs. In case of single value provided, it's broadcast to apply on all artists.
        In case of a {tag:attr_value} dictionary, attr_value applied according to its tat key.  Note the {tag:attr_value}.keys() must be exactly the same as 
        in artist_dic.keys(). 
        """
        for attr_name,tag_attr_value_dic in nested_attr_tag_value_dic.items():
            if not isinstance(tag_attr_value_dic,dict):
                for tag,art in artist_dic.items():
                    plt.setp(art,**{attr_name:tag_attr_value_dic})
            else:
                for tag,attr_value in tag_attr_value_dic.items():
                    plt.setp(artist_dic[tag],**{attr_name:attr_value})

    def creat_list_of_axes_by_tagseq(self,force_axs=None,tagseq=None,ncols=1,sharex=True,**subplot_kwargs):
        """
        """
        tag_list=_replace_none_by_given(tagseq,self._taglist)
        num=len(tag_list)
        axs = _build_list_of_axes_by_num(num,force_axs=force_axs,ncols=ncols,sharex=sharex,**subplot_kwargs)
        return dict(zip(tag_list,axs))

    def plot_split_axes(self,force_axs=None,force_taglist=None,force_axdic=None,ncols=1,sharex=True,tagpos='ul',unit=None,xlim=None,plotkw={},**subplot_kwargs):
        if force_axdic == None:
            axdic = self.creat_list_of_axes_by_tagseq(force_axs=force_axs,tagseq=force_taglist,ncols=ncols,sharex=sharex,**subplot_kwargs)
        else:
            axdic = force_axdic

        for tag,axt in axdic.items():
            pd_temp=self.regroup_data_by_tag([tag])
            pd_temp.plot(axt,**plotkw)
        _treat_axes_dict(axdic,tagpos=tagpos,unit=unit,xlim=xlim)
        self.axes = axdic

    def imshow(self,force_axs=None,tagseq=None):
        pass

    def plot_split_axes_byX(self,force_axs=None,sharex=True,tagpos='ul',unit=False,ylab=False,force_taglist=None,ncols=1,plotkw={},**fig_kw):
        """
        line plot with each tag in a separate axes, data with different tags have the same xdata.
        Now only unit is used, ylab is not used.
        Arguments:
            force_taglist: could be used to plot for only selected tags, or to force the sequence of tags.
        """
        if force_taglist == None:
            tag_list=self._taglist
        else:
            tag_list=force_taglist
        num=len(tag_list)
        axdic={}
        if force_axs==None:
            if ncols == 1:
                fig,axs=plt.subplots(nrows=num, ncols=1, sharex=sharex, sharey=False, squeeze=True, subplot_kw=None, **fig_kw)
            else:
                nrows=num/ncols+1
                fig,axt=plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=False, squeeze=True, subplot_kw=None, **fig_kw)
                axs=axt.flatten()[0:num]
        else:
            if num<=len(force_axs):
                axs=force_axs[0:num]
        for tag,axt in zip(tag_list,axs):
            pd_temp=self.regroup_data_by_tag([tag])
            pd_temp.plot(axt,**plotkw)
            axdic[tag]=axt
            g.Set_AxText(axt,tag,tagpos)
            if unit==True:
                if isinstance(unit,str):
                    print "forced unit is used"
                    axt.set_ylabel(unit)
                else:
                    try:
                        axt.set_ylabel(pd_temp.data[tag]['unit'])
                    except AttributeError:
                        pass
        if force_axs==None:
            return fig,axdic
        else:
            return axdic

    def to_dic(self,sharex=False):
        outdic={}
        if sharex==True:
            for tag,tag_dic in self.data.items():
                outdic[tag]=tag_dic['y']
        else:
            raise ValueError("This can only be use for sharex=True, please specify this explicitly or check data length for each tag.")
        return outdic

    def to_DataFrame(self,tagindex=False,col_name=None,df_index=None):
        """
        Change the data to DataFrame

        Parameters:
        ----------
        tagindex : True if the tags are used as index of the new DataFrame, in this case the col_name must be provided. False if the tags are used as DataFrame column names.
        """
        outdic=self.to_dic(sharex=True)
        if tagindex == False:
            if df_index == None:
                index=self.data.values()[0]['x'] #as sharex is True, so we can pick x data for any tag and use it as index of the DataFrame.
            else:
                index=None
            df = pa.DataFrame(outdic,index=index) #as sharex is True, so we can pick x data for any tag and use it as index of the DataFrame.

        else:
            if col_name == None:
                raise ValueError("In case of using tag as DataFrame index, new columns names must be provided! col_name could not be None!")
            else:
                df=exio.DataFrame_from_Dic_Key_as_Index(outdic, columns=col_name)
        return df

    def set_new_tags(self, old_new_tag_tuple_list):
        """
        Change the old tag to new tag according to old_new_tag_tuple_list

        Parameters:
        -----------
        old_new_tag_tuple_list: [(oldtag1,newtag1),(oldtag2,newtag2)]

        Notes:
        ------
        1. In-place operation.
        """
        for (oldtag,newtag) in old_new_tag_tuple_list:
            self.data[newtag] = self.data[oldtag]
            del self.data[oldtag]

    def to_pd(self,filename):
        """
        pickle the self.data to filename; the surfix '.pd' could be used as indication as pdata
        """
        pb.pfdump(self.data,filename)

def pmerge(*pdlist):
    pdnew=Pdata()
    for pd in pdlist:
        pdnew.data.update(pd.data)
    return pdnew

def open(pdfilename):
    """
    create a pdata from "*.pd" file
    """
    pdata=Pdata()
    data=pb.pfload(pdfilename)
    pdata.data=data
    return pdata

def dic_pdata_to_pd(pd_dic,filename):
    pd_data_dic={}
    for kw,pd in pd_dic.items():
        pd_data_dic[kw]=pd.data
    pb.pfdump(pd_data_dic,filename)

def dic_pdata_from_pd(pdfilename):
    pd_data_dic=pb.pfload(pdfilename)
    pd_dic={}
    for kw,pd_data in pd_data_dic.items():
        pd=Pdata()
        pd.data=pd_data
        pd_dic[kw]=pd
    return pd_dic


def plot_dic_as_pdata_sharex_noerror(indic,force_axs=None,sharex=True,tagpos='ul',unit=False,ylab=False,force_taglist=None,ncols=1,plotkw={},**fig_kw):
    pd = Pdata()
    pd.add_entry_sharex_noerror_by_dic(indic)
    pd.plot_split_axes_byX(force_axs=force_axs,sharex=sharex,tagpos=tagpos,unit=unit,ylab=ylab,force_taglist=force_taglist,ncols=ncols,plotkw={},**fig_kw)
    return pd

def read_csv(filepath_or_buffer,sep=',',header=0,index_col=None,index_func=None,names=None,date_parser=None,skiprows=None,na_values=None,df_func=None,
             force_sharex=None,**kwds):
    """
    A wrapper of pa.read_csv and Pdata.Pdata.add_entry_sharex_noerror_by_dic, to read csv directly to Pdata.Pdata

    Parameters:
    -----------
    force_sharex: use force_sharex to overwrite the default sharex of Pdata that are the index column of the csv/DataFrame object.

    Notes
    -----
    1. the index_col (index of the pandas DataFrame) will be treated as shared xaxis of the Pdata.Pdata object
    2. the (column) names of the pandas.DataFrame object will be used as tags in the Pdata.Pdata object.

    See also
    --------
    from_DataFrame
    """
    df = pa.read_csv(filepath_or_buffer, sep=sep, dialect=None, header=header, index_col=index_col, names=names, skiprows=skiprows, na_values=na_values, keep_default_na=True,                     thousands=None, comment=None, parse_dates=False, keep_date_col=False, dayfirst=False, date_parser=date_parser, nrows=None, iterator=False, 
                     chunksize=None, skip_footer=0, converters=None, verbose=False, delimiter=None, encoding=None, squeeze=False, **kwds)

    return from_DataFrame(df,df_func=df_func,index_func=index_func,force_sharex=force_sharex)

def from_DataFrame(df,df_func=None,index_func=None,force_sharex=None):
    """
    Create a sharex Pdata.Pdata object from pandas DataFrame

    Parameters:
    -----------
    df_func: function that applies on DataFrame before feeding data into Pdata.
    index_func: index function that will be applied before using the DataFrame index as shared xaxis of the Pdata object, this is useful as sometimes DataFrame index could
        be a bit strange and not readily compatible with matplotlib functions.
    force_sharex: In case index_func could not achieve the object to transform the index to desired sharex xaxis, force_sharex is used to force write the Pdata shared xaxis.
    """
    pd = Pdata()
    if df_func != None:
        df = df_func(df)
    if force_sharex == None:
        if index_func == None:
            pd.add_entry_sharex_noerror_by_dic(df,x=df.index)
        else:
            pd.add_entry_sharex_noerror_by_dic(df,x=index_func(df.index))
    else:
        pd.add_entry_sharex_noerror_by_dic(df,x=force_sharex)
    return pd

def _build_list_of_axes_by_num(num,force_axs=None,ncols=None,**kwargs):
    if force_axs==None:
        if ncols == 1:
            fig,axs=plt.subplots(nrows=num, ncols=1, **kwargs)
        else:
            if num%ncols == 0:
                nrows=num/ncols
            else:
                nrows=num/ncols+1
            fig,axt=plt.subplots(nrows=nrows, ncols=ncols, **kwargs)
            axs=axt.flatten()[0:num]
    else:
        if num<=len(force_axs):
            axs=force_axs[0:num]
        else:
            raise ValueError("given force_axs length is smaller than required.")
    return axs

def _treat_axes_dict(axdic,tagpos='ul',unit=None,xlim=None):
    """
    Label tag and unit for a dictionary of axes with tag as keys.
    """
    for tag,axt in axdic.items():
        g.Set_AxText(axt,tag,tagpos)
        if unit !=None :
            if isinstance(unit,str):
                print "forced unit is used"
                axt.set_ylabel(unit)
            else:
                raise ValueError("Strange unit")
        else:
            pass

        if xlim != None:
            axt.set_xlim(xlim)


class NestPdata(object):
    """
    NestPdata will receive a dictionary of Pdata object.
    """
    def __init__(self,dic_pdata):
        nestpdata_datadic = {}
        for parent_tag in dic_pdata.keys():
            pdtemp = dic_pdata[parent_tag]
            nestpdata_datadic[parent_tag] = pdtemp.data
        self.data = nestpdata_datadic
        self.parent_tags = self.data.keys()
        self.child_pdata = dic_pdata
        self.child_tags = self.child_pdata.values()[0].list_tags()

    def permuate_tag(self):
        self.data = pb.Dic_Nested_Permuate_Key(self.data)
        self.parent_tags = self.data.keys()
        self.child_pdata = {}
        for parent_tag in self.parent_tags:
            pdata_temp = Pdata()
            pdata_temp.add_entry_by_dic(**self.data[parent_tag])
            self.child_pdata[parent_tag] = pdata_temp
        self.child_tags = self.child_pdata.values()[0].list_tags()

    def creat_list_of_axes_by_tagseq(self,force_axs=None,tagseq=None,ncols=1,sharex=True,**subplot_kwargs):
        """
        """
        tag_list=_replace_none_by_given(tagseq,self.parent_tags)
        num=len(tag_list)
        axs = _build_list_of_axes_by_num(num,force_axs=force_axs,ncols=ncols,sharex=sharex,**subplot_kwargs)
        return dict(zip(tag_list,axs))

    def plot_split_parent_tag(self,force_axs=None,parent_tagseq=None,force_axdic=None,ncols=1,sharex=True,tagpos='ul',unit=None,xlim=None,plotkw={},**subplot_kwargs):
        """
        """
        if force_axdic == None:
            axdic = self.creat_list_of_axes_by_tagseq(force_axs=force_axs,tagseq=parent_tagseq,ncols=ncols,sharex=sharex,**subplot_kwargs)
        else:
            axdic = force_axdic

        for tag,axt in axdic.items():
            pd_temp = self.child_pdata[tag]
            pd_temp.plot(axt,**plotkw)

        _treat_axes_dict(axdic,tagpos=tagpos,unit=unit,xlim=xlim)
        self.axes = axdic


    def plot_split_parent_tag_stackline(self,force_axs=None,parent_tagseq=None,force_axdic=None,ncols=1,sharex=True,tagpos='ul',unit=None,xlim=None,fillkw={},plotkw={},**subplot_kwargs):
        """
        """
        if force_axdic == None:
            axdic = self.creat_list_of_axes_by_tagseq(force_axs=force_axs,tagseq=parent_tagseq,ncols=ncols,sharex=sharex,**subplot_kwargs)
        else:
            axdic = force_axdic

        for tag,axt in axdic.items():
            pd_temp = self.child_pdata[tag]
            pd_temp.plot_stackline(axes=axt,fillkw=fillkw,**plotkw)

        _treat_axes_dict(axdic,tagpos=tagpos,unit=unit,xlim=xlim)
        self.axes = axdic

    def copy(self):
        nestpd_dic = {}
        for tag,child_pd in self.child_pdata.items():
            nestpd_dic[tag] = child_pd.copy()
        nestpd = NestPdata(nestpd_dic)
        return nestpd

    @staticmethod
    def _regroup_parent_tag(nestpd,taglist):
        nestpd_dic = Dic_Extract_By_Subkeylist(nestpd.child_pdata,taglist)
        return NestPdata(nestpd_dic)

    def regroup_data_by_tag(self,taglist,keyword=True,mode='parent'):
        if mode == 'parent':
            return self._regroup_parent_tag(self,taglist)
        elif mode == 'child':
            nestpd_dic_pre = self.copy()
            nestpd_dic_pre.permuate_tag()  #change the child to parent
            nestpd = self._regroup_parent_tag(nestpd_dic_pre,taglist)
            nestpd.permuate_tag()  #change back
            return nestpd
        else:
            raise ValueError("unknown mode.")


