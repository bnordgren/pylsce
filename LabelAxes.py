#!/usr/bin/python

from collections import OrderedDict
import pandas as pa
import matplotlib as mat
import numpy as np
import g

def _replace_none_by_given(orinput,default):
    if orinput == None:
        return default
    else:
        return orinput

def _propagate(tags,arg,itearg=False):
    """
    Propagate the arg input to a (ordered)dict.

    when arg is a list with equal length of tags, the default tags
        sequence will be used.

    Parameters:
    -----------
    arg: the argument that's to be broadcasted.
    itearg: True if the arg itself is expected to be an iterable
    """
    tagnum = len(tags)
    if isinstance(arg,list):
        if itearg:
            if isinstance(arg[0],(tuple,list,np.ndarray)):
                if len(arg) != tagnum:
                    raise ValueError("""list length expected to be {0}"""
                                     .format(tagnum))
                else:
                    return OrderedDict(zip(tags,arg))
            else:
                return dict(zip(tags,[arg]*tagnum))
        else:
            if len(arg) != tagnum:
                raise ValueError("""list length expected to be {0}"""
                                 .format(tagnum))
            else:
                return OrderedDict(zip(tags,arg))
    elif isinstance(arg,dict):
        return arg
    else:
        return dict(zip(tags,[arg]*tagnum))


class LabelAxes(object):
    def __init__(self,tags=None,axl=None):
        if not isinstance(tags,(list,tuple,np.ndarray)):
            raise TypeError("tags could only be list or tuple numpy array")
        else:
            if not isinstance(axl[0],mat.axes.Axes):
                raise TypeError('value must be mat.axes.Axes type!')
            else:
                if len(tags) > len(axl):
                    raise ValueError("length of axl less than tags")
                else:
                    self.data = OrderedDict(zip(tags,axl))
                    self.tags = tags
                    self.tagnum = len(self.tags)
                    self.axl = list(axl)[0:self.tagnum]


    def _get_tags(self,key):
        #normal key:
        if isinstance(key,str):
            return [key]
        elif isinstance(key,int):
            return [self.tags[key]]
        else:
            if isinstance(key, slice):
                return self.tags[key]
            elif isinstance(key, list):
                if len(np.unique(map(type,key))) > 1:
                    raise TypeError("input list must be single type")
                else:
                    if isinstance(key[0],str):
                        return key
                    elif isinstance(key[0],int):
                        return [self.tags[index-1] for index in key]
                    else:
                        raise TypeError("slice not understood.")
            else:
                raise TypeError("slice not understood.")

    def __getitem__(self,key):
        subtags = self._get_tags(key)
        subvals = [self.data[stag] for stag in subtags]
        if len(subtags) == 1:
            return subvals[0]
        else:
            return LabelAxes(tags=subtags,axl=subvals)

    def __len__(self):
        return self.tagnum

    def __repr__(self):
        return '\n'.join([repr(self.__class__),"tags:",','.join(self.tags)])


    def iteritems(self):
        for tag in self.tags:
            yield tag,self.data[tag]

    def set_ylabel(self, ylabel, fontdict=None, labelpad=None, **kwargs):
        """
        """
        ylabdic = self._propagate(ylabel)
        for tag,ylabel in ylabdic.items():
            self.data[tag].set_ylabel(ylabel,
                                       fontdict=fontdict,
                                       labelpad=labelpad,
                                       **kwargs)
    def set_xlabel(self, xlabel, fontdict=None, labelpad=None, **kwargs):
        """
        """
        xlabdic = self._propagate(xlabel)
        for tag,xlabel in xlabdic.items():
            self.data[tag].set_xlabel(xlabel,
                                       fontdict=fontdict,
                                       labelpad=labelpad,
                                       **kwargs)

    def set_xlim(self,xlim,**kwargs):
        xlimdict = self._propagate(xlim,itearg=True)
        for tag,xlim in xlimdict.items():
            self.data[tag].set_xlim(xlim,**kwargs)

    def set_ylim(self,ylim,**kwargs):
        ylimdict = self._propagate(ylim,itearg=True)
        for tag,ylim in ylimdict.items():
            self.data[tag].set_ylim(ylim,**kwargs)

    def set_axis_bgcolor(self,color):
        colordict = self._propagate(color,itearg=False)
        for tag,c in colordict.items():
            self.data[tag].set_axis_bgcolor(c)

    def add_label(self,label=None,pos='ul',ftdic={'size':12},**kwargs):
        """
        cf g.Set_AxText
        """
        for tag,ax in self.iteritems():
            text = _replace_none_by_given(label,tag)
            g.Set_AxText(ax,text,pos=pos,ftdic=ftdic,**kwargs)

    def _propagate(self,arg,itearg=False):
        """
        Propagate the arg input to a (ordered)dict.

        when arg is a list with equal length of tags, the default tags
            sequence will be used.

        Parameters:
        -----------
        arg: the argument that's to be broadcasted.
        itearg: True if the arg itself is expected to be an iterable
        """
        if isinstance(arg,list):
            if itearg:
                if isinstance(arg[0],(tuple,list,np.ndarray)):
                    if len(arg) != self.tagnum:
                        raise ValueError("""list length expected to be {0}"""
                                         .format(self.tagnum))
                    else:
                        return OrderedDict(zip(self.tags,arg))
                else:
                    return dict(zip(self.tags,[arg]*self.tagnum))
            else:
                if len(arg) != self.tagnum:
                    raise ValueError("""list length expected to be {0}"""
                                     .format(self.tagnum))
                else:
                    return OrderedDict(zip(self.tags,arg))
        elif isinstance(arg,dict):
            return arg
        else:
            return dict(zip(self.tags,[arg]*self.tagnum))



    def call(self,funcname,*args,**kwargs):
        """
        Call a mat.axes.Axes method with identical args and kwargs for
            all tags.
        """
        for tag,ax in self.data.items():
            ax.__getattribute__(funcname)(*args,**kwargs)

    def apply(self,func,copy=False):
        """
        Apply function that applies on axes object.
        """
        if not copy:
            for tag,ax in self.data.items():
                func(ax)
        else:
            newaxs = [func(ax) for ax in self.data.values()]
            return LabelAxes(self.tags,newaxs)

    def duplicate_twinx(self):
        twinaxl = [host.twinx() for host in self.axl]
        return LabelAxes(tags=self.tags,axl=twinaxl)

    def build_icecore(self,num=3,keys=None,direction='vertical'):
        """
        Build an icecore like axes by num.

        Notes:
        ------
        1. The sequence is always from bottom to top.
        """
        newaxs = [g.Axes_Replace_by_IceCore(ax,num,direction=direction)
                    for ax in self.axl]
        newaxs = np.array(newaxs)
        if keys == None:
            keys = range(newaxs.shape[1])

        dic = OrderedDict()
        for i,key in enumerate(keys):
            dic[key] = LabelAxes(tags=self.tags,axl=newaxs[:,i])
        return dic

def apply_list_lax(laxlist,func,copy=False):
    """
    Apply function to list of LabelAxes object
    """
    if not copy:
        for lax in laxlist:
            lax.apply(func)
    else:
        return [lax.apply(func,True) for lax in laxlist]


class LabelAxes2D(object):
    """
    Initiate by a OrderedDict of LabelAxes object.
    """
    def __init__(self,laxdic):
        if not isinstance(laxdic,OrderedDict):
            raise TypeError("must be OrderedDict of LabelAxes objects")
        else:
            if not isinstance(laxdic.values()[0],LabelAxes):
                raise TypeError('dict values must be LabelAxes objects')
            else:
                self.child_lax = laxdic
                self.parent_tags = laxdic.keys()
                self.child_tags = laxdic.values()[0].tags


    def __repr__(self):
        return '\n'.join([repr(self.__class__),
                          "parent_tags:",','.join(self.parent_tags),
                          "child_tags:",','.join(self.child_tags)])

    def __getitem__(self,key):
        return self.child_lax[key]

    def iteritems(self):
        for ptag in self.parent_tags:
            yield ptag,self.child_lax[ptag]

    def set_xlim(self,xlim,**kwargs):
        for ptag,plax in self.child_lax.items():
            plax.set_xlim(xlim,**kwargs)

    def set_ylim(self,ylim,**kwargs):
        for ptag,plax in self.child_lax.items():
            plax.set_ylim(ylim,**kwargs)

    def set_axis_bgcolor(self,color):
        colordict = _propagate(self.parent_tags,color,itearg=False)
        for ptag in self.parent_tags:
            self.child_lax[ptag].set_axis_bgcolor(colordict[ptag])

    def add_parent_label(self,pos='ouc',ftdic={'size':12},**kwargs):
        for ptag,lax in self.iteritems():
            lax.add_label(label=ptag,pos=pos,ftdic=ftdic,**kwargs)

    def add_child_label(self,ptag=None,pos='ouc',color='m',ftdic={'size':12},**kwargs):
        ptag = _replace_none_by_given(ptag,self.parent_tags[-1])
        self.child_lax[ptag].add_label(pos=pos,ftdic=ftdic,color=color,**kwargs)





