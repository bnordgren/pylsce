#!/usr/bin/python

from collections import OrderedDict
import pandas as pa
import matplotlib as mat
import numpy as np

class LabelAxes(object):
    def __init__(self,tags=None,axl=None):
        if not isinstance(tags,(list,tuple)):
            raise TypeError("tags could only be list or tuple")
        else:
            if not isinstance(axl[0],mat.axes.Axes):
                raise TypeError('value must be mat.axes.Axes type!')
            else:
                self.data = OrderedDict(zip(tags,axl))
                self.tags = tags
                self.axl = list(axl)
                self.tagnum = len(self.tags)


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
                        return [self.tags[index] for index in key]
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

    def __repr__(self):
        return '\n'.join([repr(self.__class__),"tags:",','.join(self.tags)])



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

    def _propagate(self,arg,itearg=False):
        """
        Propagate the arg input to a (ordered)dict.

        when arg is a list with equal length of tags, the default tags
            sequence will be used.

        Parameters:
        -----------
        arg: the argument that's to be broadcasted.
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

    def apply(self,func):
        """
        Apply function that applies on axes object.
        """
        for tag,ax in self.data.items():
            func(ax)



