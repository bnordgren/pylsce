#!/usr/bin/python

import g


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
