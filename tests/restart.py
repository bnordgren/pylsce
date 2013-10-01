'''
This file contains some code to support the analysis of restart files.
Restart files contain different variables with different dimension names.
'''

import netCDF4 as nc

def compare_summary(first, second, threshold=1.e-6) :
    retval = True
    differences = []
    for pft in range(13) : 
        difference = first[:,pft].sum() - second[:,pft].sum()
        differences.append(difference)
        retval = retval and (abs(difference) <= threshold)
    return (retval, differences)

def load_sechiba(year, variable='veget_max') : 
    sechiba = nc.Dataset('sechiba_restart_%d.nc' % year)
    return sechiba.variables[variable]

def load_stomate(year, variable='PFTpresent') : 
    stomate = nc.Dataset('stomate_restart_%d.nc' % year)
    return stomate.variables[variable]

def sumByPFT(data):
    totals = []
    for pft in range(13) : 
        totals.append(data[:,pft].sum())
    return totals

