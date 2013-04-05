from osgeo import ogr
import matplotlib as mat
import pandas as pa
import numpy as np
import Pdata

def _check_df_inner_out_ring_validity(df):
    if not isinstance(df,pa.DataFrame):
        raise TypeError("first argument must be DataFrame")
    else:
        if 'inner_ring' not in df.columns or 'out_ring' not in df.columns:
            raise NameError("inner_ring/out_ring not dataframe column")
        else:
            pass


def get_layer_attribute_table(layer,feature_range=None):
    if feature_range == None:
        select_range = range(layer.GetFeatureCount())
    else:
        select_range = feature_range
    data_list = []
    index_list = []

    for i in select_range:
        feature = layer.GetFeature(i)
        data_list.append(feature.items())
        index_list.append(i)
    return pa.DataFrame(data_list,index=index_list)

def get_line_from_linearring(ring):
    """
    Get the list of vertices in the linearring.

    Returns:
    --------
    a list of tuples which represent the vertices.
    """
    geometry_name = ring.GetGeometryName()
    if geometry_name != 'LINEARRING':
        raise TypeError("The type is {0}".format(geometry_name))
    else:
        num_point = ring.GetPointCount()
        d = [ring.GetPoint(i) for i in range(num_point)]
        line = [(a,b) for (a,b,c) in d]
        return line


def get_linelist_from_polygon(polygon):
    """
    Get a "polygon" from the polygon object.

    Returns:
    --------
    A list of sublist, each sublist is a list of vertices from
        a linearring object.
    """
    geometry_name = polygon.GetGeometryName()
    if geometry_name != 'POLYGON':
        raise TypeError("the type is {0}".format(geometry_name))
    else:
        geocount = polygon.GetGeometryCount()
        linelist = []
        for i in range(geocount):
            ring = polygon.GetGeometryRef(i)
            line = get_line_from_linearring(ring)
            if line != []:
                linelist.append(line)
    if len(linelist) == 1:
        return (linelist,None)
    else:
        return (linelist[0:1],linelist[1:])

def get_lines_from_multipolygon(mpolygon):
    """
    Get lines from MultiPolygon object.

    Returns:
    --------
    """
    geometry_name = mpolygon.GetGeometryName()
    polygon_num = mpolygon.GetGeometryCount()
    if geometry_name != 'MULTIPOLYGON':
        raise TypeError("the type is {0}".format(geometry_name))
    else:
        out_ring_list,inner_ring_list = [],[]
        for i in range(polygon_num):
            polygon = mpolygon.GetGeometryRef(i)
            (out_ring,inner_ring) = get_linelist_from_polygon(polygon)
            for subring in out_ring:
                out_ring_list.append(subring)
            if inner_ring != None:
                if len(inner_ring) == 1:
                    inner_ring_list.append(inner_ring[0])
                else:
                    for subring in inner_ring:
                        inner_ring_list.append(subring)
            else:
                pass
    if inner_ring_list == []:
        return (out_ring_list,None)
    else:
        return (out_ring_list,inner_ring_list)


def get_geometry_from_feature(feature):
    """
    Get geometry from feature.
    """
    georef = feature.GetGeometryRef()
    geometry_name = georef.GetGeometryName()
    if geometry_name == 'POLYGON':
        return get_linelist_from_polygon(georef)
    elif geometry_name == 'MULTIPOLYGON':
        return get_lines_from_multipolygon(georef)
    else:
        raise TypeError("input feature type is {0}".format(geometry_name))


def transform_layer_geometry_to_ring_dataframe(layer,feature_range=None):
    data_list = []
    index_list = []

    if feature_range == None:
        select_range = range(layer.GetFeatureCount())
    else:
        select_range = feature_range

    for i in select_range:
        feature = layer.GetFeature(i)
        out_ring_list,inner_ring_list = get_geometry_from_feature(feature)
        data_list.append({'out_ring':out_ring_list, 'inner_ring':inner_ring_list})
        index_list.append(i)
    return pa.DataFrame(data_list,index=index_list)

def dataframe_of_ring_change_projection(df,m):
    _check_df_inner_out_ring_validity(df)
    dfnew = df.copy()
    for name in ['inner_ring','out_ring']:
        for i in dfnew.index:
            if dfnew[name][i] == None:
                pass
            else:
                ddt = dfnew[name][i]
                dfnew[name][i] = map(lambda templist:map(lambda x:m(*x),templist),ddt)
    return dfnew

def group_dataframe_of_ring(df,groupby):
    """
    group the inner_ring,out_ring dataframe by key.
    """
    _check_df_inner_out_ring_validity(df)
    grp = df.groupby(groupby)

    def merge_list(inlist):
        outlist = []
        for first_level_list in inlist:
            if first_level_list == None:
                pass
            else:
                for sublist in first_level_list:
                    outlist.append(sublist)
        return outlist

    dfdic = {}
    for name in ['inner_ring','out_ring']:
        dfdic[name] = grp[name].apply(merge_list)
    return pa.DataFrame(dfdic)


def get_geometry_type_from_feature(feature):
    georef = feature.GetGeometryRef()
    geometry_name = georef.GetGeometryName()
    return geometry_name

def get_geometry_count_from_feature(feature):
    georef = feature.GetGeometryRef()
    geometry_count = georef.GetGeometryCount()
    return geometry_count

def Add_Linearring_to_Axes(ax,ring,facecolor='0.7',edgecolor='k',
                           transfunc=None,
                           **kwargs):
    """
    Add Linearring to Axes.

    Parameters:
    -----------
    transfunc: functions used for spatial transformation, they should receive
        tuple as parameter and return tuple.
    """
    if transfunc == None:
        ringnew = ring
    else:
        ringnew = [transfunc(t) for t in ring]

    poly = mat.collections.PolyCollection([ringnew],
                                          facecolor=facecolor,
                                          edgecolor=edgecolor,
                                          **kwargs)
    ax.add_collection(poly)

def Add_Polygon_to_Axes(ax,list_of_rings,facecolor='0.7',edgecolor='k',
                        inner_ring_facecolor='w',inner_ring_edgecolor='k',
                        inner_ring_kwargs={},
                        transfunc=None,
                        **kwargs):
    if len(list_of_rings) == 0:
        raise ValueError("input list_of_rings has length 0")
    elif len(list_of_rings) == 1:
        Add_Linearring_to_Axes(ax,list_of_rings[0],facecolor=facecolor,
                               edgecolor=edgecolor,
                               transfunc=transfunc,**kwargs)
    else:
        Add_Linearring_to_Axes(ax,list_of_rings[0],facecolor=facecolor,
                               edgecolor=edgecolor,
                               transfunc=transfunc,
                               **kwargs)
        for ring in list_of_rings[1:]:
            Add_Linearring_to_Axes(ax,ring,facecolor=inner_ring_facecolor,
                                   edgecolor=inner_ring_edgecolor,
                                   transfunc=transfunc,
                                   **inner_ring_kwargs)

def Add_MultiPolygon_to_Axes(ax,list_of_polygon,
                             facecolor='0.7',
                             edgecolor='k',
                             inner_ring_facecolor='w',
                             inner_ring_edgecolor='k',
                             inner_ring_kwargs={},
                             transfunc=None,
                             **kwargs):
    if len(list_of_polygon) == 0:
        raise ValueError("input list_of_polygon has length 0")
    else:
        for list_of_rings in list_of_polygon:
            Add_Polygon_to_Axes(ax,list_of_rings,
                                facecolor=facecolor,
                                edgecolor=edgecolor,
                                inner_ring_facecolor=inner_ring_facecolor,
                                inner_ring_edgecolor=inner_ring_edgecolor,
                                inner_ring_kwargs={},
                                transfunc=transfunc,
                                **kwargs)

def Add_Feature_to_Axes_Polygon(ax,feature,transfunc=None):
    geotype = get_geometry_type_from_feature(feature)
    geometry_vertices = get_geometry_from_feature(feature)
    if geotype == 'POLYGON':
        Add_Polygon_to_Axes(ax,geometry_vertices,inner_ring_facecolor='w',
                            transfunc=transfunc)
    elif geotype == 'MULTIPOLYGON':
        Add_MultiPolygon_to_Axes(ax,geometry_vertices,inner_ring_facecolor='w',
                                 transfunc=transfunc)
    else:
        raise ValueError("geometry type not polygon!")



def dataframe_build_geoindex_from_lat_lon(df,lat_name='lat',
                                          lon_name='lon',
                                          lat=None,lon=None):
    """
    Build a geoindex column for the dataframe "df", by check each
        latitude/longitude pairs (lat_name/lon_name) falling in which
        grid cell of the grid as specified by the vectors of lat/lon.
        The latitude/longitude pairs falling outside the grid will
        have geoindex values as np.nan.

    Parameters:
    -----------
    df: input dataframe.
    lat_name/lon_name: the latitude/longitude field name of the dataframe.
    lat/lon: the latitude/longitude vectors used to compose the grid.
    """
    df['geoindex'] = [(None,None)]*len(df.index)
    for i in df.index:
        vlat = df[lat_name][i]
        vlon = df[lon_name][i]
        try:
            df['geoindex'][i] = gnc.find_index_by_point(lat_range,
                                                        lon_range,(vlat,vlon))
        except ValueError:
            df['geoindex'][i] = np.nan
    return df

def mdata_by_geoindex_dataframe(df,shape=(360,720),mask=None):
    """
    Transfer the geoindexed dataframe into Pdata.Mdata for plotting
        or writing out ot NetCDF file.

    Parameters:
    ----------
    shape: the shape of array to be constructed, limited to 2D array.
    mask: the mask that's to be applied.

    Notes:
    ------
    1. the df.index must be tuples.
    """
    ydic = {}
    for name in df.columns.tolist():
        data = np.ones(shape)*np.nan
        for index,value in df[name].iterkv():
            if not isinstance(index,tuple):
                raise TypeError("index {0} not tuple".format(index))
            else:
                data[index]=value
        if mask != None:
            data = np.ma.masked_array(data,mask=mask)
        ydic[name] = data
    return Pdata.Mdata.from_dict_of_array(ydic)


