from osgeo import ogr
import matplotlib as mat


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
        return linelist

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
        polygon_list = []
        for i in range(polygon_num):
            polygon = mpolygon.GetGeometryRef(i)
            linelist = get_linelist_from_polygon(polygon)
            polygon_list.append(linelist)
        return polygon_list


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

def get_geometry_type_from_feature(feature):
    georef = feature.GetGeometryRef()
    geometry_name = georef.GetGeometryName()
    return geometry_name


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



