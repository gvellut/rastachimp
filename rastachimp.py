# -*- coding: utf-8 -*-

from rtree import Rtree
from shapely.geometry import Polygon, LineString, shape, mapping
from shapely.prepared import prep
from shapely.ops import unary_union, linemerge
        
def simplify(polygons_f_with_value, distance):
    polygons_f, values = zip(*polygons_f_with_value)
    # fiona polygons to shapely polygons
    polygons = list(map(shape, polygons_f))
    lines = _flatten([_polygon_contour(p) for p in polygons ])
    union = unary_union(lines)
    edges = linemerge(union)
    # simpl_edges are mapped 1-1 with original edges
    simpl_edges = edges.simplify(distance, True)
    faces = _build_topology(polygons, edges)
    simpl_polygons = _build_geometry(faces, simpl_edges)
    # add back original values and suppr degenerate polygons
    simpl_polygons_with_values = filter(lambda p : not _is_degenerate(p[0]), zip(simpl_polygons, values))
    # shapely polygons to fiona polygons
    return map(lambda x : (mapping(x[0]), x[1]), simpl_polygons_with_values)
    
def _build_spatial_index(geometries):
    spatial_index = Rtree()
    for g,geom in enumerate(geometries):
        spatial_index.add(g,geom.bounds)
    return spatial_index
    
# builds topology based on polygons and edges
def _build_topology(polygons, edges):
    edge_si = _build_spatial_index(edges)
    faces = []
    for polygon in polygons:
        rings = [polygon.exterior] + list(polygon.interiors)
        preps = map(prep, rings)
        face = []
        for (ring,pring) in zip(rings,preps):
            indexes  = list(edge_si.intersection(ring.bounds))
            indexes  = list(filter(lambda i: pring.contains(edges[i]),indexes))
            # sequence of oriented arcs for current ring
            topo_ring = _build_topo_ring(list(map(lambda i: edges[i],indexes)),indexes)
            face.append(topo_ring)
        faces.append(face)
    return faces
    
def _build_topo_ring(edges, edge_global_indexes):
    # suppose edges are sufficient to build a ring
    
    # second element of result tuples indicates if the topo_ring has the
    # same orientation as edge or not 
    # init : arbitrary direction for first edge  
    result, node_xy = [(0,True)], edges[0].coords[-1]
    
    indexes = set(range(1, len(edges)))
    while len(indexes) > 0:
        # search for edge that starts or ends at node_xy
        next_ = None
        for idx in indexes:
            if node_xy in (edges[idx].coords[0],edges[idx].coords[-1]):
                next_ = idx
                break
    
        if next_ is None:
            raise Exception('_build_topo_ring: cannot find the next edge')
    
        if node_xy == edges[next_].coords[0]:
            result.append((next_,True))
            node_xy = edges[next_].coords[-1]
        else:
            result.append((next_,False))
            node_xy = edges[next_].coords[0]
    
        indexes.remove(next_)
    
    return list(map(lambda r: (edge_global_indexes[r[0]],r[1]), result))
    
def _build_geometry(faces, edges):
    polygons = []
    for face in faces:
        rings = []
        for topo_ring in face:
            ring = _build_geom_ring(topo_ring, edges)
            rings.append(ring)
        polygons.append(Polygon(rings[0], rings[1:]))
    return polygons
    
def _build_geom_ring(ring, edges):
    coords = []
    for edge_index, is_same_order in ring:
        if is_same_order:
            coords.extend(edges[edge_index].coords[:-1])
        else:
            # reverse
            coords.extend(edges[edge_index].coords[:0:-1])
    # close
    coords.append(coords[0])
    return coords    
       
def _polygon_contour(polygon):
    return [ LineString(tuple(polygon.exterior.coords)) ] + [ LineString(tuple(geom.coords)) for geom in polygon.interiors ]
    
def _is_degenerate(polygon):
    return polygon.area == 0

def _flatten(l):
    return [item for subl in l for item in subl]
 
