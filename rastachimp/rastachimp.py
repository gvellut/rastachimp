from collections import defaultdict
import logging
from operator import itemgetter

from numba import jit
import numpy as np
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon
from shapely.ops import linemerge, unary_union
from shapely.prepared import prep

from .utils import build_spatial_index, flatten, set_geom, take

logger = logging.getLogger(__package__)


def apply_topological_op(features, topo_op, **params):
    """
    Apply topo_op to the graph formed by the polygon features.
    The input features are an array of tuples (<geom>,<data>). The geoms can be a mix
    of polygons and multipolygons
    The lines that belong to multiple polygons are merged before being passed
    to topo_op
    """
    geoms, data = zip(*features)
    multi, polygons = _decompose_multi(geoms)
    # edges is a single MultiLineString
    edges = _decompose(polygons)
    faces = _build_topology(polygons, edges)

    # topoop_edges must be mapped 1-1 with the output edges
    topoop_edges, faces, edges = topo_op(faces, edges, **params)

    # topoop_polygons are mapped 1-1 with the polygons
    topoop_polygons = _build_geometry(faces, topoop_edges)
    topoop_multi = _recompose_multi(multi, topoop_polygons)
    topoop_features = zip(topoop_multi, data)
    return topoop_features


def smooth_chaikin(features, iterations=2, keep_border=False):
    """
    Smooth the edges using the Chaikin method
    This op doesn't offer the same guarantees as simplify_dp : It is possible for
    edges to cross after this op is applied
    """
    sm_features = apply_topological_op(
        features, _smooth_chaikin_edges, iterations=iterations, keep_border=keep_border
    )
    return sm_features


def _smooth_chaikin_edges(
    faces, edges: MultiLineString, iterations, keep_border
) -> MultiLineString:

    if keep_border:
        border_index = _detect_border(faces)

    xy_edge_index = _index_edges_by_xy(edges)

    sm_edges = []
    for edge_index, edge in enumerate(edges):
        if keep_border and edge_index in border_index:
            # nothing to do
            sm_edges.append(edge)
            continue

        # the 2nd condition is so that there is no smoothing when
        # the edge touches another edge at a single point.
        # By construction, this point is always the first (and last)
        # of coords
        smooth_around = edge.is_ring and len(xy_edge_index[edge.coords[0]]) == 1
        if smooth_around:
            # so that the smoothing also happens between the last and first
            # segment
            coords = np.array(list(edge.coords) + [edge.coords[1]])
        else:
            coords = np.array(edge.coords)

        # cf https://stackoverflow.com/a/47255374
        for _ in range(iterations):
            L = coords.repeat(2, axis=0)
            R = np.empty_like(L)
            R[0] = L[0]
            R[2::2] = L[1:-1:2]
            R[1:-1:2] = L[2::2]
            R[-1] = L[-1]
            coords = L * 0.75 + R * 0.25

        coords = coords.tolist()
        if smooth_around:
            # points have been added in the loop
            index_close = 2 ** iterations
            coords = coords[index_close - 1 : -index_close]
        sm_edges.append(LineString(coords))

    # faces and edges have not been modified
    return MultiLineString(sm_edges), faces, edges


def densify(features, max_distance: float = None, n: int = 1, keep_border=False):
    """
    Add points to the boundaries of the polygons
    Can be useful to control the smoothing performed by the Chaikin method
    """
    ds_features = apply_topological_op(
        features,
        _densify_edges,
        max_distance=max_distance,
        n=n,
        keep_border=keep_border,
    )
    return ds_features


def _densify_edges(
    faces, edges: MultiLineString, max_distance, n, keep_border
) -> MultiLineString:

    if keep_border:
        border_index = _detect_border(faces)

    ds_edges = []
    for edge_index, edge in enumerate(edges):
        if keep_border and edge_index in border_index:
            # nothing to do
            ds_edges.append(edge)
            continue

        coords = np.array(edge.coords)
        if max_distance:
            dists = _point_distance(coords)
            ns = np.ceil(dists / max_distance).astype(np.int64)
        else:
            # np points => np-1 segments
            ns = np.full(len(coords) - 1, n + 1)

        # + 1 to contain the last point
        n_ds_points = np.sum(ns) + 1
        ds_coords = np.zeros((n_ds_points, 2))
        _densify_edge(coords, ns, ds_coords)
        ds_edges.append(LineString(zip(ds_coords[:, 0], ds_coords[:, 1])))

    return MultiLineString(ds_edges), faces, edges


@jit(nopython=True)
def _densify_edge(input, ns, output):
    curr_j = 0
    for i in range(len(input) - 1):
        n = ns[i]
        start = input[i]
        stop = input[i + 1]
        diff = stop - start
        for j in range(n):
            output[curr_j + j] = start + diff * j / n
        curr_j += n
    # the last point is not processed in the loop: stays the same in the output
    output[-1] = input[-1]


def simplify_dp(features, distance, keep_border=False):
    """
    Simplify polygon features using the Douglas-Peucker method.
    This op simplifies edges shared by multiple polygons in the same way. It will
    also prevent edges from crossing each other.
    """

    simpl_features = apply_topological_op(
        features, _simplify_dp_edges, distance=distance, keep_border=keep_border
    )

    # remove degenerate polygons (no area)
    f_simpl_features = []
    for f in simpl_features:
        # if polgygon is degenerate: This will make it empty
        # if mulitp: It will remove degenerate member polygons
        geom = f[0].buffer(0)
        if geom.is_empty:
            # degenerate
            continue
        f_simpl_features.append(set_geom(f, geom))

    return f_simpl_features


def _simplify_dp_edges(faces, edges: MultiLineString, distance, keep_border):
    if keep_border:
        border_index = _detect_border(faces)
        # simplify from Shapely does not change polylines composed of only 2 points
        # so the edges at the border are cut into 2 point polylines
        faces, edges = _segmentize_border(faces, edges, border_index)

    # e.simplify does DP (from Shapely / GEOS)
    return edges.simplify(distance, True), faces, edges


def _decompose_multi(geoms):
    multi = []
    polygons = []
    for geom in geoms:
        if geom.type == "MultiPolygon":
            polygons.extend(geom.geoms)
            multi.append(
                ("MultiPolygon", list(range(len(polygons) - len(geom), len(polygons))))
            )
        else:
            # Polygon
            polygons.append(geom)
            multi.append(("Polygon", [len(polygons) - 1]))
    return multi, polygons


def _recompose_multi(multi, polygons):
    rec_multi = []
    for gtype, poly_indices in multi:
        if gtype == "MultiPolygon":
            rec_multi.append(MultiPolygon(take(polygons, poly_indices)))
        else:
            # a single geoemtry
            rec_multi.append(polygons[poly_indices[0]])
    return rec_multi


def _decompose(polygons):
    lines = flatten((_polygon_contour(p) for p in polygons))
    union = unary_union(lines)
    if isinstance(union, LineString):
        edges = MultiLineString([union])
    else:
        edges = linemerge(union)
    return edges


def _build_topology(polygons, edges):
    """
    Decompose polygons in an ordered sequence of edges
    """
    edge_si = build_spatial_index(edges)
    faces = []
    for polygon in polygons:
        rings = [polygon.exterior] + list(polygon.interiors)
        preps = map(prep, rings)
        face = []
        for ring, pring in zip(rings, preps):
            indexes = list(edge_si.intersection(ring.bounds))
            indexes = list(filter(lambda i: pring.contains(edges[i]), indexes))
            # sequence of oriented edges for current ring
            topo_ring = _build_topo_ring(
                list(map(lambda i: edges[i], indexes)), indexes
            )
            face.append(topo_ring)
        faces.append(face)
    return faces


def _build_topo_ring(edges, edge_global_indexes):
    # Second element of tuple indicates if the topo_ring has the same orientation
    # as the edge
    # init : arbitrary direction for first edge
    result, node_xy = [(0, True)], edges[0].coords[-1]

    indexes = set(range(1, len(edges)))
    while len(indexes) > 0:
        # tries to find edge that starts or ends with node_xy
        next_ = None
        candidates_next = []
        for idx in indexes:
            if node_xy in (edges[idx].coords[0], edges[idx].coords[-1]):
                candidates_next.append(idx)

        if len(candidates_next) == 1:
            next_ = candidates_next[0]
        elif len(candidates_next) == 2:
            # can happen if the ring boundary touches itself at a single point
            # (not valid but output by rasterio / GDAL in polygonize)
            # one of the edges forms a closed ring (but not considered a hole)
            # we take that edge as next
            logger.warning("Invalid polygon: Self-touching contour")
            if (
                edges[candidates_next[0]].coords[0]
                == edges[candidates_next[0]].coords[-1]
            ):
                next_ = candidates_next[0]
            elif (
                edges[candidates_next[1]].coords[0]
                == edges[candidates_next[1]].coords[-1]
            ):
                next_ = candidates_next[1]
        # else : should not happen (invalid polygon)

        # if the input polygons are valid and not degenerate
        # it should be always possible to find a next_
        assert next_ is not None

        # TODO have more explicit error messages in case of invalid inputs

        if node_xy == edges[next_].coords[0]:
            result.append((next_, True))
            node_xy = edges[next_].coords[-1]
        else:
            result.append((next_, False))
            node_xy = edges[next_].coords[0]

        indexes.remove(next_)

    return list(map(lambda r: (edge_global_indexes[r[0]], r[1]), result))


def _detect_border(faces):
    edge_counter = defaultdict(int)
    for if_, face in enumerate(faces):
        for ir, ring in enumerate(face):
            for ie, (edge_index, is_same_order) in enumerate(ring):
                edge_counter[edge_index] += 1
    return set([edge_index for edge_index, c in edge_counter.items() if c == 1])


def _segmentize_border(faces, edges, border_index):
    """
    Cuts edges at the border (part of a single polygon)
    into its constituent segments
    """

    # to simple linestrings
    edges = list(edges.geoms)

    borders = []
    for if_, face in enumerate(faces):
        for ir, ring in enumerate(face):
            for ie, (edge_index, is_same_order) in enumerate(ring):
                if edge_index in border_index:
                    borders.append((if_, ir, ie, edge_index, is_same_order))

    # the edges at the border are part of only a sing ring
    # however, a ring can have multiple edges at the border
    # we sort so that the greatest "ie" is processed first
    # => when splice is performed, it won't disturb the next ones
    borders = sorted(borders, reverse=True, key=itemgetter(2))
    for border in borders:
        if_, ir, ie, edge_index, order = border
        coords = edges[edge_index].coords
        segments = []
        for i in range(len(coords) - 1):
            # same order as the origin edge
            segments.append(LineString([coords[i], coords[i + 1]]))
        if not order:
            segments = list(reversed(segments))
        num_indices_added = len(segments) - 1
        indices_added = list(
            zip(
                range(len(edges), len(edges) + num_indices_added),
                [order] * num_indices_added,
            )
        )
        # so that there are no holes and no need to change the existing indices
        # first output segment takes the index of its edge
        edges[edge_index] = segments[0]
        # add new segments at the end
        edges.extend(segments[1:])
        # the index of the first segment (at "ie") stays the same (edge_index)
        # the others are spliced after in the index list of the ring
        splice_index = ie + 1
        faces[if_][ir][splice_index:splice_index] = indices_added

    return faces, MultiLineString(edges)


def _index_edges_by_xy(edges):
    rindex = defaultdict(set)
    for i, edge in enumerate(edges):
        rindex[edge.coords[0]].add(i)
        rindex[edge.coords[-1]].add(i)
    return rindex


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
            # revert
            coords.extend(edges[edge_index].coords[:0:-1])
    # close the contour
    coords.append(coords[0])
    return coords


def _polygon_contour(polygon):
    return [LineString(tuple(polygon.exterior.coords))] + [
        LineString(tuple(geom.coords)) for geom in polygon.interiors
    ]


def _point_distance(coords):
    # the coordinates are each on a row with 2 elements
    diff = np.diff(coords, axis=0)
    dists = np.sqrt(diff[:, 0] ** 2 + diff[:, 1] ** 2)
    return dists
