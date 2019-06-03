from rtree import Rtree
from shapely.geometry import shape


def set_geom(f, geom):
    return (geom, f[1])


def as_shapely(fs):
    return (set_geom(f, shape(f[0])) for f in fs)


def extract_geoms(fs):
    return [f[0] for f in fs]


def build_spatial_index(geometries):
    spatial_index = Rtree()
    for g, geom in enumerate(geometries):
        spatial_index.add(g, geom.bounds)
    return spatial_index


def take(array, indices):
    out = []
    for i in indices:
        out.append(array[i])
    return out


def flatten(l):
    return [item for subl in l for item in subl]
