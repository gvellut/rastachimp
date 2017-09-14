# RastaChimp

RastaChimp is a simple utility library to perform topological simplification of polygons vectorized from rasters, as returned by the [rasterio.features.shapes function](https://mapbox.github.io/rasterio/topics/features.html). 

[Shapely](https://toblerity.org/shapely/manual.html) (and GEOS) can be used to simplify polygons but will create gaps and overlaps between   polygons that share a boundary (even with preserve_topology=true). So instead, RastaChimp simplifies the edges that compose the polygons (still using Shapely) and rebuilds simplified polygons from them. It works well as long as the output of rasterio.features.shapes is in pixel coordinates (default) or the transform to world coordinates does not have a rotation component, which is usually the case. Because of limited floating precision, it may not work well in the more general case, although it could if the polygons are noded (for example, by using the [bpol tool of v.clean in GRASS](https://grass.osgeo.org/grass73/manuals/v.clean.html)).

### Install

- Copy the rastachimp.py file into your project
- Install recent versions of:
    - [Shapely](https://pypi.python.org/pypi/Shapely)
    - [Rtree](https://pypi.python.org/pypi/Rtree/)
    - [Rasterio](https://pypi.python.org/pypi/rasterio) is not strictly needed, as long as the input to the rastachimp.simplify function is in the same shape as what is returned from rasterio.features.shapes

### Usage

The only function to use from the package is ```simplify```. As inputs, it takes the tuples (geometry, value) as output from rasterio.features.shapes and the tolerance (in the same unit as the input geometries). It returns an iterator of tuples (geometry, value), jsut like rasterio.features.shapes.

```python
from rasterio import features
from rastachimp import simplify

# 2D array with classification values
image = np.array(...., dtype=np.uint8)

# no transform => in pixel coordinates
shapes = features.shapes(image)

tolerance = 2
# simpl is an interator of tuples (geometry, value)
simpl = simplify(shapes, tolerance)

```
### Issues

Report issues at https://github.com/gvellut/rastachimp/issues

### TODO

- Provide a pip-able install package

