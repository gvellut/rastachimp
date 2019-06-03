# RastaChimp

RastaChimp is a utility library to perform topological processing of polygons vectorized from rasters, as returned by the [rasterio.features.shapes function](https://mapbox.github.io/rasterio/topics/features.html):

- Simplification, using Douglas-Peucker: [Shapely](https://toblerity.org/shapely/manual.html) (and GEOS) can be used to simplify polygons but will create gaps and overlaps between polygons that share a boundary (even with preserve_topology=true)
- Smoothing, using the Chaikin method
- Densify

RastaChimp decomposes the polygons into edges (using Shapely), applies an operation to them then rebuilds transformed polygons. It works well as long as the output of rasterio.features.shapes is in pixel coordinates (default) or the transform to world coordinates does not have a rotation component, which is usually the case. Because of limited floating precision, it may not work well in the more general case, although it could if the polygons are noded (for example, by using the [bpol tool of v.clean in GRASS](https://grass.osgeo.org/grass73/manuals/v.clean.html)).

### Install

`pip install rastachimp`

### Usage

The function to perform simplification is ```simplify_dp```. As inputs, it takes the tuples (geometry, value) as output from rasterio.features.shapes and the tolerance (in the same unit as the input geometries). It returns an iterator of tuples (geometry, value), just like rasterio.features.shapes.

```python
from rasterio import features
from rastachimp import as_shapely, simplify_dp

# 2D array with classification values
image = np.array(...., dtype=np.uint8)

# no transform => in pixel coordinates
shapes = features.shapes(image)

# convert to Shapely geometry (features.shapes returns a GeoJSON dict)
shapes = as_shapely(shapes)

tolerance = 2
# simpl is an interator of tuples (geometry, value)
simpl = simplify_dp(shapes, tolerance)

```

### Issues

Report issues at https://github.com/gvellut/rastachimp/issues

### TODO

- Better detection and handling of bad inputs
- Other generalization algorithms (Visvalingam)
- Other smoothing algorithm (Bezier Curves)

