import fiona
import numpy as np
from pytest import approx
from rasterio import features
from shapely.geometry import LineString, mapping

from rastachimp import as_shapely, simplify_dp, smooth_chaikin
from rastachimp.rastachimp import _densify_edges, _point_distance


def test_simplify_basic():
    shapes = _sample_data_basic()
    # convert to Shapely geometry
    simpl = list(simplify_dp(as_shapely(shapes), 2, True))

    assert 6 == len(simpl)

    count = {2: 0, 1: 0, 5: 0}
    for _, value in simpl:
        count[value] += 1
    assert count[2] == 2
    assert count[1] == 3
    assert count[5] == 1

    _to_shp(simpl, "simpl_dp.shp")


def test_smooth_basic():
    shapes = _sample_data_basic()

    simpl = list(simplify_dp(as_shapely(shapes), 0.9, True))
    smooth = list(smooth_chaikin(simpl, 5, True))

    assert 6 == len(smooth)

    count = {2: 0, 1: 0, 5: 0}
    for _, value in smooth:
        count[value] += 1
    assert count[2] == 2
    assert count[1] == 3
    assert count[5] == 1

    _to_shp(smooth, "smooth_chaikin.shp")


def test_densify():
    # 12 points
    coords = np.arange(24).reshape(12, 2)
    edge = LineString(coords)

    # subdivide once (n=1)
    ds_edges, _, _ = _densify_edges(None, [edge], None, n=1, keep_border=False)
    ds_coords = ds_edges[0].coords

    # 2 new segments for each input
    assert len(ds_coords) == 23
    # check that coordinates of input poins are still the same
    for i in range(0, 23, 2):
        assert tuple(ds_coords[i]) == tuple(coords[i // 2])

    # check new segment is colinear with the original input
    np_coords = np.array(ds_coords)
    for i in range(1, 23, 2):
        # points at i are new, those at i+1 and i-1 are same as
        # input
        x1, y1 = np_coords[i] - np_coords[i - 1]
        x2, y2 = np_coords[i + 1] - np_coords[i - 1]
        assert x1 * y2 - x2 * y1 == approx(0)

    # all segments have the same length ~2.82
    ds_edges, _, _ = _densify_edges(
        None, [edge], max_distance=1, n=None, keep_border=False
    )
    ds_coords = ds_edges[0].coords

    # for max_distance=1 => 3 new segments for each input
    assert len(ds_coords) == 34
    dists = _point_distance(ds_coords)
    for i in range(len(dists)):
        assert dists[i] < 1


def _sample_data_basic():
    image = np.array(
        [
            [2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1],
            [2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1],
            [2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1],
            [2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1],
            [5, 5, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1],
            [5, 5, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1],
            [1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1],
            [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1],
            [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        ],
        dtype=np.uint8,
    )

    shapes = features.shapes(image)
    return shapes


def _to_shp(transf, dst):
    meta = {
        "schema": {"geometry": "Polygon", "properties": {"VALUE": "int:10"}},
        "crs": "EPSG:3857",
    }
    with fiona.open(dst, "w", driver="Shapefile", **meta) as sink:
        for s in transf:
            f = {"geometry": mapping(s[0]), "properties": {"VALUE": int(s[1])}}
            sink.write(f)
