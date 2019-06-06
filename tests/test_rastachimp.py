import numpy as np
from pytest import approx
from rasterio import features
from shapely.geometry import LineString

from rastachimp import as_shapely, simplify_dp
from rastachimp.rastachimp import _densify_edges, _point_distance


def test_simplify_basic():
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
    # convert to Shapely geometry
    simpl = list(simplify_dp(as_shapely(shapes), 2.0, True))

    assert 6 == len(simpl)

    count = {2: 0, 1: 0, 5: 0}
    for _, value in simpl:
        count[value] += 1
    assert count[2] == 2
    assert count[1] == 3
    assert count[5] == 1


def test_densify():
    # 12 points
    coords = np.arange(24).reshape(12, 2)
    edge = LineString(coords)

    # subdivise 1 fois (n=1)
    ds_edges, _, _ = _densify_edges(None, [edge], None, n=1, keep_border=False)
    ds_coords = ds_edges[0].coords

    # 2 nouveaux segments par ancien
    assert len(ds_coords) == 23
    # vérifier que les coordonnées des anciens points
    # sont toujours les mêmes
    for i in range(0, 23, 2):
        assert tuple(ds_coords[i]) == tuple(coords[i // 2])

    # vérifier la colinéarité d'un nouveau segment
    # par rapport au segment dont il est issu
    np_coords = np.array(ds_coords)
    for i in range(1, 23, 2):
        # les points à i sont nouveaux, ceux à i+1 et i-1 sont les
        # originaux
        x1, y1 = np_coords[i] - np_coords[i - 1]
        x2, y2 = np_coords[i + 1] - np_coords[i - 1]
        assert x1 * y2 - x2 * y1 == approx(0)

    # tous les segments ont la même longueur ~2.82
    ds_edges, _, _ = _densify_edges(
        None, [edge], max_distance=1, n=None, keep_border=False
    )
    ds_coords = ds_edges[0].coords

    # pour max_distance=1 => 3 nouveaux segments par ancien
    assert len(ds_coords) == 34
    dists = _point_distance(ds_coords)
    for i in range(len(dists)):
        assert dists[i] < 1
