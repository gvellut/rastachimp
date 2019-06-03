import numpy as np
from rasterio import features

from rastachimp import as_shapely, simplify_dp


def test_basic():
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
    simpl = list(simplify_dp(as_shapely(shapes), 2.0))

    print(np.dtype(image.dtype).kind)

    assert 6 == len(simpl)

    count = {2: 0, 1: 0, 5: 0}
    for _, value in simpl:
        count[value] += 1
    assert count[2] == 2
    assert count[1] == 3
    assert count[5] == 1
