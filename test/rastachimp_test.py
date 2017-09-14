import numpy as np
from rasterio import features
from rastachimp import simplify
import unittest

class TestSimplify(unittest.TestCase):
    def test_basic(self):
        image = np.array(
              [[2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1],
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
               [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]], 
                            dtype=np.uint8)

        shapes = features.shapes(image)
        simpl = list(simplify(shapes, 2.0))
        
        print(np.dtype(image.dtype).kind)
        
        self.assertEqual(6, len(simpl))
        
        count = {2:0, 1:0, 5:0}
        for _,value in simpl:
            count[value] += 1
        self.assertEqual(count[2],2)
        self.assertEqual(count[1],3)
        self.assertEqual(count[5],1)

if __name__ == '__main__':
    unittest.main()