import unittest
from astar import astar, reconstruct_path, euclidean_distance
from draw_astar import draw_exploration

class TestAStar(unittest.TestCase):
    def setUp(self):
        self.graph = {
            (0, 0): {(0, 1): 1, (1, 0): 1},
            (0, 1): {(0, 0): 1, (0, 2): 1},
            (0, 2): {(0, 1): 1, (1, 2): 1},
            (1, 0): {(1, 1): 1, (0, 0): 1},
            (1, 1): {(1, 0): 1, (1, 2): 1, (2, 1): 1},
            (1, 2): {(1, 1): 1, (0, 2): 1, (2, 2): 1},
            (2, 1): {(2, 0): 1, (1, 1): 1, (2, 2): 1},
            (2, 0): {(2, 1): 1},
            (2, 2): {(1, 2): 1, (2, 1): 1}
        }

    def test_astar_shortest_path(self):
        start, end = (0, 0), (2, 2)
        came_from, cost_so_far = astar(self.graph, start, end)
        path = reconstruct_path(came_from, start, end)
        self.assertEqual(path, [(0, 0), (1, 0), (1, 1), (1, 2), (2, 2)])

    def test_astar_no_path(self):
        start, end = (0, 0), (3, 3)
        came_from, cost_so_far = astar(self.graph, start, end)
        self.assertNotIn(end, came_from)

    def test_astar_draw(self):
        start, end = (0, 0), (2, 2)
        came_from, cost_so_far = astar(self.graph, start, end)
        #path = reconstruct_path(came_from, start, end)
        draw_exploration(start, end, came_from)

    def test_euclidean_distance(self):
        self.assertAlmostEqual(euclidean_distance((0, 0), (3, 4)), 5.0)
        self.assertAlmostEqual(euclidean_distance((1, 1), (1, 1)), 0.0)

if __name__ == "__main__":
    unittest.main()
