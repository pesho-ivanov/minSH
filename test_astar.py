import unittest
import os

from astar import astar, reconstruct_path, euclidean_distance
import utils

class TestFastaFunctions(unittest.TestCase):
    def setUp(self):
        self.sequence_length = 100
        self.description = "Test sequence"
        self.file_name = "test_sequence.fasta"
        self.random_sequence = utils.generate_random_sequence(self.sequence_length)
        utils.save_fasta_file(self.file_name, self.description, self.random_sequence)

    def tearDown(self):
        if os.path.exists(self.file_name):
            os.remove(self.file_name)

    def test_read_fasta_file(self):
        read_sequence = utils.read_fasta_file(self.file_name)
        self.assertEqual(self.random_sequence, read_sequence)

    def test_save_fasta_file(self):
        with open(self.file_name, "r") as fasta_file:
            description_line = fasta_file.readline().strip()
            self.assertEqual(f"> {self.description}", description_line)

            sequence = ""
            for line in fasta_file:
                sequence += line.strip()
            self.assertEqual(self.random_sequence, sequence)

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
        came_from, cost_so_far = astar(start, end)
        path = reconstruct_path(came_from, start, end)
        self.assertEqual(path, [(0, 0), (1, 0), (1, 1), (1, 2), (2, 2)])

    def test_astar_no_path(self):
        start, end = (0, 0), (3, 3)
        came_from, cost_so_far = astar(start, end)
        self.assertNotIn(end, came_from)

    def test_astar_draw(self):
        start, end = (0, 0), (2, 2)
        came_from, cost_so_far = astar(start, end)
        #path = reconstruct_path(came_from, start, end)
        utils.draw_exploration(start, end, came_from)

    def test_euclidean_distance(self):
        self.assertAlmostEqual(euclidean_distance((0, 0), (3, 4)), 5.0)
        self.assertAlmostEqual(euclidean_distance((1, 1), (1, 1)), 0.0)

if __name__ == "__main__":
    unittest.main()
